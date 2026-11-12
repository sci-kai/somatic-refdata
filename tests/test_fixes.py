"""Regressions of the pre-release review (v0.1.0): one test per fixed defect."""

import http.server
import json
import os
import threading

import pytest

from somatic_refdata import build, catalog, cli, dataset, net, profiles
from somatic_refdata.cache import Cache
from somatic_refdata.dataset import ConfigError

CONTENT = bytes(range(256)) * 4  # 1024 bytes


class FlakyHandler(http.server.BaseHTTPRequestHandler):
    """Serves CONTENT; the first full request announces all bytes but sends only 400
    (connection closed early); Range requests are answered correctly."""

    calls = []

    def do_GET(self):  # noqa: N802 (http.server API)
        self.calls.append(self.headers.get("Range"))
        rng = self.headers.get("Range")
        if rng:
            start = int(rng.split("=")[1].split("-")[0])
            if start >= len(CONTENT):
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{len(CONTENT)}")
                self.end_headers()
                return
            body = CONTENT[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(CONTENT) - 1}/{len(CONTENT)}")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(CONTENT)))
        self.end_headers()
        self.wfile.write(CONTENT[:400])
        self.close_connection = True

    def log_message(self, *args):
        pass


@pytest.fixture
def server():
    FlakyHandler.calls = []
    httpd = http.server.HTTPServer(("127.0.0.1", 0), FlakyHandler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}/file.bin"
    httpd.shutdown()


def test_truncated_download_is_resumed(server, tmp_path):
    part = tmp_path / "file.bin.part"
    net.get_file(server, part, retries=3, delay=0)
    assert part.read_bytes() == CONTENT
    assert FlakyHandler.calls == [None, "bytes=400-"]


def test_partial_file_of_wrong_size_is_not_accepted_on_416(server, tmp_path):
    part = tmp_path / "file.bin.part"
    part.write_bytes(b"x" * 2000)  # larger than the file: 416, then a fresh download
    net.get_file(server, part, retries=4, delay=0)
    assert part.read_bytes() == CONTENT


def test_output_does_not_write_through_links_into_the_cache(tmp_path):
    cached = tmp_path / "cache" / "big.bam"
    cached.parent.mkdir()
    cached.write_bytes(b"whole BAM")
    config = {"subset": {"contigs": None}, "data": {}}
    builder = build.Builder(config, {"base_url": "https://x"}, None, tmp_path / "ds")
    builder.output("bam/a.bam", {"v": 1}, lambda path: build.link(cached, path))
    assert (tmp_path / "ds/bam/a.bam").samefile(cached)
    # another recipe (e.g. BAM windows): the slice replaces the link, not the cached file
    builder.output("bam/a.bam", {"v": 2}, lambda path: path.write_bytes(b"slice"))
    assert cached.read_bytes() == b"whole BAM"
    assert (tmp_path / "ds/bam/a.bam").read_bytes() == b"slice"


def test_manifest_written_atomically(tmp_path):
    manifest = build.Manifest(tmp_path)
    manifest.save()
    assert json.loads((tmp_path / build.MANIFEST).read_text())["files"] == {}
    assert not list(tmp_path.glob("*.tmp"))


def test_run_command_makes_a_local_pipeline_absolute(tmp_path, monkeypatch):
    (tmp_path / "pipe").mkdir()
    monkeypatch.chdir(tmp_path)
    profile = {**profiles.PROFILE_DEFAULTS, "name": "p", "revision": "v1"}
    cmd = profiles.run_command(profile, "params.yaml", "out", "pipe", None, None, [], [])
    assert cmd[2] == str(tmp_path / "pipe")
    assert "-r" not in cmd  # a local folder has no revision


def test_whole_bams_with_a_subset_are_refused():
    config = {"subset": {"contigs": ["chr22"], "bam_subsample": None, "bam_windows": None}}
    with pytest.raises(ConfigError, match="subset.contigs with data.bam but no"):
        build.check_bam_subset(config, [{"bams": {"T": "x"}}])
    build.check_bam_subset(config, [{"bams": {}}])
    build.check_bam_subset({"subset": {**config["subset"], "bam_windows": ["chr22:1-9"]}}, [])


def test_profile_index_prefers_the_recorded_index():
    unit = {"bams": {"T": "https://s/WGS_IL_T_1.bwa.dedup.bam"}, "bais": {}}
    assert profiles.index_of(unit, "T").endswith(".bam.bai")
    unit["bais"] = {"T": "https://s/WGS_IL_T_1.bwa.dedup.bai"}
    assert profiles.index_of(unit, "T") == "https://s/WGS_IL_T_1.bwa.dedup.bai"
    assert profiles.index_of(unit, "N") is None
    assert profiles.index_of({"bams": {"T": "bam/x.cram"}}, "T") == "bam/x.cram.crai"


def write_config(tmp_path, text):
    path = tmp_path / "ds.yaml"
    path.write_text(text)
    return path


def test_selection_options_of_seqc2_on_hg008_are_an_error(tmp_path):
    config = dataset.load(write_config(tmp_path, "source: hg008\n"))
    with pytest.raises(ConfigError, match="--pairs: not available for source hg008"):
        dataset.override_select(config, pairs="T-ND")


def test_defaults_are_not_changed_by_overrides(tmp_path):
    before = json.dumps(dataset.DEFAULTS, sort_keys=True)
    config = dataset.load(write_config(tmp_path, "{}\n"))
    dataset.override_select(config, technologies="WES", centers="all")
    assert json.dumps(dataset.DEFAULTS, sort_keys=True) == before


def test_md5_list_rejects_error_pages(tmp_path, monkeypatch):
    cache = Cache(tmp_path / "cache")
    url = "https://example.org/list.md5"
    page = '<?xml version="1.0"?><Error>Service unavailable</Error>'
    monkeypatch.setattr(net, "get_text", lambda u: page)
    with pytest.raises(RuntimeError, match="can't be fetched"):
        catalog.md5_list(cache, url)
    assert not cache.path(url).exists()
    # a cached list stays in use when the refresh gets an error page
    cache.path(url).parent.mkdir(parents=True)
    cache.path(url).write_text(f"{'a' * 32}  x.vcf.gz\n")
    assert catalog.md5_list(cache, url, refresh=True) == [f"{'a' * 32}  x.vcf.gz"]


def hg008_facts(config):
    """dataset.json of a local HG008 dataset built with config."""
    return {
        "name": "hg",
        "source": "hg008",
        "source_name": "HG008",
        "access": "local",
        "resources": {"reference": "reference/ref.fa.gz"},
        "vcf_samples": {"tumor": None, "normal": None},
        "config": {k: config[k] for k in ("source", "select", "data", "subset")},
        "units": [
            {
                "pair": "T-ND",
                "center": "NYGC",
                "replicate": None,
                "technology": "illumina-nygc",
                "series": "illumina-nygc",
                "platform": "Illumina NovaSeq 6000, PCR-free 2x150",
                "aligner": "bwa-mem",
                "build": "GRCh38",
                "library_type": "wgs",
                "samples": {"T": "HG008-T", "N": "HG008-N-D"},
                "coverage": {"T": 161, "N": 118},
                "target_bed": None,
                "vcfs": [
                    {
                        "aligner": "bwa-mem",
                        "caller": "muTect2",
                        "path": "vcf/T-ND_illumina-nygc.bwa-mem.muTect2.vcf.gz",
                        "tumor_id": "HG008-T",
                        "normal_id": "HG008-N",
                    }
                ],
                "bams": {},
            }
        ],
    }


def test_samplesheet_of_a_local_hg008_dataset(tmp_path, capsys):
    path = write_config(tmp_path, "source: hg008\n")
    outdir = tmp_path / "out"
    outdir.mkdir()
    config = dataset.load(path)
    (outdir / "dataset.json").write_text(json.dumps(hg008_facts(config)))
    argv = ["samplesheet", str(path), "--outdir", str(outdir), "--profile"]
    cli.main([*argv, "variantinterpretation"])
    sheet = (outdir / "variantinterpretation" / "samplesheet.csv").read_text()
    assert "HG008_T-ND_illumina-nygc_bwa-mem_muTect2" in sheet
    # a dataset built with another selection is refused
    other = write_config(tmp_path, "source: hg008\nselect: {callers: []}\ndata: {vcf: false}\n")
    with pytest.raises(SystemExit):
        cli.main(["samplesheet", str(other), "--outdir", str(outdir), "--profile", "variantcomp"])
    assert "built with another selection" in capsys.readouterr().err


def test_clear_errors_instead_of_tracebacks(tmp_path, capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["plan", str(tmp_path / "missing.yaml")])
    assert exit_info.value.code == 1
    assert "ERROR: FileNotFoundError" in capsys.readouterr().err


def test_verify_without_a_built_dataset_fails(tmp_path, capsys):
    path = write_config(tmp_path, "{}\n")
    with pytest.raises(SystemExit) as exit_info:
        cli.main(["verify", str(path), "--outdir", str(tmp_path / "never")])
    assert exit_info.value.code == 1
    assert "run 'somref fetch" in capsys.readouterr().err


def test_default_cache_folder(monkeypatch):
    from somatic_refdata.cache import default_cache_dir

    monkeypatch.delenv("SOMREF_CACHE", raising=False)
    assert str(default_cache_dir()).endswith(os.path.join(".cache", "somref"))
