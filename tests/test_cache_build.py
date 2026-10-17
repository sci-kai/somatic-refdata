"""Download cache (local HTTP server), manifest of outputs, interval merging."""

import functools
import hashlib
import http.server
import threading

import pytest

from seqc2_toolkit import build, net
from seqc2_toolkit.cache import Cache


@pytest.fixture()
def server(tmp_path):
    """A local HTTP server for the folder tmp_path/www (supports Range requests? no)."""
    www = tmp_path / "www"
    www.mkdir()
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(www))
    handler.log_message = lambda *a: None
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield www, f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def md5(data):
    return hashlib.md5(data).hexdigest()


def test_cache_download_verify_and_reuse(tmp_path, server, monkeypatch):
    www, url = server
    (www / "a.vcf.gz").write_bytes(b"content")
    cache = Cache(tmp_path / "cache")
    path = cache.fetch(f"{url}/a.vcf.gz", md5(b"content"))
    assert path.read_bytes() == b"content"
    assert cache.status(f"{url}/a.vcf.gz", md5(b"content")) == "ok"
    # a second fetch must not download (network disabled)
    monkeypatch.setattr(net, "get_file", lambda *a, **k: pytest.fail("downloaded again"))
    assert cache.fetch(f"{url}/a.vcf.gz", md5(b"content")) == path


def test_cache_redownloads_corrupt_file(tmp_path, server):
    www, url = server
    (www / "a.txt").write_bytes(b"good")
    cache = Cache(tmp_path / "cache")
    path = cache.fetch(f"{url}/a.txt", md5(b"good"))
    path.write_bytes(b"bad!")
    assert cache.status(f"{url}/a.txt", md5(b"good")) == "missing"
    assert cache.fetch(f"{url}/a.txt", md5(b"good")).read_bytes() == b"good"


def test_cache_md5_mismatch_after_download(tmp_path, server):
    www, url = server
    (www / "a.txt").write_bytes(b"other")
    cache = Cache(tmp_path / "cache")
    with pytest.raises(RuntimeError, match="md5 mismatch"):
        cache.fetch(f"{url}/a.txt", md5(b"expected"))
    assert not cache.path(f"{url}/a.txt").exists()


def test_manifest_skip_rebuild_prune(tmp_path):
    root = tmp_path / "ds"
    root.mkdir()
    manifest = build.Manifest(root)
    (root / "x.txt").write_text("x")
    manifest.record("x.txt", {"source": "s1"})
    assert manifest.up_to_date("x.txt", {"source": "s1"}, verify=False)
    # another recipe or a changed file: rebuild
    assert not manifest.up_to_date("x.txt", {"source": "s2"}, verify=False)
    (root / "x.txt").write_text("changed")
    assert not manifest.up_to_date("x.txt", {"source": "s1"}, verify=False)
    # outputs not produced in this run are removed
    fresh = build.Manifest(root)
    fresh.prune()
    assert not (root / "x.txt").exists() and fresh.entries == {}


def test_merge_intervals():
    intervals = [(0, 0, 10, "chr1"), (0, 5, 20, "chr1"), (0, 20, 30, "chr1"), (1, 0, 5, "chr2")]
    assert build.merge(intervals) == [["chr1", 0, 30], ["chr2", 0, 5]]
