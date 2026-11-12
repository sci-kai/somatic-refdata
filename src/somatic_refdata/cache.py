"""Download cache shared by all datasets: every remote file is downloaded once.

Layout: <cache>/<host>/<path of the URL>, e.g.
  ~/.cache/somref/ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/.../WGS_IL_1.bwa.muTect2.vcf.gz

<cache>/registry.json records for every cached file its size, modification time and md5,
and whether the md5 was checked against the server's md5 list. A file whose size and
modification time are unchanged is not hashed again (unless verify=True).
"""

import hashlib
import json
import logging
import os
from datetime import datetime
from pathlib import Path
from urllib.parse import urlparse

from . import net

log = logging.getLogger("somref")


def default_cache_dir():
    """Cache folder: $SOMREF_CACHE, else ~/.cache/somref."""
    return Path(os.environ.get("SOMREF_CACHE") or Path.home() / ".cache" / "seqc2")


def md5sum(path, chunk=1 << 20):
    """md5 hex digest of a file."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            digest.update(block)
    return digest.hexdigest()


def file_state(path):
    """(size, modification time in ns) of a file: detects changes without hashing."""
    stat = Path(path).stat()
    return stat.st_size, stat.st_mtime_ns


class Cache:
    """Downloads remote files once, verifies them against known md5 sums."""

    def __init__(self, root):
        self.root = Path(root).expanduser().resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.registry_path = self.root / "registry.json"
        self.registry = {}
        if self.registry_path.exists():
            self.registry = json.loads(self.registry_path.read_text())

    def path(self, url):
        """Local path of a URL in the cache."""
        parsed = urlparse(url)
        return self.root / parsed.netloc / parsed.path.lstrip("/")

    def save(self):
        """Write the registry (after every change, so an interrupted run keeps it)."""
        tmp = self.registry_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.registry, indent=1, sort_keys=True))
        tmp.replace(self.registry_path)

    def _record(self, key, path, md5, verified):
        """Remember size, modification time and md5 of a cached file."""
        size, mtime = file_state(path)
        self.registry[key] = {
            "size": size,
            "mtime_ns": mtime,
            "md5": md5,
            "md5_verified": verified,
            "checked": datetime.now().isoformat(timespec="seconds"),
        }
        self.save()

    def status(self, url, md5=None):
        """'ok' (cached, md5 fine), 'unverified' (cached, no md5 known) or 'missing'.

        Uses the registry only (no hashing): for plans and summaries.
        """
        path = self.path(url)
        entry = self.registry.get(str(path.relative_to(self.root)))
        if not path.exists() or entry is None:
            return "missing"
        if (entry["size"], entry["mtime_ns"]) != file_state(path):
            return "missing"
        if md5 is not None and entry["md5"] != md5:
            return "missing"
        return "ok" if entry["md5_verified"] else "unverified"

    def fetch(self, url, md5=None, verify=False):
        """Local path of url: downloaded if missing or corrupt, md5-checked if md5 is known.

        md5: expected md5 (from the server's md5 lists) or None (file is only recorded).
        verify: hash cached files again even if size and modification time are unchanged.
        """
        path = self.path(url)
        key = str(path.relative_to(self.root))
        entry = self.registry.get(key)
        if path.exists():
            # unchanged since it was checked: no need to hash a large file again
            unchanged = entry is not None and (entry["size"], entry["mtime_ns"]) == file_state(path)
            if unchanged and not verify and (md5 is None or entry["md5"] == md5):
                log.debug("cached: %s", key)
                return path
            # changed, unknown or verification requested: compare with the expected md5
            observed = md5sum(path)
            if md5 is None or observed == md5:
                self._record(key, path, observed, md5 is not None)
                log.info("cached (md5 %s): %s", "ok" if md5 else "recorded", path.name)
                return path
            log.warning(
                "md5 mismatch for cached %s (%s != %s): downloading again", key, observed, md5
            )
            path.unlink()
        # download into <file>.part (resumable), check, then move into place
        path.parent.mkdir(parents=True, exist_ok=True)
        part = path.with_name(path.name + ".part")
        log.info("download: %s", url)
        net.get_file(url, part)
        observed = md5sum(part)
        if md5 is not None and observed != md5:
            part.unlink()
            raise RuntimeError(
                f"md5 mismatch after download of {url}: expected {md5}, got {observed}; "
                "run the command again (the server may have delivered a broken file)"
            )
        part.rename(path)
        self._record(key, path, observed, md5 is not None)
        return path
