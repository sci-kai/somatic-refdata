"""HTTPS access: CA bundle, resumable downloads, small text files, remote file facts.

Downloads use Python (urllib), not curl: conda builds of curl/htslib sometimes keep an
unrelocated build-time CA path (`curl: (77) error adding trust anchors from file:
/home/conda/feedstock_root/...`). The CA bundle found here is also exported as
CURL_CA_BUNDLE / SSL_CERT_FILE for samtools/bcftools reading remote files.
"""

import logging
import os
import shutil
import ssl
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

log = logging.getLogger("somref")

# CA bundles of common Linux distributions (Debian/Ubuntu, RHEL/Fedora, openSUSE, Alpine)
SYSTEM_CA_BUNDLES = [
    "/etc/ssl/certs/ca-certificates.crt",
    "/etc/pki/tls/certs/ca-bundle.crt",
    "/etc/ssl/ca-bundle.pem",
    "/etc/ssl/cert.pem",
]
NETWORK_ERRORS = (urllib.error.URLError, TimeoutError, ConnectionError, ssl.SSLError)

_context = None  # SSL context, set by setup_tls()


def find_ca_bundle(explicit=None):
    """First existing CA bundle: explicit, env vars, certifi, system paths (or None)."""
    candidates = [explicit] if explicit else []
    # environment variables only count if they point to a file (a broken conda
    # placeholder path is skipped)
    candidates += [os.environ.get(v) for v in ("CURL_CA_BUNDLE", "SSL_CERT_FILE")]
    candidates.append(os.environ.get("REQUESTS_CA_BUNDLE"))
    try:
        import certifi

        candidates.append(certifi.where())
    except ImportError:
        pass
    candidates += SYSTEM_CA_BUNDLES
    for path in candidates:
        if path and Path(path).is_file():
            return str(path)
    if explicit:
        sys.exit(f"ERROR: --ca-bundle {explicit} does not exist")
    return None


def setup_tls(explicit=None):
    """Create the SSL context and export the CA bundle for htslib/libcurl."""
    global _context
    bundle = find_ca_bundle(explicit)
    if bundle is None:
        log.warning("no CA bundle found (set --ca-bundle or install certifi); using defaults")
        _context = ssl.create_default_context()
        return None
    # htslib reads CURL_CA_BUNDLE for https:// inputs (samtools faidx/view on URLs)
    os.environ["CURL_CA_BUNDLE"] = bundle
    os.environ["SSL_CERT_FILE"] = bundle
    _context = ssl.create_default_context(cafile=bundle)
    log.debug("CA bundle: %s", bundle)
    return bundle


def context():
    """The SSL context (created with the default CA search on first use)."""
    if _context is None:
        setup_tls()
    return _context


def total_size(headers, status):
    """Full size of the remote file from a response: Content-Range (206, 416) or
    Content-Length (200); None if the server doesn't say."""
    content_range = headers.get("Content-Range") or ""
    if "/" in content_range and content_range.rsplit("/", 1)[1].strip().isdigit():
        return int(content_range.rsplit("/", 1)[1])
    length = headers.get("Content-Length")
    if status == 200 and length and length.isdigit():
        return int(length)
    return None


def get_file(url, part, retries=5, delay=10):
    """Download url into the file part; resumes a partial file (HTTP Range); retries.

    A body shorter than the size the server announced (connection closed early) is an
    incomplete attempt: the next attempt continues at its end. Raises RuntimeError if the
    download fails after all attempts.
    """
    for attempt in range(1, retries + 1):
        # bytes already downloaded by an earlier, interrupted attempt
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, context=context(), timeout=300) as resp:
                # 206: the server continues at offset; 200: it sends the whole file
                mode = "ab" if offset and resp.status == 206 else "wb"
                total = total_size(resp.headers, resp.status)
                with open(part, mode) as fh:
                    shutil.copyfileobj(resp, fh, 1 << 20)
            size = part.stat().st_size
            if total is None or size == total:
                return
            if size > total:
                # more than the file has: a stale partial file of another version
                part.unlink()
                reason = f"{size} bytes, more than the {total} of the file; starting again"
            else:
                reason = f"incomplete: {size} of {total} bytes"
        except urllib.error.HTTPError as err:
            # 416: nothing left to download, if the partial file has the full size
            if err.code == 416 and offset:
                total = total_size(err.headers or {}, 416)
                if total is None or total == offset:
                    return
                part.unlink()
                reason = f"partial file of {offset} bytes, the file has {total}; starting again"
                log.warning("download of %s: %s", url, reason)
                continue
            # other client errors (404, 403, ...) won't go away by retrying
            if 400 <= err.code < 500:
                raise RuntimeError(f"download of {url} failed: HTTP {err.code}") from err
            reason = f"HTTP {err.code}"
        except NETWORK_ERRORS as err:
            reason = str(getattr(err, "reason", err))
        log.warning("download of %s failed (%s); attempt %d/%d", url, reason, attempt, retries)
        if attempt < retries:
            time.sleep(delay)
    raise RuntimeError(f"download of {url} failed after {retries} attempts")


def get_text(url):
    """Content of a small text file at url, or None if it can't be fetched."""
    try:
        with urllib.request.urlopen(url, context=context(), timeout=120) as resp:
            return resp.read().decode()
    except NETWORK_ERRORS as err:
        log.warning("could not fetch %s: %s", url, getattr(err, "reason", err))
        return None


def head(url):
    """Size (bytes) and Last-Modified of the file at url: (size or None, date or None)."""
    try:
        request = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(request, context=context(), timeout=60) as resp:
            size = resp.headers.get("Content-Length")
            return (int(size) if size else None, resp.headers.get("Last-Modified"))
    except NETWORK_ERRORS:
        return (None, None)
