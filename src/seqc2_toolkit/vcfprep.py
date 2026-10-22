"""Optional changes of the somatic VCFs of a local dataset (dataset config section 'vcf').

  vcf:
    sample_order: null    sample columns in this order, e.g. [NORMAL, TUMOR] (default: as in
                          the SEQC2 files, TUMOR then NORMAL)
    format_dp: false      add FORMAT/DP = sum of FORMAT/AD per sample to records without
                          FORMAT/DP (the MuTect2 3.8 records of SEQC2 have none)

Both are for workflows that expect these properties (e.g. filters that index the samples
by position, read-depth fields). Records that already have FORMAT/DP are kept unchanged; a
sample with a missing AD value gets DP '.'. The header gets a FORMAT/DP definition (if
missing) and a line ##seqc2_toolkit_vcf=... that documents the changes.

bcftools was considered: 'bcftools view -s' reorders samples, but '+fill-tags' would
overwrite an existing FORMAT/DP and has no per-record "only if missing" rule.
"""

import gzip
import logging

from . import __version__

log = logging.getLogger("seqc2")

DP_HEADER = (
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth: sum of FORMAT/AD '
    '(added by seqc2-toolkit)">'
)


def transform_lines(lines, sample_order=None, format_dp=False, source="<vcf>"):
    """(output lines, stats) of VCF lines (without newline) with the changes applied.

    stats: records, dp_added, dp_existing, ad_missing.

    >>> out, st = transform_lines(["##fileformat=VCFv4.2",
    ...     "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\tT\\tN",
    ...     "1\\t5\\t.\\tA\\tG\\t.\\tPASS\\t.\\tGT:AD\\t0/1:7,3\\t0/0:9"], ["N", "T"], True)
    >>> out[-1].split("\\t")[8:]
    ['GT:AD:DP', '0/0:9:9', '0/1:7,3:10']
    >>> st["dp_added"], st["records"]
    (1, 1)
    """
    out = []
    stats = {"records": 0, "dp_added": 0, "dp_existing": 0, "ad_missing": 0}
    order = None
    has_dp_header = False
    for line in lines:
        if line.startswith("##"):
            # remember whether the header already defines FORMAT/DP
            has_dp_header = has_dp_header or line.startswith("##FORMAT=<ID=DP,")
            out.append(line)
            continue
        if line.startswith("#CHROM"):
            cols = line.split("\t")
            samples = cols[9:]
            wanted = list(sample_order) if sample_order else samples
            # the requested order must be a permutation of the samples in the file
            if sorted(samples) != sorted(wanted):
                raise RuntimeError(
                    f"{source}: samples {samples} don't match vcf.sample_order {wanted}; "
                    "fix vcf.sample_order in the dataset config"
                )
            order = [samples.index(name) for name in wanted]
            if format_dp and not has_dp_header:
                out.append(DP_HEADER)
            changes = []
            if sample_order:
                changes.append(f"sample order {' '.join(wanted)}")
            if format_dp:
                changes.append("FORMAT/DP set to sum(FORMAT/AD) where missing")
            # unstructured header line (no <...>), so that htslib doesn't parse it
            out.append(f"##seqc2_toolkit_vcf=seqc2-toolkit {__version__}; {'; '.join(changes)}")
            out.append("\t".join(cols[:9] + [samples[i] for i in order]))
            continue
        if order is None:
            raise RuntimeError(f"{source}: record before the #CHROM header line")
        cols = line.split("\t")
        stats["records"] += 1
        keys = cols[8].split(":")
        # sample columns in the requested order
        values = [cols[9 + i].split(":") for i in order]
        if format_dp and "DP" in keys:
            # keep an existing FORMAT/DP untouched
            stats["dp_existing"] += 1
        elif format_dp:
            ad_index = keys.index("AD") if "AD" in keys else None
            if ad_index is None:
                stats["ad_missing"] += 1
            for sample_values in values:
                # trailing FORMAT fields may be omitted in a sample column: pad with '.'
                sample_values += ["."] * (len(keys) - len(sample_values))
                ad = sample_values[ad_index].split(",") if ad_index is not None else ["."]
                # DP is the sum of all allele depths; '.' if any depth is missing
                dp = str(sum(int(x) for x in ad)) if all(x.isdigit() for x in ad) else "."
                sample_values.append(dp)
            keys.append("DP")
            stats["dp_added"] += 1
        out.append("\t".join(cols[:8] + [":".join(keys)] + [":".join(v) for v in values]))
    if order is None:
        raise RuntimeError(f"{source}: no #CHROM header line found")
    return out, stats


def transform_file(src, plain_dest, sample_order=None, format_dp=False):
    """Write the changed VCF (bgzipped or plain src) as plain text to plain_dest; stats."""
    with gzip.open(src, "rt") as handle:
        lines, stats = transform_lines(
            (line.rstrip("\n") for line in handle), sample_order, format_dp, str(src)
        )
    plain_dest.write_text("\n".join(lines) + "\n")
    log.info(
        "%s: %d records; FORMAT/DP added to %d, already present in %d, AD missing in %d",
        plain_dest.name,
        stats["records"],
        stats["dp_added"],
        stats["dp_existing"],
        stats["ad_missing"],
    )
    return stats
