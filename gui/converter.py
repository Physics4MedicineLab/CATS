"""
Convert between CATS CSV files and BED files
while preserving metadata and comments.
"""
import os
import sys
from time import (
    localtime,
    strftime
)
from io import StringIO
from typing import List, Tuple, Dict, Set, Union

import pandas as pd

# The nine mandatory BED columns in their canonical order.
# Do NOT change this list unless you know exactly what you are doing.
MANDATORY_BED: List[str] = [
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRgb",
]

def _dedup(seq: List[str]) -> List[str]:
    """Return `seq` with duplicates removed while preserving order.

    The first occurrence of each element is kept; subsequent duplicates
    are discarded.
    """
    seen: Set[str] = set()
    out: List[str] = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def read_bed(path: str) -> Tuple[List[str], List[str], List[List[str]]]:
    """Read a BED file and return `(comments, header, rows)`.

    Parameters
    ----------
    path:
        Path to the BED file.

    Returns
    -------
    tuple
        * `comments`: Every line starting with `# ` (without the trailing
          newline).
        * `header`:   Column names taken from the first comment line that
          starts with `# ` and contains at least one tab. Duplicates are
          removed while preserving order.
        * `rows`:     The data rows as lists of strings.
    """
    comments: List[str] = []
    rows: List[List[str]] = []
    header: List[str] | None = None

    with open(path, newline="") as fh:
        for ln in fh:
            ln_clean = ln.rstrip("\r\n")  # Strip newline characters

            # Comment line
            if ln_clean.startswith("#"):
                comments.append(ln_clean)

                # The first comment line that contains a TAB defines the header
                if header is None and ln_clean.startswith("# ") and "\t" in ln_clean:
                    header = _dedup([col.strip() for col in ln_clean[2:].split("\t")])
                continue

            # Data line
            parts = ln_clean.split("\t")

            # If no header found yet, create a generic one based on width (it should not happen)
            if header is None:
                header = [f"col_{i + 1}" for i in range(len(parts))]

            # If the row is wider than the header, fold overflow into last col
            if len(parts) > len(header):
                parts = parts[: len(header) - 1] + ["\t".join(parts[len(header) - 1 :])]

            # Pad short rows so that every row matches header length
            # (it should not happen, since len(header) - len(parts) should be 0)
            parts.extend(["."] * (len(header) - len(parts)))
            rows.append(parts)

    if header is None:
        raise ValueError("Cannot determine BED header – no suitable '# <header>' line found.")

    header = _dedup(header)
    return comments, header, rows


def read_csv(path: str) -> Tuple[List[str], pd.DataFrame]:
    """
    Read a CSV file while preserving leading comment lines.
    """
    comments: List[str] = []
    data_lines: List[str] = []

    with open(path, newline="") as fh:
        for ln in fh:
            if ln.startswith("#"):
                comments.append(ln)
            else:
                data_lines.append(ln)

    if not data_lines:
        raise ValueError("CSV has no data rows after comments.")

    df = pd.read_csv(StringIO("".join(data_lines)), dtype=str, keep_default_na=False)

    # Drop duplicate columns if needed
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]

    return comments, df


def _split_index(idx: str) -> Tuple[str, int]:
    """Split an index string of the form 'chrom:position' into its parts."""
    chrom, pos = idx.split(":")
    return chrom, int(pos)


def _pam_row(
    *,
    idx: str,
    seq: str,
    transcript_id: str,
    gene_name: str,
    strand: Union[str, None],
    extra: Dict[str, str],
    prefix: str,
    color: str,
) -> Dict[str, str]:
    """Create a dictionary representing a single BED row for a PAM site."""
    chrom, pos = _split_index(idx)
    start = pos - 1
    end = start + (len(seq) if seq else 1)

    return {
        "chrom": chrom,
        "chromStart": start,
        "chromEnd": end,
        "name": f"{prefix}{transcript_id} | {gene_name} | {extra.get('Biotype') or extra.get('Regions', '.')}",
        "score": "0",
        "strand": strand or ".",
        "thickStart": start,
        "thickEnd": end,
        "itemRgb": color,
        **extra,
    }

def csv_to_bed(df: pd.DataFrame) -> Tuple[List[str], pd.DataFrame]:
    """Convert a CATS CSV *DataFrame* to a BED-like *DataFrame*.

    The returned BED frame includes the mandatory nine columns followed by
    any extra columns that were present in `df`.
    """
    extra_cols: List[str] = [col for col in df.columns if col not in MANDATORY_BED]
    header_cols: List[str] = MANDATORY_BED + extra_cols  # Already unique

    rows: List[Dict[str, str]] = []
    for _, row in df.iterrows():
        extra: Dict[str, str] = {col: (row[col] or ".") for col in extra_cols}

        # Variant position
        if row.get("Variant position"):
            chrom, span = row["Variant position"].split(":")
            start_s, end_s = span.split("-")
            start, end = int(start_s), int(end_s)
            rows.append(
                {
                    "chrom": chrom,
                    "chromStart": start - 1,
                    "chromEnd": end,
                    "name": f"Variant:{row.get('Variant name', '.')}",
                    "score": "0",
                    "strand": row.get("Strand", "."),
                    "thickStart": start - 1,
                    "thickEnd": end,
                    "itemRgb": "255,0,0",
                    **extra,
                }
            )

        # PAM positions
        if row.get("Matched seq index"):
            rows.append(
                _pam_row(
                    idx=row["Matched seq index"],
                    seq=row.get("Matched seq", ""),
                    transcript_id=row.get("Transcript ID", "."),
                    gene_name=row.get("Gene Name", "."),
                    strand=row.get("Strand", "."),
                    extra=extra,
                    prefix="PAM:",
                    color="0,0,255",
                )
            )
        elif row.get("First seq index") and row.get("Second seq index"):
            rows.extend(
                [
                    _pam_row(
                        idx=row["First seq index"],
                        seq=row.get("First seq", ""),
                        transcript_id=row.get("Transcript ID", "."),
                        gene_name=row.get("Gene Name", "."),
                        strand=row.get("Strand", "."),
                        extra=extra,
                        prefix="PAM 1:",
                        color="0,0,255",
                    ),
                    _pam_row(
                        idx=row["Second seq index"],
                        seq=row.get("Second seq", ""),
                        transcript_id=row.get("Transcript ID", "."),
                        gene_name=row.get("Gene Name", "."),
                        strand=row.get("Strand", "."),
                        extra=extra,
                        prefix="PAM 2:",
                        color="0,255,0",
                    ),
                ]
            )

    if not rows:
        raise ValueError("No convertible rows found in CSV – check column names.")

    bed_df = pd.DataFrame(rows, columns=header_cols).drop_duplicates()
    return header_cols, bed_df

def bed_to_csv(header: List[str], rows: List[List[str]]) -> pd.DataFrame:
    """Extract the non-mandatory columns from a BED file into a CSV *DataFrame*."""
    header = _dedup(header)
    if len(header) <= len(MANDATORY_BED):
        raise ValueError(
            "BED file has no columns beyond the mandatory nine – nothing to write to CSV."
        )

    df_full = pd.DataFrame(rows, columns=header)

    # Drop duplicate columns if any
    if df_full.columns.duplicated().any():
        df_full = df_full.loc[:, ~df_full.columns.duplicated()]

    extra_cols = [col for col in df_full.columns if col not in MANDATORY_BED]
    return df_full[extra_cols]

def _write_comments(fh, comments: List[str], to_csv: bool = False) -> None:
    """Write comment lines to `fh`.

    If `to_csv` is True, the last comment line (usually the BED header) is
    suppressed because the CSV format has its own header.
    """
    comments_to_write = comments[:-1] if to_csv else comments
    end_char = "\n" if to_csv else ""
    for ln in comments_to_write:
        fh.write(ln + end_char)

def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("Usage: python converter.py <file.(csv|bed)>")

    in_path = os.path.abspath(sys.argv[1])
    if not os.path.isfile(in_path):
        sys.exit(f"File not found: {in_path}")

    stem, ext = os.path.splitext(in_path)
    ext = ext.lower()
    if ext not in (".csv", ".bed"):
        sys.exit("Unsupported extension – use .csv or .bed")

    out_path = f"{stem}{'.bed' if ext == '.csv' else '.csv'}"

    if ext == ".csv": # .csv -> .bed
        comments, df_csv = read_csv(in_path)
        header_cols, df_bed = csv_to_bed(df_csv)

        with open(out_path, "w") as fh:
            _write_comments(fh, comments)
            fh.write("# " + "\t".join(header_cols) + "\n")
            df_bed.to_csv(fh, sep="\t", header=False, index=False, na_rep=".")

    else:  # .bed -> .csv
        comments, header, rows = read_bed(in_path)
        df_csv = bed_to_csv(header, rows)

        with open(out_path, "w") as fh:
            _write_comments(fh, comments, to_csv=True)
        df_csv.drop_duplicates(inplace=True)
        df_csv.to_csv(out_path, mode="a", index=False)

    print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tCONVERTER - File converted in: {out_path}")


if __name__ == "__main__":
    main()
