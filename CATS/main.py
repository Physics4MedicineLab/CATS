import argparse
import gzip
import itertools
import os
import re
import tempfile
from time import (
    localtime,
    strftime
)
from typing import Union

import gffutils
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
from parsley import makeGrammar

# Mapping of IUPAC notation to regex pattern
IUPAC_MAP = {
    "W": "[AT]",
    "N": "[ACGT]",
    "*": "[ACGT]",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[CG]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    ".": "[-.]",
    "-": "[-.]"
}

class SequenceVariant:
    """Parsed variant sequence in hgvs format"""
    def __init__(self, ac, gene, var_type, posedit):
        self.ac = ac
        self.gene = gene
        self.var_type = var_type
        self.posedit = posedit

    def __repr__(self):
        return (f"SequenceVariant(ac={self.ac!r}, gene={self.gene!r}, "
                f"var_type={self.var_type!r}, posedit={self.posedit!r})")

class PosEdit:
    """Position/Edit"""
    def __init__(self, pos, edit, uncertain=False):
        self.pos = pos
        self.edit = edit
        self.uncertain = uncertain

    def __repr__(self):
        return (f"PosEdit(pos={self.pos!r}, edit={self.edit!r}, "
                f"uncertain={self.uncertain!r})")

class NARefAlt:
    """Nucleic Acid Reference/Alternative"""
    def __init__(self, ref, alt, mut_type=None):
        self.ref = ref
        self.alt = alt
        self.mut_type = mut_type  # e.g., "substitution", "deletion", etc.

    def __repr__(self):
        return f"NARefAlt(ref={self.ref!r}, alt={self.alt!r}, mut_type={self.mut_type!r})"

grammar = r"""
hgvs_variant = accn:ac opt_gene_expr:gene ':' var_type:vt '.' posedit:posedit -> SequenceVariant(ac, gene, vt, posedit)

# Variant type letter
var_type = <('c'|'g'|'m'|'n'|'p'|'r')>

# Posedit options
posedit = substitution | delins | deletion | insertion | duplication

substitution = pos:pos dna_subst:edit -> PosEdit(pos, edit)
delins       = pos:pos dna_delins:edit -> PosEdit(pos, edit)
deletion     = pos:pos dna_del:edit   -> PosEdit(pos, edit)
insertion    = pos:pos dna_ins:edit   -> PosEdit(pos, edit)
duplication  = pos:pos dna_dup:edit   -> PosEdit(pos, edit)

pos = <digit+>:x -> int(x)

# Nucleic Acid Ref/Alt
dna_subst  = dna:ref '>' dna:alt    -> NARefAlt(ref, alt, "substitution")
dna_delins = 'delins' dna_seq:seq    -> NARefAlt(None, seq, "delins")
dna_del    = 'del' !('ins') (dna_seq:seq | ->None):seq -> NARefAlt(seq, None, "deletion")
dna_ins    = 'ins' dna_seq:seq       -> NARefAlt(None, seq, "insertion")
dna_dup    = 'dup' (dna_seq:seq | ->None):seq -> NARefAlt(seq, seq, "duplication")

dna_seq = <(('A'|'C'|'G'|'T'))+>:x -> x
dna = dna_seq

# Accession
accn = <letter (letterOrDigit | ('-' | '_'))* ('.' digit+)?>

# Gene
opt_gene_expr = (paren_gene | ->None):gene -> gene
paren_gene = '(' gene_symbol:sym ')' -> sym
gene_symbol = <letter (letterOrDigit | ('-' | '_'))+>

# Basic Tokens
digit = :x ?(x in "0123456789")
letter = :x ?(x.isalpha())
letterOrDigit = :x ?(x.isalnum())
"""

def create_sequence_pattern(sequence):
    """
    Create a regex pattern for a given sequence using `IUPAC_MAP`.
    """
    pattern = ""
    for symbol in sequence:
        pattern += IUPAC_MAP.get(symbol, symbol)
    return pattern


def is_gzipped(file_path: str) -> bool:
    """
    Check if a file is gzipped by reading its magic number.
    """
    with open(file_path, 'rb') as f:
        magic_number = f.read(2)
    return magic_number == b'\x1f\x8b'


def get_pathogenic_variants_from_api(snv: bool, verbose=True, gene_list=None):
    """
    Fetch pathogenic variants from ClinVar E-Utils API.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = base_url + "esearch.fcgi"
    summary_url = base_url + "esummary.fcgi"
    db = "clinvar"

    gene_term = ''
    if gene_list:
        gene_term = '(' + ' OR '.join(f'"{gene}"[Gene Name]' for gene in gene_list) + ')'

    if snv:
        term = ('pathogenic AND (("clinsig pathogenic"[Properties] '
                'OR "clinsig pathogenic low penetrance"[Properties] '
                'OR "clinsig established risk allele"[Properties]) '
                'AND ("single nucleotide variant"[Type of variation]))')
    else:
        term = ('pathogenic AND (("clinsig pathogenic"[Properties] '
                'OR "clinsig pathogenic low penetrance"[Properties] '
                'OR "clinsig established risk allele"[Properties]))')

    if gene_term:
        term = f'{term} AND {gene_term}'

    if verbose:
        count_params = {"db": db, "term": term, "retmode": "json", "retmax": 0}
        count_response = requests.get(search_url, params=count_params)
        count_response.raise_for_status()
        total_count = int(count_response.json()['esearchresult']['count'])
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tTotal number of variants: {total_count}.")
        if total_count == 0:
            return {}

    search_params = {"db": db, "term": term, "retmode": "json", "retmax": 999999}
    search_response = requests.get(search_url, params=search_params)
    search_response.raise_for_status()
    variant_ids = search_response.json()['esearchresult']['idlist']
    pathogenic_variants = {}

    for variant_id in variant_ids:
        summary_params = {"db": db, "id": variant_id, "retmode": "json"}
        summary_response = requests.get(summary_url, params=summary_params)
        summary_response.raise_for_status()
        variant_info = summary_response.json()['result'].get(variant_id, {})
        variation_title = variant_info.get('title', '')
        for variant in variant_info.get("variation_set", []):
            variant_type = variant.get("variant_type", "")
            for loc in variant.get('variation_loc', []):
                if loc['assembly_name'] == 'GRCh38':
                    chrom = "chr" + loc['chr']
                    pos_start = int(loc['start'])
                    pos_stop = int(loc['stop'])
                    pathogenic_variants.setdefault(chrom, []).append(
                        (pos_start, pos_stop, variation_title, variant_type)
                    )
    return pathogenic_variants


def extract_transcripts(db):
    """
    Extract transcripts from GFF/GTF feature database.
    """
    transcripts = {}
    for transcript in db.features_of_type('transcript'):
        if "MANE_Select" in transcript.attributes.get('tag', [None]):
            transcripts[transcript.id] = {
                'id': transcript.id,
                'chrom': transcript.chrom,
                'start': transcript.start,
                'end': transcript.end,
                'strand': transcript.strand,
                'gene_id': transcript.attributes['gene_id'][0],
                'gene_name': transcript.attributes.get('gene_name', [None])[0],
                'transcript_type': transcript.attributes.get('transcript_type', [None])[0],
            }
    return transcripts


def process_record_w_transcripts(args):
    """
    Worker function for parallel searches in parse_fasta() when transcripts exist.
    """
    (record, seq1_pattern, seq2_pattern, transcripts,
     pathogenic_variants, window_size, num_bases, 
     gene_list, single_seq) = args

    result = []
    record_transcript_id = record.id.split('|')[0]
    record_info = transcripts.get(record_transcript_id)
    if not record_info:
        return result
    seq = str(record.seq)[::-1] if record_info["strand"] == "-" else str(record.seq)
    seq1_regex = re.compile(seq1_pattern, re.IGNORECASE)

    if gene_list and record_info['gene_name'] not in gene_list:
        return result

    # SINGLE-SEQUENCE SEARCH
    if single_seq:
        for match in seq1_regex.finditer(seq):
            start_idx = max(0, match.start() - num_bases)
            end_idx = min(len(seq), match.start() + len(match.group()) + num_bases)
            transcript_parts = record.id.split("|")
            base_dict = {
                "Transcript ID": f"{transcript_parts[0]}|{transcript_parts[1]}",
                "Gene Name": transcript_parts[5],
                "Biotype": '|'.join(transcript_parts[7:]),
                "Strand": record_info["strand"],
                "Sequence": seq[start_idx:end_idx][::-1] if record_info["strand"] == "-" else str(record.seq),
                "Matched seq": match.group(),
                "Matched seq index": f"{record_info['chrom']}:{record_info['start'] + match.start()}"
            }
            if not pathogenic_variants:
                result.append(base_dict)
            else:
                if record_info["chrom"] in pathogenic_variants:
                    for var_start, var_stop, var_title, _ in pathogenic_variants[record_info["chrom"]]:
                        seq_region_start = record_info['start'] + start_idx
                        seq_region_end = record_info['start'] + end_idx

                        region_start = record_info['start'] + match.start()
                        region_end = region_start + len(match.group())

                        if var_start >= region_start and var_stop <= region_end:
                            var_distance = 0
                        else:
                            var_distance = min([
                                abs(region_start - var_start),
                                abs(region_end - var_start),
                                abs(region_start - var_stop),
                                abs(region_end - var_stop),
                            ]) + 1

                        if (var_start <= seq_region_end) and (var_stop >= seq_region_start):
                            variant_dict = base_dict.copy()
                            variant_dict.update({
                                "Variant position": f"{record_info['chrom']}:{var_start}-{var_stop}",
                                "Variant distance": var_distance,
                                "Variant name": var_title
                            })
                            result.append(variant_dict)
    # DOUBLE-SEQUENCE SEARCH
    else:
        seq2_regex = re.compile(seq2_pattern, re.IGNORECASE)
        matches_seq1 = list(seq1_regex.finditer(seq))
        matches_seq2 = list(seq2_regex.finditer(seq))
        for idx1, idx2 in itertools.product(matches_seq1, matches_seq2):
            if abs(idx1.start() - idx2.start()) <= window_size:
                start_idx = max(0, min(idx1.start(), idx2.start()) - num_bases)
                end_idx = min(len(seq), max(idx1.start(), idx2.start()) + num_bases + 1)
                transcript_parts = record.id.split("|")
                base_dict = {
                    "Transcript ID": f"{transcript_parts[0]}|{transcript_parts[1]}",
                    "Gene Name": transcript_parts[5],
                    "Biotype": '|'.join(transcript_parts[7:]),
                    "Strand": record_info["strand"],
                    "Sequence": seq[start_idx:end_idx][::-1] if record_info["strand"] == "-" else str(record.seq),
                    "First seq": idx1.group(),
                    "First seq index": f"{record_info['chrom']}:{record_info['start'] + idx1.start()}",
                    "Second seq": idx2.group(),
                    "Second seq index": f"{record_info['chrom']}:{record_info['start'] + idx2.start()}"
                }
                if not pathogenic_variants:
                    result.append(base_dict)
                else:
                    if record_info["chrom"] in pathogenic_variants:
                        for var_start, var_stop, var_title, _ in pathogenic_variants[record_info["chrom"]]:
                            seq_region_start = record_info['start'] + start_idx
                            seq_region_end = record_info['start'] + end_idx

                            region1_start = record_info['start'] + idx1.start()
                            region1_end = region1_start + len(idx1.group())
                            region2_start = record_info['start'] + idx2.start()
                            region2_end = region2_start + len(idx2.group())

                            if (var_start >= region1_start and var_stop <= region1_end) or (var_start >= region2_start and var_stop <= region2_end):
                                var_distance = 0
                            else:
                                var_distance = min([
                                    abs(region1_start - var_start),
                                    abs(region1_end - var_start),
                                    abs(region2_start - var_start),
                                    abs(region2_end - var_start),
                                    abs(region1_start - var_stop),
                                    abs(region1_end - var_stop),
                                    abs(region2_start - var_stop),
                                    abs(region2_end - var_stop),
                                ]) + 1

                            if (var_start <= seq_region_end) and (var_stop >= seq_region_start):
                                variant_dict = base_dict.copy()
                                variant_dict.update({
                                    "Variant position": f"{record_info['chrom']}:{var_start}-{var_stop}",
                                    "Variant distance": var_distance,
                                    "Variant name": var_title
                                })
                                result.append(variant_dict)
    return result


def process_record_w_transcripts_pc(args):
    """
    Worker function for parallel searches in parse_fasta() for protein-coding transcripts.
    """
    (record, seq1_pattern, seq2_pattern, transcripts,
     pathogenic_variants, window_size, num_bases, 
     gene_list, single_seq) = args

    result = []
    record_transcript_id = record.id.split('|')[0]
    record_info = transcripts.get(record_transcript_id)
    if not record_info:
        return result
    seq = str(record.seq)[::-1] if record_info["strand"] == "-" else str(record.seq)
    seq1_regex = re.compile(seq1_pattern, re.IGNORECASE)

    if gene_list and record_info['gene_name'] not in gene_list:
        return result

    if single_seq:
        for match in seq1_regex.finditer(seq):
            start_idx = max(0, match.start() - num_bases)
            end_idx = min(len(seq), match.start() + len(match.group()) + num_bases)
            transcript_parts = record.id.split("|")
            base_dict = {
                "Transcript ID": f"{transcript_parts[0]}|{transcript_parts[1]}",
                "Gene Name": transcript_parts[5],
                "Regions": '|'.join(transcript_parts[7:]),
                "Strand": record_info["strand"],
                "Sequence": seq[start_idx:end_idx][::-1] if record_info["strand"] == "-" else str(record.seq),
                "Matched seq": match.group(),
                "Matched seq index": f"{record_info['chrom']}:{record_info['start'] + match.start()}"
            }
            if not pathogenic_variants:
                result.append(base_dict)
            else:
                if record_info["chrom"] in pathogenic_variants:
                    for var_start, var_stop, var_title, _ in pathogenic_variants[record_info["chrom"]]:
                        seq_region_start = record_info['start'] + start_idx
                        seq_region_end = record_info['start'] + end_idx

                        region_start = record_info['start'] + match.start()
                        region_end = region_start + len(match.group())

                        if var_start >= region_start and var_stop <= region_end:
                            var_distance = 0
                        else:
                            var_distance = min([
                                abs(region_start - var_start),
                                abs(region_end - var_start),
                                abs(region_start - var_stop),
                                abs(region_end - var_stop),
                            ]) + 1

                        if (var_start <= seq_region_end) and (var_stop >= seq_region_start):
                            variant_dict = base_dict.copy()
                            variant_dict.update({
                                "Variant position": f"{record_info['chrom']}:{var_start}-{var_stop}",
                                "Variant distance": var_distance,
                                "Variant name": var_title
                            })
                            result.append(variant_dict)
    else:
        seq2_regex = re.compile(seq2_pattern, re.IGNORECASE)
        matches_seq1 = list(seq1_regex.finditer(seq))
        matches_seq2 = list(seq2_regex.finditer(seq))
        for idx1, idx2 in itertools.product(matches_seq1, matches_seq2):
            if abs(idx1.start() - idx2.start()) <= window_size:
                start_idx = max(0, min(idx1.start(), idx2.start()) - num_bases)
                end_idx = min(len(seq), max(idx1.start(), idx2.start()) + num_bases + 1)
                transcript_parts = record.id.split("|")
                base_dict = {
                    "Transcript ID": f"{transcript_parts[0]}|{transcript_parts[1]}",
                    "Gene Name": transcript_parts[5],
                    "Regions": '|'.join(transcript_parts[7:]),
                    "Strand": record_info["strand"],
                    "Sequence": seq[start_idx:end_idx][::-1] if record_info["strand"] == "-" else str(record.seq),
                    "First seq": idx1.group(),
                    "First seq index": f"{record_info['chrom']}:{record_info['start'] + idx1.start()}",
                    "Second seq": idx2.group(),
                    "Second seq index": f"{record_info['chrom']}:{record_info['start'] + idx2.start()}"
                }
                if not pathogenic_variants:
                    result.append(base_dict)
                else:
                    if record_info["chrom"] in pathogenic_variants:
                        for var_start, var_stop, var_title, _ in pathogenic_variants[record_info["chrom"]]:
                            seq_region_start = record_info['start'] + start_idx
                            seq_region_end = record_info['start'] + end_idx

                            region1_start = record_info['start'] + idx1.start()
                            region1_end = region1_start + len(idx1.group())
                            region2_start = record_info['start'] + idx2.start()
                            region2_end = region2_start + len(idx2.group())

                            if (var_start >= region1_start and var_stop <= region1_end) or (var_start >= region2_start and var_stop <= region2_end):
                                var_distance = 0
                            else:
                                var_distance = min([
                                    abs(region1_start - var_start),
                                    abs(region1_end - var_start),
                                    abs(region2_start - var_start),
                                    abs(region2_end - var_start),
                                    abs(region1_start - var_stop),
                                    abs(region1_end - var_stop),
                                    abs(region2_start - var_stop),
                                    abs(region2_end - var_stop),
                                ]) + 1

                            if (var_start <= seq_region_end) and (var_stop >= seq_region_start):
                                variant_dict = base_dict.copy()
                                variant_dict.update({
                                    "Variant position": f"{record_info['chrom']}:{var_start}-{var_stop}",
                                    "Variant distance": var_distance,
                                    "Variant name": var_title
                                })
                                result.append(variant_dict)
    return result


def process_record_no_transcripts(args):
    """
    Worker function for parallel searches in parse_fasta() when no transcripts exist.
    """
    (record, seq1_pattern, seq2_pattern, transcripts,
     pathogenic_variants, window_size, num_bases, 
     gene_list, single_seq) = args

    result = []
    seq = str(record.seq)
    seq1_regex = re.compile(seq1_pattern, re.IGNORECASE)

    if single_seq:
        for match in seq1_regex.finditer(seq):
            start_idx = max(0, match.start() - num_bases)
            end_idx = min(len(seq), match.start() + len(match.group()) + num_bases)
            result.append({
                "Transcript ID": record.id,
                "Sequence": seq[start_idx:end_idx],
                "Matched seq": match.group(),
                "Seq index": f"{match.start()+1}"
            })
    else:
        seq2_regex = re.compile(seq2_pattern, re.IGNORECASE)
        matches_seq1 = list(seq1_regex.finditer(seq))
        matches_seq2 = list(seq2_regex.finditer(seq))
        for idx1, idx2 in itertools.product(matches_seq1, matches_seq2):
            if abs(idx1.start() - idx2.start()) <= window_size:
                start_idx = max(0, min(idx1.start(), idx2.start()) - num_bases)
                end_idx = min(len(seq), max(idx1.start(), idx2.start()) + num_bases + 1)
                result.append({
                    "Transcript ID": record.id,
                    "Sequence": seq[start_idx:end_idx],
                    "First seq": idx1.group(),
                    "First seq index": f"{idx1.start()+1}",
                    "Second seq": idx2.group(),
                    "Second seq index": f"{idx2.start()+1}"
                })
    return result


def parse_fasta(fasta_file: str,
                seq1_pattern: str,
                seq2_pattern: Union[str, None],
                pathogenic_variants: Union[dict, None],
                gtf_db: Union[gffutils.FeatureDB, None],
                window_size: int,
                num_bases: int,
                gene_list: Union[set, None],
                is_pc: bool) -> list:
    """
    Parse the (possibly gzipped) FASTA and search for matches.
    A progress bar is displayed to report the status.
    """
    result = []
    if pathogenic_variants == {}:
        return result

    open_func = gzip.open if is_gzipped(fasta_file) else open
    transcripts = extract_transcripts(gtf_db) if gtf_db else None

    with open_func(fasta_file, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))

    single_seq = (seq2_pattern is None or seq2_pattern.strip() == "")

    # Select the appropriate worker function
    if transcripts:
        process_record = process_record_w_transcripts_pc if is_pc else process_record_w_transcripts
    else:
        process_record = process_record_no_transcripts

    for record in tqdm(records, desc="Parsing records"):
        result.extend(process_record((
                    record,
                    seq1_pattern,
                    seq2_pattern,
                    transcripts,
                    pathogenic_variants,
                    window_size,
                    num_bases,
                    gene_list,
                    single_seq
                )))

    return result


def check_sequences(seqs: list) -> list:
    """
    Check up to two sequences for validity.
    """
    chs_allowed = set("ACGTWN*RYSKMBDHV.-")
    seq1 = seqs[0].upper()
    if not set(seq1) <= chs_allowed:
        raise ValueError(f"Sequence {seq1} has unknown base(s).")
    seq2 = seqs[1]
    if seq2 is not None and seq2.strip():
        seq2 = seqs[1].upper()
        if not set(seq2) <= chs_allowed:
            raise ValueError(f"Sequence {seq2} has unknown base(s).")
    else:
        seq2 = None
    return [seq1, seq2]


def parse_hgvs_notation(variation_title: str, start: int, end: int, parser):
    """
    Parse the HGVS notation from the variation title.
    """
    try:
        variant: SequenceVariant = parser(variation_title.split(" ")[0]).hgvs_variant()
        accession = variant.ac
        edit = variant.posedit.edit
        mut_type = variant.posedit.edit.mut_type
        pos_start = start
        pos_end = end
        ref, alt = edit.ref, edit.alt
    except:
        return None

    return {
        'mutation_type': mut_type,
        'pos_start': pos_start,
        'pos_end': pos_end,
        'change': edit,
        'alt': alt
    }


def apply_variant_to_sequence(transcript_seq, transcript_begin, variant_info):
    """
    Apply the variant to the transcript sequence (substitution, deletion, insertion, duplication).
    """
    seq = list(transcript_seq)
    mutation_type = variant_info['mutation_type']
    pos_start = variant_info['pos_start']
    pos_end = variant_info['pos_end']
    alt = variant_info['alt']
    idx = pos_start - transcript_begin
    if idx < 0 or idx >= len(seq):
        return None
    if mutation_type == 'substitution':
        seq[idx] = alt
    elif mutation_type == 'deletion':
        del seq[idx:pos_end]
    elif mutation_type == 'insertion':
        seq.insert(idx, alt)
    elif mutation_type == 'duplication':
        seq.insert(idx + 1, seq[idx])
    elif mutation_type == 'delins':
        del seq[idx:pos_end]
        seq.insert(idx, alt)
    else:
        return None
    return ''.join(seq)


def run_cats(fasta_file,
             seq1,
             seq2,
             output,
             gtf_file=None,
             window_size=5,
             num_bases=25,
             pathogenicity=False,
             snv=False,
             gene_list=None,
             variant_window=None):
    """
    Main analysis function.
    """
    print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tChecking and setting parameters.")

    ext = os.path.splitext(output)[1].lower()
    if ext not in [".csv", ".tsv"]:
        raise ValueError("Output extension not recognised: possible choices are 'csv', 'tsv'.")
    if os.path.isfile(output):
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: WARNING\tOutput file {output} exists. It will be overwritten.")

    seq1, seq2 = check_sequences([seq1, seq2])
    seq1_pattern = create_sequence_pattern(seq1)
    seq2_pattern = None if not seq2 else create_sequence_pattern(seq2)
    if seq2_pattern:
        window_size = max(window_size, len(seq1), len(seq2))

    snv = bool(snv)
    pathogenic = pathogenicity or snv
    variant_window = variant_window or num_bases

    if gene_list:
        if os.path.isfile(gene_list):
            with open(gene_list, 'r') as f:
                gene_list = {line.strip() for line in f if line.strip()}
        else:
            gene_list = {g.strip() for g in gene_list.split(';') if g.strip()}
    else:
        gene_list = None

    gtf_dbfn = None
    if fasta_file == "human":
        fasta_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode.v47.transcripts.fa.gz")
        gtf_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode.v47.annotation.gtf.gz")
        gtf_dbfn = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode_v47.db")
        is_pc = False
    elif fasta_file == "human_pc":
        fasta_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode.v47.pc_transcripts.fa.gz")
        gtf_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode.v47.annotation.gtf.gz")
        gtf_dbfn = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/human/gencode_v47_pc.db")
        is_pc = True
    elif fasta_file == "mouse":
        fasta_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode.vM36.transcripts.fa.gz")
        gtf_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode.vM36.annotation.gtf.gz")
        gtf_dbfn = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode_vM36.db")
        is_pc = False
    elif fasta_file == "mouse_pc":
        fasta_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode.vM36.pc_transcripts.fa.gz")
        gtf_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode.vM36.annotation.gtf.gz")
        gtf_dbfn = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../db/mouse/gencode_vM36_pc.db")
        is_pc = True
    elif not gtf_file and os.path.isfile(fasta_file):
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: WARNING\tWith a custom FASTA and no specified '--gtf', CATS will parse only the sequences; pathogenicity and gene_list will be ignored.")
    elif gtf_file and os.path.isfile(fasta_file):
        gtf_dbfn = os.path.splitext(os.path.abspath(gtf_file))[0] + ".db"
        is_pc = False
    elif not os.path.isfile(fasta_file):
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    elif not os.path.isfile(gtf_file):
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")

    if gtf_file:
        if not os.path.exists(gtf_dbfn):
            print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tCreating GFF database.")
            gffutils.create_db(
                gtf_file,
                dbfn=gtf_dbfn,
                force=True,
                keep_order=True,
                merge_strategy='merge',
                sort_attribute_values=True,
                disable_infer_genes=True,
                disable_infer_transcripts=True
            )
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tConnecting to GFF database.")
        gtf_db = gffutils.FeatureDB(gtf_dbfn)
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tExtracting transcripts.")
        transcripts = extract_transcripts(gtf_db)
    else:
        pathogenic = False
        temp_fasta_file = fasta_file
        gtf_db = None
        gene_list = None
        is_pc = False

    if pathogenic:
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tRetrieving pathogenic variants from ClinVar API.")
        pathogenic_variants = get_pathogenic_variants_from_api(snv=snv, verbose=True, gene_list=gene_list)
    else:
        pathogenic_variants = None

    if pathogenic:
        hgvs_parser = makeGrammar(grammar, globals())
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tCreating temporary FASTA file with variants applied.")
        temp_fasta = tempfile.NamedTemporaryFile(mode='w+', delete=False)
        temp_fasta_file = temp_fasta.name
        open_func = gzip.open if is_gzipped(fasta_file) else open
        with open_func(fasta_file, 'rt') as handle, temp_fasta as temp_fasta_handle:
            records = list(SeqIO.parse(handle, 'fasta'))
            for record in records:
                record_transcript_id = record.id.split('|')[0]
                record_info = transcripts.get(record_transcript_id)
                if record_info is None:
                    continue
                if gene_list is not None and record_info['gene_name'] not in gene_list:
                    continue
                chrom = record_info['chrom']
                if pathogenic_variants is None or chrom not in pathogenic_variants:
                    continue
                for variant in pathogenic_variants[chrom]:
                    variant_start, variant_end, variation_title, variant_type = variant
                    if (variant_start <= record_info['end'] and variant_end >= record_info['start']):
                        variant_info = parse_hgvs_notation(variation_title, variant_start, variant_end, hgvs_parser)
                        if variant_info is None:
                            continue
                        seq = str(record.seq)[::-1] if record_info["strand"] == "-" else str(record.seq)
                        modified_seq = apply_variant_to_sequence(seq, record_info['start'], variant_info)
                        if modified_seq is None:
                            continue
                        modified_seq_oriented = modified_seq[::-1] if record_info["strand"] == "-" else modified_seq
                        new_record = record[:]
                        new_record.seq = Seq(modified_seq_oriented)
                        new_record.id += f'|variant_{variant_type.replace(" ", "-")}_{variation_title.split(" ")[0]}'
                        new_record.description += f" with variant {variation_title}"
                        SeqIO.write(new_record, temp_fasta_handle, 'fasta')
    else:
        if gene_list is not None:
            print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tCreating temporary FASTA file for specified genes.")
            temp_fasta = tempfile.NamedTemporaryFile(mode='w+', delete=False)
            open_func = gzip.open if is_gzipped(fasta_file) else open
            with open_func(fasta_file, 'rt') as handle, temp_fasta as temp_fasta_handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    record_transcript_id = record.id.split('|')[0]
                    record_info = transcripts.get(record_transcript_id)
                    if record_info and record_info['gene_name'] in gene_list:
                        SeqIO.write(record, temp_fasta_handle, 'fasta')
            temp_fasta_file = temp_fasta.name
        else:
            temp_fasta_file = fasta_file

    print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tStarting FASTA parsing process.")
    result_sequences = parse_fasta(
        fasta_file=temp_fasta_file,
        seq1_pattern=seq1_pattern,
        seq2_pattern=seq2_pattern,
        pathogenic_variants=pathogenic_variants,
        gtf_db=gtf_db,
        window_size=window_size,
        num_bases=num_bases,
        gene_list=gene_list,
        is_pc=is_pc
    )

    if not result_sequences:
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: WARNING\tNo matching sequences found.")
    else:
        res = pd.DataFrame(result_sequences).drop_duplicates()
        if pathogenic and "Regions" in res.columns:
            res = res[res.apply(lambda row: row["Variant name"].split(" ")[0] in row["Regions"], axis=1)]
        if pathogenic and "Biotype" in res.columns:
            res = res[res.apply(lambda row: row["Variant name"].split(" ")[0] in row["Biotype"], axis=1)]
        if variant_window != num_bases:
            res = res[res["Variant distance"] <= variant_window]
        os.makedirs(os.path.dirname(output), exist_ok=True)
        with open(output, 'w') as f:
            f.write("# Used settings:\n")
            f.write(f"#\tFASTA:          {fasta_file}\n")
            f.write(f"#\tSeq1:           {seq1}\n")
            if seq2:
                f.write(f"#\tSeq2:           {seq2}\n")
            f.write(f"#\tWindow size:    {window_size}\n")
            f.write(f"#\tNum. bases:     {num_bases}\n")
            f.write(f"#\tPathogenicity:  {pathogenic}\n")
            f.write(f"#\tSNV:            {snv}\n")
            if gene_list:
                f.write(f"#\tGene(s):        {gene_list}\n")
            if pathogenic:
                f.write(f"#\tVariant window: {variant_window}\n")

        if ext == ".csv":
            res.to_csv(output, mode="a", index=False)
        elif ext == ".tsv":
            res.to_csv(output, mode="a", index=False, sep="\t")
        print(f"{strftime('%Y-%m-%d %H:%M:%S', localtime())}: INFO\tResults saved in {output}.")

    if pathogenic or gene_list is not None:
        try:
            os.remove(temp_fasta_file)
        except Exception:
            pass


def main():
    parser = argparse.ArgumentParser(
        description="Parse a FASTA file and find sequences containing one (or two) specified sequences of interest."
    )
    parser.add_argument("--fasta", "-f", dest="fasta_file", required=True,
                        help="Path to the FASTA file or use 'human', 'mouse', 'human_pc' or 'mouse_pc' keyword to access corresponding transcripts.")
    parser.add_argument("--seq1", "-1", dest="seq1", required=True,
                        help="First sequence of interest.")
    parser.add_argument("--seq2", "-2", dest="seq2", required=False,
                        help="(Optional) Second sequence of interest. If omitted, only seq1 will be searched.")
    parser.add_argument("--output", "-o", dest="output", required=True,
                        help="Output file name. Possible extensions: 'csv', 'tsv'")
    parser.add_argument("--window-size", "-w", dest="window_size", default=5, type=int,
                        help="Size of the window around the sequences (for double-sequence mode). Default is 5.")
    parser.add_argument("--num-bases", "-n", dest="num_bases", default=25, type=int,
                        help="Number of preceding and subsequent bases for each output sequence (Default is 25)")
    parser.add_argument("--gtf", "-g", dest="gtf_file", required=False, default=None,
                        help="Path to the GTF file for annotation")
    parser.add_argument("--pathogenicity", "-p", dest="pathogenicity", action='store_true',
                        help="Retrieve only sequences containing potentially pathogenic variants (ClinVar).")
    parser.add_argument("--single-nucleotide-variant", "-snv", dest="snv", action='store_true',
                        help="Retrieve only sequences associated with SNVs from ClinVar. Implies --pathogenicity.")
    parser.add_argument("--gene-list", "-gl", dest="gene_list", required=False,
                        help=("Path to a txt file containing a list of gene names to be analyzed, "
                              "OR a semicolon-separated list of gene names (e.g. 'HBB;HTT')."))
    parser.add_argument("--variant-window", "-vw", dest="variant_window", required=False, type=int,
                        help="Maximum distance between the mutation and the found PAM sequence. Implies --pathogenicity.")
    args = parser.parse_args()

    run_cats(
        fasta_file=args.fasta_file,
        seq1=args.seq1,
        seq2=args.seq2,
        gtf_file=args.gtf_file,
        output=args.output,
        window_size=args.window_size,
        num_bases=args.num_bases,
        pathogenicity=args.pathogenicity,
        snv=args.snv,
        gene_list=args.gene_list,
        variant_window=args.variant_window,
    )

if __name__ == "__main__":
    main()