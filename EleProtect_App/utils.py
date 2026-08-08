# utils.py

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import pandas as pd
import io
import re
from docx import Document

# =========================================================
# CONSTANTS
# =========================================================

HOTSPOTS = [175, 245, 248, 249, 273, 282]

# Realistic mutation frequency weights (literature-based approximation)
TP53_MUTATION_FREQUENCY = {
    175: 0.052,
    245: 0.041,
    248: 0.067,
    249: 0.031,
    273: 0.062,
    282: 0.038
}

TP53_HOTSPOT_INFO = {
    175: "Structural hotspot mutation",
    245: "DNA-binding destabilization",
    248: "DNA contact mutation",
    249: "Structural alteration",
    273: "DNA contact hotspot",
    282: "Structural hotspot"
}

# =========================================================
# SEQUENCE CLEANING
# =========================================================

def clean_sequence(seq: str) -> str:
    s = str(seq).replace("\r", "").replace("\n", "").strip()
    s = re.sub(r"[^A-Za-z]", "", s)
    return s.upper()


def is_nucleotide(seq: str) -> bool:
    if len(seq) == 0:
        return False
    valid = set("ACGTN")
    return sum(1 for c in seq if c in valid) / len(seq) >= 0.75


def translate_if_needed(seq: str) -> str:
    seq = clean_sequence(seq)
    if is_nucleotide(seq):
        try:
            return str(Seq(seq).translate(to_stop=True))
        except Exception:
            return str(Seq(seq).translate())
    return seq


# =========================================================
# FILE PARSING
# =========================================================

def parse_fasta_text(text: str) -> dict:
    sequences = {}
    cur_name = None
    cur_seq = []

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if cur_name:
                sequences[cur_name] = "".join(cur_seq)
            cur_name = line[1:].strip().split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)

    if cur_name:
        sequences[cur_name] = "".join(cur_seq)

    if not sequences:
        sequences["Sequence_1"] = re.sub(r"\s+", "", text)

    return sequences


def parse_docx_bytes(b: bytes) -> dict:
    doc = Document(io.BytesIO(b))
    text = "\n".join(p.text for p in doc.paragraphs)

    if ">" in text:
        return parse_fasta_text(text)

    cleaned = re.sub(r"\s+", "", text)
    if cleaned:
        return {"docx_seq": cleaned}
    return {}


def parse_excel_bytes(b: bytes) -> dict:
    df = pd.read_excel(io.BytesIO(b), sheet_name=0, header=None, engine="openpyxl")
    sequences = {}

    for i, row in df.iterrows():
        values = row.dropna().astype(str).tolist()
        if values:
            sequences[f"row_{i+1}"] = "".join(values)

    return sequences


def parse_upload(uploaded_file):
    filename = uploaded_file.name.lower()
    content = uploaded_file.read()

    if filename.endswith((".fasta", ".fa")) or content.strip().startswith(b">"):
        text = content.decode("utf-8", errors="ignore")
        return parse_fasta_text(text)

    elif filename.endswith(".txt"):
        text = content.decode("utf-8", errors="ignore")
        if ">" in text:
            return parse_fasta_text(text)
        return {"txt_seq": re.sub(r"\s+", "", text)}

    elif filename.endswith(".docx"):
        return parse_docx_bytes(content)

    elif filename.endswith((".xls", ".xlsx")):
        return parse_excel_bytes(content)

    else:
        try:
            text = content.decode("utf-8", errors="ignore")
            if ">" in text:
                return parse_fasta_text(text)
            return {"uploaded_seq": re.sub(r"\s+", "", text)}
        except:
            return {}


# =========================================================
# BLOSUM62 SCORING
# =========================================================

# Load BLOSUM62 matrix once
BLOSUM62 = substitution_matrices.load("BLOSUM62")

def blosum_score(a, b):
    if a == "-" or b == "-":
        return -4
    try:
        return BLOSUM62[a, b]
    except KeyError:
        return -4


# =========================================================
# ALIGNMENT & HOTSPOT MAPPING
# =========================================================

def align_and_map(human_seq: str, query_seq: str, hotspots=None) -> pd.DataFrame:

    if hotspots is None:
        hotspots = HOTSPOTS

    if not human_seq or not query_seq:
        return pd.DataFrame()

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    alignment = next(iter(aligner.align(human_seq, query_seq)))

    aligned_human = alignment[0]
    aligned_query = alignment[1]

    results = []
    human_pos = 0

    for i in range(len(aligned_human)):

        if aligned_human[i] != "-":
            human_pos += 1

        if human_pos in hotspots:

            score = blosum_score(aligned_human[i], aligned_query[i])

            results.append({
                "Codon": human_pos,
                "Human": aligned_human[i],
                "Query": aligned_query[i],
                "Conserved": aligned_human[i] == aligned_query[i],
                "BLOSUM62_Score": score
            })

    return pd.DataFrame(results)


# =========================================================
# ADVANCED FEATURE ENGINEERING
# =========================================================

def build_model_features(df_hotspots: pd.DataFrame) -> dict:

    if df_hotspots.empty:
        return {
            "conservation_score": 0,
            "mutation_frequency": 0,
            "retrogene_variability": 0,
            "sequence_identity": 0,
            "mean_blosum": 0,
            "std_blosum": 0,
            "damaging_fraction": 0,
            "hotspot_count": 0,
            "weighted_mutation_burden": 0
        }

    identity = df_hotspots["Conserved"].mean()
    mean_blosum = df_hotspots["BLOSUM62_Score"].mean()
    std_blosum = df_hotspots["BLOSUM62_Score"].std()

    damaging_fraction = (
        (df_hotspots["BLOSUM62_Score"] < 0).sum()
        / len(df_hotspots)
    )

    freq_vals = df_hotspots["Codon"].map(
        TP53_MUTATION_FREQUENCY
    ).fillna(0)

    mutation_frequency = freq_vals.mean()

    weighted_mutation_burden = (
        freq_vals * abs(df_hotspots["BLOSUM62_Score"])
    ).mean()

    conservation_score = mean_blosum / 11
    retrogene_variability = 1 - conservation_score

    return {
        "conservation_score": float(conservation_score),
        "mutation_frequency": float(mutation_frequency),
        "retrogene_variability": float(retrogene_variability),
        "sequence_identity": float(identity),
        "mean_blosum": float(mean_blosum),
        "std_blosum": float(std_blosum),
        "damaging_fraction": float(damaging_fraction),
        "hotspot_count": len(df_hotspots),
        "weighted_mutation_burden": float(weighted_mutation_burden)
    }


# =========================================================
# HOTSPOT ANNOTATION
# =========================================================

def annotate_hotspots(df_map: pd.DataFrame) -> pd.DataFrame:

    if df_map.empty:
        return df_map

    df_map["Annotation"] = df_map["Codon"].map(TP53_HOTSPOT_INFO)
    return df_map