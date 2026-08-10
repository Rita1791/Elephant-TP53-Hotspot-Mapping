# Reviewer summary

## Question

How do canonical human TP53 cancer-mutation hotspot residues map onto curated
elephant TP53 and TP53-related protein sequences?

## Canonical method

The main script validates protein sequences, identifies human UniProt P04637,
performs explicit global pairwise alignment, maps six human TP53 coordinates,
and reports exact residue identity separately from gaps and unmapped positions.
The human reference is excluded from comparative denominators.

## Reproduce

```bash
python -m pip install -r requirements.txt
python scripts/TP53_Comparative_Analysis.py
pytest -q
```

## Primary evidence

- `data/processed/TP53_clean.fasta`
- `results/tp53_hotspot_mapping.csv`
- `results/tp53_hotspot_identity_summary.csv`
- `results/tp53_comparative_features.csv`
- `results/tp53_summary.json`
- `figures/tp53_hotspot_conservation.png`

## Terminology

“Hotspot conservation” in this repository means exact amino-acid identity at
an alignment-mapped human residue. It is not a formal phylogenetic conservation
score. Exploratory clustering is descriptive and is not a validated biological
classifier.

## Main limitations

- Results depend on sequence selection, annotation quality and alignment
  parameters.
- Sequence identity does not prove shared function or retrogene activity.
- The analysis does not explain the mechanism of elephant cancer resistance or
  independently resolve Peto's paradox.
- No clinical, prognostic or therapeutic claim is supported.

## Research lineage

The associated publication is linked by DOI in `CITATION.cff`. The broader
mammalian TP53 conservation study is a companion project with a separate
dataset and statistical question, not a duplicate result set.

