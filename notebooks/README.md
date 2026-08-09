# 🧬 TP53 Computational Analysis Notebooks

This directory contains the computational notebooks and supporting scripts used in the **Comparative In-Silico Analysis of TP53 Mutation Between Humans and Elephants** project.

The notebooks document the computational workflow from sequence-level comparative analysis through exploratory machine-learning and deep-learning representation.

---

## 🔬 Research Context

TP53 is a conserved tumour-suppressor gene involved in genome integrity, DNA-damage response, cell-cycle regulation, apoptosis, and cancer biology.

This project investigates TP53 sequence conservation and divergence between humans and elephants, with particular interest in the computational characterization of TP53 sequences and mutation-associated regions.

The computational framework combines:

- comparative protein-sequence analysis;
- sequence alignment;
- evolutionary conservation analysis;
- TP53 hotspot mapping;
- exploratory machine learning;
- feature prioritization;
- deep-learning representation;
- reproducible computational workflows.

> **Scientific scope:** The analyses in this repository are computational and exploratory. Sequence conservation or similarity is not, by itself, evidence of equivalent biological function, cancer resistance, or clinical protection.

---

# 📂 Notebook Contents

| File | Purpose |
|---|---|
| `TP53_Comparative_Analysis.ipynb` | Comparative analysis of human and elephant TP53 sequences |
| `EleProtect_Feature_Prioritization.py` | Feature-prioritization workflow associated with the EleProtect analysis |
| `TP53_Deep_Learning_Architecture_and_Predictive_Modelling.ipynb` | Exploratory deep-learning representation and predictive-model architecture |
| `requirements.txt` | Reproducible Python environment specification |

---

# 🧬 1. TP53 Comparative Analysis

### `TP53_Comparative_Analysis.ipynb`

This notebook forms the core comparative sequence-analysis component of the project.

The workflow focuses on:

1. loading TP53 protein sequences;
2. validating sequence inputs;
3. comparing human and elephant TP53 sequences;
4. examining sequence-level conservation;
5. identifying divergent residues;
6. supporting downstream hotspot interpretation;
7. exporting reproducible computational results.

### Input Data

The notebook uses processed sequence data from:

```text
data/processed/TP53_clean.fasta
