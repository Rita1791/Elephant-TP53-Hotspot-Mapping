# 📓 Computational Notebooks

This directory contains the interactive computational notebooks
associated with the Elephant TP53 comparative research project.

The notebooks provide an inspectable record of selected
computational analyses and support exploration of the sequence,
comparative, and EleProtect-related components of the project.

---

## 🧬 Notebook Overview

| Notebook | Purpose |
|---|---|
| `TP53_AI_Pipeline.ipynb` | Computational TP53 analysis workflow |
| `EleProtect_ML_Analysis.ipynb` | Machine-learning analysis associated with the EleProtect component |

---

# 🔬 `TP53_AI_Pipeline.ipynb`

This notebook contains the computational workflow associated
with the TP53 analysis pipeline.

The notebook should be interpreted together with:

- `data/raw/`
- `data/processed/`
- `data/Database/`
- `scripts/`
- `Results/`
- `figures/`

The exact analytical procedures, parameters, and outputs are
defined by the executable notebook cells.

---

# 🐘 `EleProtect_ML_Analysis.ipynb`

This notebook contains computational analysis associated with
the EleProtect research/application component.

Its outputs should be interpreted in the context of the
underlying TP53 research dataset and the EleProtect prototype.

---

# 🔄 Relationship to the Repository

The computational workflow is organised conceptually as:

```text
Public Biological Resources
          │
          ▼
       data/raw/
          │
          ▼
     Data Processing
          │
          ▼
    data/processed/
          │
          ▼
    Structured Dataset
     data/Database/
          │
          ▼
       Notebooks
          │
          ▼
        Results
          │
          ▼
       Figures
