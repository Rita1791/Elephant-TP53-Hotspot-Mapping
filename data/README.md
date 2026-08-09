# 🧬 Research Data

This directory contains the biological sequence resources,
reference metadata, and processed comparative datasets used
in the TP53 comparative genomics analysis.

The data architecture follows a simple principle:

> **Raw biological resources are separated from derived
> computational datasets.**

---

## 📂 Data Architecture

```text
data/
│
├── raw/
│   │
│   ├── Blast/
│   ├── ncbi_dataset/
│   ├── african_tp53_hits.fasta
│   ├── asian_tp53_hits.fasta
│   ├── elephant_assembly_data_report.jsonl
│   ├── elephant_tp53_RTG9.fasta
│   ├── elephant_tp53_RTG9_translate.fasta
│   ├── elephant_tp53_canonical.fasta
│   ├── human_tp53.fasta
│   └── provenance.txt
│
└── processed/
    ├── human_elephant_tp53_pair.fasta
    ├── human_elephant_tp53_pair_2seq.fasta
    └── human_elephant_tp53_retrogene_comparison.fasta
