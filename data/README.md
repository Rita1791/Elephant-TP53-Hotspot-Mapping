# 🧬 Research Data

This directory contains the biological sequence resources,
processed comparative datasets, structured TP53 research
data, and provenance information supporting the project.

The data architecture separates:

> **Raw biological resources → Processed datasets → Structured research database**

---

# 📂 Data Architecture

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
├── processed/
│   ├── TP53_all_sequences.fasta
│   ├── TP53_clean.fasta
│   ├── human_elephant_tp53_pair.fasta
│   └── human_elephant_tp53_retrogene_comparison.fasta
│
├── Database/
│   ├── tp53_elephant_database.csv
│   ├── tp53_elephant_database.json
│   └── README.md
│
└── README.md
