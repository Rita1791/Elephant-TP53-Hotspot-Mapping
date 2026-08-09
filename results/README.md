# Results

This directory contains the computational outputs generated during the comparative analysis of **TP53 between humans and elephants**, including sequence alignment, hotspot mapping, phylogenetic analysis, feature extraction, and exploratory machine-learning analyses.

The files are organized according to the major analytical stages of the project. Together, these outputs provide a traceable link between the input datasets, computational workflows, visualizations, and the biological observations reported in the associated research.

---

## 📁 Results Structure

```text
results/
│
├── MEGA/
│   ├── TP53_MSA.mas
│   ├── TP53_MSA.meg
│   └── TP53_tree.nwk
│
├── ML/
│   ├── TP53_hotspot_analysis.csv
│   ├── basic_features.csv
│   ├── identity_barplot.png
│   ├── top10_tp53_like.csv
│   ├── tp53_features_with_similarity.csv
│   └── tp53_features_with_similarity_clustered.csv
│
├── MSA/
│   ├── TP53_MSA.fasta
│   └── TP53_MSA_logo.png
│
├── MSA_1and2/
│   ├── aligned.fasta
│   ├── alignment.png
│   ├── 2_alignment.png
│   ├── alignment.jvp
│   └── clustalo.aln
│
├── MSA_3/
│   ├── TP53_3seq_alignment.fasta
│   ├── TP53_3seq_alignment.aln
│   └── TP53_3seq_MSA.png
│
└── phylogeny/
    ├── TP53_tree.nwk
    ├── TP53_tree.png
    ├── TP53_tree.emf
    ├── Guide_Tree.tree
    └── Phylogenetic_Tree.phylotree
