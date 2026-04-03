# Elephant TP53 Hotspot Mapping

## Overview
This repository presents a reproducible computational framework for analyzing evolutionary conservation of TP53 mutation hotspots across human and elephant genomes.

## Scientific Context
TP53 is a critical tumor suppressor gene frequently mutated in human cancers. Despite increased body mass and lifespan, elephants exhibit reduced cancer incidence (Peto’s paradox). This study investigates whether canonical TP53 mutation hotspots are evolutionarily conserved and how retrogene copies contribute to sequence variability.

## Research Objectives
- Map canonical human TP53 mutation hotspots (R175, G245, R248, R249, R273, R282)
- Evaluate conservation across elephant canonical sequences
- Assess variability across TP53 retrogenes
- Develop a composite scoring framework for hotspot prioritization

## Methodological Framework
### Data Acquisition
- UniProt (P04637 – Human TP53)
- NCBI protein database (elephant TP53 sequences)

### Computational Pipeline
1. Sequence preprocessing and quality control
2. Multiple sequence alignment (MAFFT)
3. Residue-level hotspot mapping
4. Conservation analysis
5. Composite scoring model:
   - Human Mutation Frequency (HMF)
   - Cross-Species Conservation (CSC)
   - Retrogene Variability Index (RVI)
6. Phylogenetic analysis (MEGA)

## Key Findings
- Canonical TP53 hotspot residues are strongly conserved across elephant species
- Retrogene copies exhibit position-specific variability
- High-priority residues identified: R175, G245, R248, R273
- Evidence supports evolutionary constraint within the DNA-binding domain

## Repository Structure
