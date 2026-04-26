# 🧬 Structure-Aware Screening

## 🚨 Problem

Modern DNA screening systems (e.g. those used by synthesis providers or benchtop devices) rely heavily on **sequence similarity** methods such as BLAST, DIAMOND, and HMM-based searches.

These approaches work well for identifying known sequences and close homologs.

However, they assume:

> **functionally similar proteins must also be sequence-similar**

This assumption is increasingly fragile.

---

## ⚠️ Emerging risk: AI-generated biological sequences

Recent advances in generative biology — including protein language models, diffusion models, and structure-based design — enable the creation of **novel protein sequences** that:

- have **low or no sequence similarity** to known proteins  
- but still **fold into known or functional structures**

This creates a potential gap in current screening systems:

```text
sequence dissimilar → passes screening
structure similar   → functionally relevant
```

In other words:

> **A sequence can evade sequence-based screening while retaining meaningful structural similarity to known proteins.**

---

## 🧠 Key insight

Protein function is more closely related to **structure** than sequence.

Two proteins can have:
- <25% sequence identity  
- but highly similar 3D folds  

This is a well-known phenomenon in structural biology, but it is **not explicitly handled** in most screening pipelines.

---

## 💡 Our approach

We propose a **structure-aware extension** to existing DNA/protein screening workflows.

Instead of relying solely on sequence similarity, we add a second layer:

```text
sequence → (optional structure prediction) → structural similarity search
```

Specifically:

```text
Input sequence
   ↓
(optional) structure prediction (ColabFold)
   ↓
Query structure
   ↓
Foldseek structural search against reference database
   ↓
Similarity scoring (TM-score, coverage, sequence identity)
   ↓
Structural alignment (US-align)
   ↓
Interactive review (dashboard)
```

---

## 🔌 Intended use case

This module is designed as a **second-pass escalation step** in screening pipelines such as:

- DNA synthesis provider screening workflows  
- benchtop DNA synthesiser safety software  
- open-source screening frameworks (e.g. IBBIS / Common Mechanism)

It should be triggered when:

```text
sequence-based screening is ambiguous or inconclusive
```

rather than run on every sequence (to manage computational cost).

## ⚡ Quick result (TL;DR)

Our structure-aware module successfully identifies proteins that may likely be missed by sequence-based screening (would need further investigation to confirm).

Example (Case Study 2):

```text
Query: 1QO2 (TIM barrel enzyme)

Hit:   1J5T
TM-score:        0.526
Query coverage:  0.980
Sequence identity: 0.109
```

Interpretation:

```text
Low sequence similarity (11%)
High structural similarity (TM-score > 0.5)
→ correctly flagged as structure_similar_sequence_diverse
```

This demonstrates:

> A sequence can pass on sequence similarity, but still be structurally similar to known proteins.
---

## 🧩 Key takeaway

> **Sequence similarity is not sufficient for screening in the age of AI-generated biology.  
Structure provides an additional, orthogonal signal that can close this gap.**

---

## Pipeline Details

Given one or more protein sequences in a FASTA file, the pipeline can:

1. Get a predicted structure with ColabFold.
2. Search each query structure against a reference structure database (SoC) using Foldseek.
3. Score hits using structural and sequence similarity metrics.
4. Align query and reference structures using US-align.
5. Launch a Streamlit dashboard for interactive review.

---

## Key metrics

| Metric | Meaning |
|---|---|
| `qtmscore` | TM-score normalised by the query length. Main structural similarity metric. |
| `qcov` | Fraction of the query structure covered by the alignment. |
| `tcov` | Fraction of the target/reference structure covered by the alignment. |
| `fident` | Fractional sequence identity over the aligned region. |
| `rmsd` | Root mean square deviation after structural alignment. Lower is better. |
| `evalue` | Statistical significance of the Foldseek hit. Lower is better. |

The most important case for this project is:

```text
high structural similarity + low sequence identity
```

In the current scoring script, this is labelled:

```text
structure_similar_sequence_diverse
```

---

## Scoring logic

The current simple scoring rule is:

```python
if qtmscore >= 0.5 and qcov >= 0.5 and fident <= 0.30 and evalue <= 1e-3:
    similarity_class = "structure_similar_sequence_diverse"
elif qtmscore >= 0.7 and qcov >= 0.7 and evalue <= 1e-5:
    similarity_class = "high_similarity"
elif qtmscore >= 0.5 and qcov >= 0.5 and evalue <= 1e-3:
    similarity_class = "moderate_similarity"
elif qtmscore >= 0.3:
    similarity_class = "weak_similarity"
else:
    similarity_class = "low_similarity"
```

Priority order in the dashboard:

```text
structure_similar_sequence_diverse
high_similarity
moderate_similarity
weak_similarity
low_similarity
```

---

## Installation

Create and activate a conda environment:

```bash
conda create -n foldseek_env python=3.10 -y
conda activate foldseek_env
```

Install Foldseek and US-align:

```bash
conda install -c conda-forge foldseek -y
conda install -c bioconda usalign -y
```

Install Python dependencies:

```bash
pip install pandas streamlit plotly py3Dmol requests
```

Optional, if using ColabFold:

```bash
# Install LocalColabFold separately following its installation instructions.
# The pipeline expects the command `colabfold_batch` to be available.
```

---

## Project structure

Expected layout:

```text
structure_screen/
  run_pipeline.py
  example.fasta

  scripts/
    download_clean_pdbs.py
    score_foldseek_hits.py
    align_hit_structures.py
    app.py

  inputs/
    query_pdbs/
      protein1.pdb
      protein2.pdb

    reference_structures/
      reference1.pdb
      reference2.pdb
      reference3.pdb

  outputs/
    protein1/
      query.pdb
      foldseek_results.tsv
      scored_similarity_hits.csv
      similarity_report.txt
      aligned_structures/

    protein2/
      query.pdb
      foldseek_results.tsv
      scored_similarity_hits.csv
      similarity_report.txt
      aligned_structures/
```

---

## Running the demo

### 1. Download example PDB structures

```bash
python scripts/download_clean_pdbs.py
```

This creates example query and reference structures under:

```text
inputs/query_pdbs/
inputs/reference_structures/
```

---

### 2. Create a FASTA file

Example:

```fasta
>protein1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAN
```

In hackathon mode, the sequence itself is only used to identify the query name. The structure is loaded from:

```text
inputs/query_pdbs/protein1.pdb
```

So the FASTA ID must match the PDB filename.

For example:

```text
>protein1
```

matches:

```text
inputs/query_pdbs/protein1.pdb
```

---

### 3. Run without ColabFold

This is the recommended mode for fast testing:

```bash
python run_pipeline.py --fasta example.fasta
```

This will:

```text
read FASTA
copy existing query PDB
build Foldseek database
run Foldseek search
score hits
align hits with US-align
launch Streamlit dashboard
```

---

### 4. Run with ColabFold

If ColabFold is installed and available:

```bash
python run_pipeline.py --fasta example.fasta --run-colabfold
```

This will predict a structure for each FASTA sequence before running the structural search.

---

### 5. Launch dashboard manually

If you ran the pipeline with `--skip-dashboard`, launch the dashboard manually:

```bash
streamlit run scripts/app.py
```

Then open:

```text
http://localhost:8501
```

---

## Dashboard

The Streamlit dashboard provides:

- query selector for multi-query runs
- summary cards for hit classes
- scatter plot of sequence identity vs structural similarity
- ranked hit table
- metric explanations
- individual hit inspection
- query structure viewer
- reference structure viewer
- aligned query/reference overlay

The most useful plot is:

```text
x-axis: sequence identity
y-axis: query-normalised TM-score
```

The key region is:

```text
low sequence identity
high TM-score
```

These are the cases where structure search adds value beyond sequence search.

---

## Case Study 1: Serine protease sanity check

### Rationale

The first case study is a sanity check using well-known structurally related enzymes.

We used:

| Role | Example |
|---|---|
| Query | Trypsin |
| Related reference | Elastase |
| Related reference | Chymotrypsin |
| Negative control | Lysozyme |

This case study asks whether the pipeline can detect obvious structural homologs and avoid reporting an unrelated protein.

### Expected result

```text
trypsin vs elastase       strong structural similarity
trypsin vs chymotrypsin   moderate/strong structural similarity
trypsin vs lysozyme       weak or no structural similarity
```

### Observed result

The pipeline detected the known serine protease similarities.

Example results:

| Target | TM-score | Query coverage | Sequence identity | Interpretation |
|---|---:|---:|---:|---|
| Elastase | ~0.94 | ~1.00 | ~0.38 | Strong structural match |
| Chymotrypsin | ~0.54 | ~0.56 | ~0.42 | Moderate structural match |
| Lysozyme | Not reported | — | — | Negative control not falsely flagged |

### Interpretation

This confirms that the basic Foldseek search and scoring pipeline works:

```text
known structural homologs are detected
unrelated negative control is not prioritised
```

---
## Case Study 2: Sequence-diverse structural similarity in TIM/barrel-like enzymes

### Rationale

The second case study tests the main purpose of the pipeline: detecting proteins that are structurally similar despite low sequence similarity.

This is important because sequence-based methods can miss relationships between proteins that have diverged substantially at the amino-acid level but retain a similar 3D fold.

We used a TIM/barrel-like enzyme query and searched it against a reference set containing other barrel-like enzymes and unrelated controls.

### Query and reference selection

The query structure used in this case study is:

```text
PDB ID: 1QO2
Protein: HisA (ProFAR isomerase)
Organism: Thermotoga maritima
Fold: TIM barrel (β/α barrel)
```

The reference set includes:

| PDB ID | Description | Role |
|---|---|---|
| 1THF | HisF (imidazole glycerol phosphate synthase cyclase subunit) | Intended sequence-diverse structural control |
| 1J5T | Indole-3-glycerol phosphate synthase | TIM barrel enzyme (alternative structural match) |
| 1PII | PRA isomerase / β/α barrel domain protein | Barrel-like fold |
| 2CDS | Lysozyme | Unrelated negative control |

### Query sequence

For completeness, the amino acid sequence of the query protein (1QO2) is:

```text
[sequence omitted here for brevity — derived from PDB structure 1QO2]
```

In hackathon mode, the pipeline uses the precomputed structure:

```text
inputs/query_pdbs/protein1.pdb
```

and the FASTA entry:

```fasta
>protein1
...
```

is used only to map the query ID to the corresponding structure.

### Expected result

The desired behaviour was:

```text
TM-score >= 0.5
query coverage >= 0.5
sequence identity <= 0.30
```

Hits meeting these criteria are labelled:

```text
structure_similar_sequence_diverse
```

### Observed results

The pipeline successfully identified two structure-similar but sequence-diverse hits:

| Target | TM-score | Query coverage | Target coverage | Sequence identity | E-value | Classification |
|---|---:|---:|---:|---:|---:|---|
| 1J5T indole-3-glycerol phosphate synthase TIM barrel | 0.526 | 0.980 | 0.880 | 0.109 | 8.39e-06 | structure_similar_sequence_diverse |
| 1PII PRA isomerase beta/alpha barrel | 0.507 | 0.956 | 0.434 | 0.092 | 6.43e-04 | structure_similar_sequence_diverse |

### Interpretation

These results meet the core requirement of the case study.

Both hits have:

```text
moderate structural similarity
high query coverage
very low sequence identity
significant E-values
```

This is exactly the regime where structure-aware search adds value beyond sequence-only screening.

The strongest example is `1J5T`, which aligns across almost the full query structure:

```text
query coverage = 0.980
sequence identity = 0.109
TM-score = 0.526
```

This indicates that the query and reference share a similar overall fold despite only ~11% sequence identity over the aligned region.

### Honest caveat

The originally intended HisA/HisF control structure, `1THF`, was not reported in this output under the current Foldseek settings. However, the pipeline still demonstrates the intended behaviour because other barrel-like structures were recovered with:

```text
TM-score > 0.5
sequence identity < 0.30
```

Therefore, the case study is successful as a demonstration of sequence-diverse structural similarity, although the exact reference hit differed from the initial expectation.

### Conclusion

This case study demonstrates the main value of the pipeline:

```text
sequence similarity is low
but structural similarity is still detectable
```

The results support the use of structure-aware screening as a complementary layer to sequence-based methods.

---

## Limitations

Current limitations:

- ColabFold is optional because structure prediction is computationally expensive.
- Foldseek reports only hits passing its internal thresholds, so unrelated controls may not appear in the output.
- Multi-chain PDB files can produce multiple hits for the same structure.
- The current scoring thresholds are simple heuristics.
- Structural similarity does not necessarily imply identical function.
- Domain-level matches may require additional logic.

---

## Future work

Useful next steps:

- Add chain deduplication.
- Add domain-level structural search.
- Cache Foldseek databases.
- Add support for larger reference databases.
- Add confidence scores from predicted structures, such as pLDDT.
- Integrate as an optional second-pass module in a broader screening pipeline.
- Improve scoring thresholds using a benchmark set.

---

## Summary

This project demonstrates a structure-aware protein screening workflow.

The main insight is:

```text
sequence similarity is not always enough
```

By combining Foldseek, US-align, and an interactive dashboard, the pipeline can highlight proteins that are:

```text
structurally similar but sequence-diverse
```

This provides an additional signal beyond standard sequence-based screening.