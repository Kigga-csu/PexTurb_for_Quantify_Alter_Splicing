# PexTurb

## Description

PexTurb quantifies **exon-pair perturbation (EPP)** from transcript annotations and isoform expression levels. For each gene, it enumerates exon pairs across annotated isoforms in a GTF file. For a target exon pair, isoforms in which at least one target exon is elongated, shortened, or absent are defined as **inclusion isoforms**, whereas isoforms preserving both target exons are defined as **exclusion isoforms**.

EPP is quantified using percent spliced in (PSI):

$$
\mathrm{PSI} = \frac{Q_{\mathrm{inclusion}}}{Q_{\mathrm{inclusion}} + Q_{\mathrm{exclusion}}},
$$

where $Q_{\mathrm{inclusion}}$ and $Q_{\mathrm{exclusion}}$ are the summed expression levels of inclusion and exclusion isoforms, respectively. PSI therefore represents the fraction of exon-pair-associated expression contributed by isoforms with at least one perturbed target exon.

Applying this calculation across exon pairs and samples produces an **exon-pair-by-sample PSI matrix**. This matrix provides a quantitative representation of alternative splicing for differential EPP analysis and downstream applications such as EPP-based splicing quantitative trait locus (epp-sQTL) mapping.

![Overview of the PexTurb method](figures/figure1.png)

**Figure 1. Overview of PexTurb.** **a**, Motivation for exon-pair-based splicing analysis. **b**, Definition of inclusion and exclusion isoforms and calculation of EPP PSI. **c**, Workflow from transcript annotation and isoform expression to a PSI matrix and downstream analyses. [View the original figure as a PDF](figures/figure1.pdf).

## 1. Installation

Create and activate a conda environment with the required Python packages:

```bash
conda create -n pexturb --override-channels -c conda-forge \
  python=3.10 numpy=1.23.3 scipy=1.9.1 statsmodels=0.13.2 -y
conda activate pexturb
```

If starting from FASTQ files, also install Salmon:

```bash
conda install --override-channels -c conda-forge -c bioconda \
  --strict-channel-priority salmon=1.4.0 -y
```

## 2. PexTurb input and output

Conceptually, PexTurb uses an isoform expression matrix **A** and an exon-pair index **B** to generate a PSI matrix. Sample group labels are additionally required for differential analysis:

```text
B = generate_index(GTF, A)
X = PexTurb(A, B, sample_groups)
```

This notation summarizes the workflow; the executable interface consists of the command-line scripts shown in Section 3. The index is generated from the GTF and uses the transcript IDs and row order in **A** to locate expression values.

### A: Isoform expression

**A** is a tab-separated isoform expression matrix with a header. The first column contains transcript IDs; each subsequent column contains expression values for one sample. Transcript IDs and sample IDs must be unique, and expression values must be finite and nonnegative. TPM or FPKM values can be used; the supplied example contains Salmon-derived TPM values.

```text
transcript_id	sample_01	sample_02	sample_03
ENST00000336777	691.889740	436.688694	687.817398
ENST00000337502	291.386979	284.906979	264.778610
```

The first header cell is a label; the other header cells are sample IDs. Version suffixes in expression-matrix transcript IDs, such as `.3` in `ENST00000336777.3`, are removed when matching transcripts to the GTF.

### B: Exon-pair index

**B** links each exon-pair feature to the relevant isoforms in **A**. Generate it using the transcriptome annotation corresponding to the expression matrix:

```bash
python 02_gtf_index/gtf_index.py \
  --tpm expression.tsv \
  --gtf annotation.gtf \
  --index index.tsv
```

The GTF should contain `gene`, `transcript`, and `exon` records. The example uses Ensembl GRCh38 release 77. The current parser supports chromosomes `1`–`22`, `X`, and `Y`, with or without the `chr` prefix.

The index is a tab-separated file with five columns:

| Column | Contents |
|---|---|
| `chr:exon_up+exon_down` | Chromosome and coordinates of the exon pair; equivalent pairs can be grouped with semicolons |
| `gene_id` | Gene identifier |
| `gene_name` | Gene name, or gene ID when a name is unavailable |
| `exon_pari_idx` | Comma-separated, zero-based expression-row indices for isoforms retaining both target exons; the spelling matches the file header |
| `denominator_idx` | Comma-separated, zero-based expression-row indices for isoforms sharing at least one of the four target exon-boundary coordinates |

An index can be reused when transcript identities, row order, and annotation are unchanged. Regenerate it after filtering or reordering transcript rows. A precomputed example index is provided at `example/data/simulated_index.txt`.

### Sample group annotation

Differential analysis requires a tab-separated file with two columns and **no header**. Each row gives a sample ID and its group. Supply exactly two groups; sample IDs must match columns in **A**.

```text
sample_01	control
sample_02	control
sample_26	treated
sample_27	treated
```

The group appearing first is group 1. The reported effect size is `delta_psi = mean(PSI in group 1) − mean(PSI in group 2)`. The expression matrix may include additional samples; the differential test uses the samples listed for the two groups.

### C: Output — X = PexTurb(A, B, sample_groups)

**X** is a differential result table containing exon-pair identifiers, gene annotations, `delta_psi`, `pvalue`, and `fdr`. The workflow also saves the exon-pair-by-sample PSI matrix and the index used to calculate it.

| Output | Contents |
|---|---|
| Exon-pair index | Features and their expression-row indices |
| PSI matrix | Exon-pair features as rows and samples as columns, preceded by three feature annotation columns |
| Differential result **X** | Effect size and statistical significance for each tested feature |

Output filenames and column definitions are explained in Section 3.3.

## 3. Example usage

### 3.1 Example input data

The repository includes a **simulated RNA-seq dataset** (`seed83`) containing 1,041 expressed transcripts from 160 simulated genes and 40 samples (eight groups with five replicates each). The example compares five `groups0` samples with five `groups5` samples, corresponding to the 1.5-fold-change simulation setting.

The example begins with the Salmon-derived TPM matrix. The original simulated FASTQ files are not included.

| File | Description |
|---|---|
| `example/data/simulated_tpm.txt` | Isoform expression matrix: 1,041 transcripts × 40 samples |
| `example/data/simulated_annotation.gtf` | GTF subset covering the 160 simulated genes |
| `example/data/anno_group0_vs_group5.txt` | Two-group sample annotation: five control and five treated samples |
| `example/data/simulated_index.txt` | Optional precomputed index containing 5,552 exon-pair features |
| `example/data/ground_truth_change_mat.txt` | Simulated transcript fold changes for benchmarking |
| `example/data/simulated_genes.txt` | Gene symbols selected for the simulation |
| `example/expected/group0_vs_group5_pexturb.txt` | Reference differential result for the example comparison |

See [example/README.md](example/README.md) for further details about the example data.

### 3.2 Running PexTurb

From the repository root, users can run the following command to analyze the example isoform expression matrix. It generates the exon-pair index, calculates the PSI matrix, and performs differential analysis:

```bash
bash run_all.sh \
  --tpm example/data/simulated_tpm.txt \
  --gtf example/data/simulated_annotation.gtf \
  --anno example/data/anno_group0_vs_group5.txt \
  --outdir example/output \
  --min-group-cnt 1
```

The three output files are written to `example/output/`. The `--min-group-cnt` argument specifies the minimum number of non-missing PSI values required **in each group** for a feature to be tested; the default is `1`.

<details>
<summary>Run the individual steps or use the precomputed index</summary>

Create the output directory and generate the index:

```bash
mkdir -p example/output

python 02_gtf_index/gtf_index.py \
  --tpm example/data/simulated_tpm.txt \
  --gtf example/data/simulated_annotation.gtf \
  --index example/output/simulated_tpm_index.txt
```

Calculate the PSI matrix:

```bash
python 03_psi/calc_psi.py \
  --tpm example/data/simulated_tpm.txt \
  --index example/output/simulated_tpm_index.txt \
  --psi example/output/simulated_tpm_psi.txt
```

To use the precomputed index, skip index generation and replace the `--index` path above with `example/data/simulated_index.txt`.

Perform the differential test:

```bash
PYTHONPATH=02_gtf_index python 04_pexturb_test/pexturb_test.py \
  --psi example/output/simulated_tpm_psi.txt \
  --anno example/data/anno_group0_vs_group5.txt \
  --output example/output/simulated_tpm_pexturb.txt \
  --min_group_cnt 1
```

The standalone test uses `--min_group_cnt` with underscores; `run_all.sh` exposes this setting as `--min-group-cnt`.

</details>

<details>
<summary>Prepare an isoform expression matrix from your own FASTQ files</summary>

Build a Salmon index using a GENCODE transcript FASTA corresponding to your GTF:

```bash
bash 01_salmon/build_salmon_index.sh \
  gencode.transcripts.fa.gz \
  salmon_index
```

For a decoy-aware index, provide a FASTA containing both transcript and decoy sequences and pass the matching decoy-name list as the third argument:

```bash
bash 01_salmon/build_salmon_index.sh \
  gentrome.fa.gz \
  salmon_index \
  decoys.txt
```

Prepare `samples.tsv` with no header. Use three tab-separated columns for paired-end samples or two for single-end samples:

```text
sample_01	/path/sample_01_R1.fastq.gz	/path/sample_01_R2.fastq.gz
sample_02	/path/sample_02.fastq.gz
```

Quantify the samples and combine their `quant.sf` files:

```bash
bash 01_salmon/run_salmon.sh \
  --samples samples.tsv \
  --index salmon_index \
  --outdir salmon_output \
  --threads 8

python 01_salmon/build_tpm_matrix.py \
  --quant-root salmon_output \
  --output expression.txt
```

By default, a transcript is removed if its TPM is below 1 in at least 50% of the samples. Add `--keep-all` to the matrix-building command to retain all transcripts shared by the samples.

Then analyze your expression matrix with a matching GTF and two-group sample annotation:

```bash
bash run_all.sh \
  --tpm expression.txt \
  --gtf annotation.gtf \
  --anno sample_groups.tsv \
  --outdir results
```

</details>

### 3.3 Results interpretation

The example generates:

```text
example/output/
├── simulated_tpm_index.txt
├── simulated_tpm_psi.txt
└── simulated_tpm_pexturb.txt
```

#### Output file explanation

| File | Explanation |
|---|---|
| `simulated_tpm_index.txt` | The exon-pair index **B**, generated from the example expression matrix and GTF |
| `simulated_tpm_psi.txt` | PSI matrix with feature identifiers, `gene_id`, and `gene_name` in the first three columns, followed by the 40 sample columns; undefined PSI values are written as `NA` |
| `simulated_tpm_pexturb.txt` | Differential result **X** for `groups0` versus `groups5`, with one row per tested exon-pair feature |

**X contains the following columns:**

| Column | Meaning |
|---|---|
| `chr:exon_up+exon_down` | Exon-pair coordinates; semicolons separate pairs represented by the same feature |
| `gene_id` | Gene identifier from the GTF |
| `gene_name` | Gene name from the GTF, or gene ID if no name is available |
| `delta_psi` | Mean PSI in group 1 minus mean PSI in group 2 |
| `pvalue` | Two-sided Welch t-test p-value after clipping PSI to [0.001, 0.999] and applying the logit transformation |
| `fdr` | Benjamini–Hochberg adjusted p-value across tested features |

The supplied reference result contains **5,552 tested features**, of which **378 have FDR < 0.05**. These counts describe the stored example result; statistical significance alone does not measure accuracy against the simulation ground truth. Features with too few non-missing values in either group are excluded from testing.

**Current output convention.** The manuscript defines EPP PSI as the perturbed-isoform fraction, as shown in Description and Figure 1. The current `calc_psi.py` instead writes the retained-pair fraction, which is its complement: `EPP PSI = 1 − output PSI`. Consequently, a positive `delta_psi` in the current output indicates greater pair retention in group 1; the corresponding EPP change has the opposite sign. This distinction applies to interpreting effect direction.

## Repository Contents

```text
PexTurb/
├── README.md
├── run_all.sh                         # index → PSI → differential test
├── 01_salmon/
│   ├── build_salmon_index.sh           # build the Salmon index
│   ├── run_salmon.sh                   # quantify paired- or single-end samples
│   └── build_tpm_matrix.py             # merge quant.sf files into a TPM matrix
├── 02_gtf_index/
│   └── gtf_index.py                    # generate exon-pair indices
├── 03_psi/
│   └── calc_psi.py                     # calculate the PSI matrix
├── 04_pexturb_test/
│   └── pexturb_test.py                 # differential analysis and FDR correction
├── example/
│   ├── README.md
│   ├── prepare_subset_gtf.py           # helper for constructing the example GTF
│   ├── data/                          # example inputs and simulation truth
│   └── expected/                      # reference differential result
└── figures/
    ├── figure1.png                     # method overview displayed in this README
    └── figure1.pdf                     # original PDF extracted from the manuscript
```

## 4. Contact

For questions or bug reports, please contact [234701007@csu.edu.cn](mailto:234701007@csu.edu.cn).
