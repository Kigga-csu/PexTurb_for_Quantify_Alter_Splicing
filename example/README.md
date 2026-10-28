# Example Dataset

This folder contains a runnable simulated example for PexTurb. It is the `seed83` simulation from the study benchmark, chosen because PexTurb performs well at all effect sizes (TPR at FPR = 5%: 0.956 at 1.1×, 0.977 at 1.25×, 0.979 at 1.5×) and the dataset is small enough to ship with the code.

## Directory structure

```
example/
├── README.md                        # this file
├── prepare_subset_gtf.py            # helper that generated simulated_annotation.gtf
├── data/
│   ├── simulated_tpm.txt
│   ├── simulated_annotation.gtf
│   ├── simulated_index.txt
│   ├── anno_group0_vs_group5.txt
│   ├── ground_truth_change_mat.txt
│   └── simulated_genes.txt
└── expected/
    └── group0_vs_group5_pexturb.txt
```

## Input files (`data/`)

| File | Format | Dimensions / size | Description |
|---|---|---|---|
| `simulated_tpm.txt` | TSV, header | 1,041 transcripts × 40 samples (474 KB) | Isoform expression matrix. First column: transcript ID (`ENST...`, no version); remaining 40 columns: sample TPM values produced by Salmon. This is the main PexTurb input. |
| `simulated_annotation.gtf` | GTF | 8,465 lines (3.2 MB) | Subset of Ensembl GRCh38 release 77 containing the 160 simulated genes (gene / transcript / exon records). Used to generate the exon-pair index. |
| `anno_group0_vs_group5.txt` | TSV, no header | 10 rows (179 B) | Sample groups for one comparison: `sample_01`–`sample_05` are group 1 (`groups0`, control), `sample_26`–`sample_30` are group 2 (`groups5`, 1.5× fold change). Each row is `sampleID<TAB>group`. |

### `simulated_tpm.txt`

```text
samples	sample_05	sample_36	sample_09	...
ENST00000336777	691.889740	436.688694	687.817398	...
ENST00000337502	291.386979	284.906979	264.778610	...
```

The first header cell is ignored by PexTurb; column names must match the sample IDs in the annotation file.

### `simulated_annotation.gtf`

Generated from the full `Homo_sapiens.GRCh38.77.sorted.gtf` with:

```bash
python prepare_subset_gtf.py \
    --tpm    data/simulated_tpm.txt \
    --gtf    /path/to/Homo_sapiens.GRCh38.77.sorted.gtf \
    --output data/simulated_annotation.gtf
```

The script finds all genes containing the 1,041 expressed transcripts and keeps every transcript/exon record of those genes, so the resulting index is identical to the one built from the full GTF.

### `anno_group0_vs_group5.txt`

```text
sample_01	groups0
sample_02	groups0
...
sample_26	groups5
...
sample_30	groups5
```

The file has **no header**. The group appearing first (`groups0`) is group 1, so `delta_psi > 0` means a higher retained-pair fraction in `groups0` with the current scripts (equivalently, a lower perturbed-pair fraction).

## Optional / auxiliary files

| File | Description |
|---|---|
| `simulated_index.txt` | Precomputed exon-pair index (5,552 events, 1.1 MB). Byte-identical to the output of `02_gtf_index/gtf_index.py` on the two input files above, so step 2 can be skipped. Format: `chr:exon_up+exon_down`, `gene_id`, `gene_name`, `exon_pari_idx`, `denominator_idx`. |
| `ground_truth_change_mat.txt` | Simulated fold change matrix: 1,038 transcripts × 8 groups. Each row is a transcript; each column is one group. Groups 0–2 are null (all 1×), groups 3–7 are 1.1×, 1.25×, 1.5×, 3× and 5× (160 transcripts changed, one per simulated gene). Not used by PexTurb; provided for benchmarking (e.g. AUROC/TPR against the truth). |
| `simulated_genes.txt` | The 160 gene symbols selected for the simulation, in one line. Provided for provenance. |
| `prepare_subset_gtf.py` | Helper script that produced `simulated_annotation.gtf` from the full GTF. Not part of the PexTurb pipeline. |

## Output file (`expected/`)

| File | Description |
|---|---|
| `group0_vs_group5_pexturb.txt` | Expected differential EPP result for `groups0` vs `groups5`, produced by the original PexTurb test. 5,552 exon pairs + header. Columns: `chr:exon_up+exon_down`, `gene_id`, `gene_name`, `delta_psi`, `pvalue`, `fdr`. |

## How to run

From the repository root:

```bash
bash run_all.sh \
    --tpm    example/data/simulated_tpm.txt \
    --gtf    example/data/simulated_annotation.gtf \
    --anno   example/data/anno_group0_vs_group5.txt \
    --outdir example/output \
    --min-group-cnt 1
```

Expected output (about 5 seconds):

```
example/output/
├── simulated_tpm_index.txt
├── simulated_tpm_psi.txt
└── simulated_tpm_pexturb.txt
```

`simulated_tpm_pexturb.txt` should match `expected/group0_vs_group5_pexturb.txt` row by row: identical `delta_psi` and a maximum p-value difference below 2.1 × 10⁻¹⁴.

## Provenance

- Simulation dataset: `/data5/lab/wangshixian/IsoMantis_simulation/wsx_test/seed83`
- Annotation: Ensembl GRCh38 release 77 (`Homo_sapiens.GRCh38.77.sorted.gtf`)
- Quantification: Salmon 1.4.0
- The original FASTQ files are about 3.2 GB per seed and are not included.
