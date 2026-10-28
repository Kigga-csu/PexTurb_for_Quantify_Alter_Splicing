# Example dataset

This simulated dataset provides a control group with baseline isoform expression and a case group with altered isoform expression.

| Group label | Role | Samples |
|---|---|---|
| `groups0` | Control | `sample_01`–`sample_05` |
| `groups5` | Case | `sample_26`–`sample_30` |

## Input files

These files correspond to the inputs described in the [main README](../README.md#2-pexturb-input-and-output).

| Input | File | Format |
|---|---|---|
| A: Isoform expression | `data/simulated_tpm.txt` | Tab-separated matrix with transcript IDs in the first column, sample IDs in the header, and TPM values in the remaining columns |
| B: GTF annotation | `data/simulated_annotation.gtf` | Matching gene, transcript, and exon annotations |
| Sample groups for `pexturb_test.py` | `data/anno_group0_vs_group5.txt` | Two tab-separated columns without a header: sample ID and group label |

The sample group file selects the control and case samples from the expression matrix. Its `groups0` rows appear first, so the reported `delta_psi` is mean(control PSI) minus mean(case PSI).

## Run the example

From the repository root:

```bash
bash run_all.sh \
  --tpm example/data/simulated_tpm.txt \
  --gtf example/data/simulated_annotation.gtf \
  --anno example/data/anno_group0_vs_group5.txt \
  --outdir example/output \
  --min-group-cnt 1
```

The workflow generates the exon-pair index, calculates PSI, and performs the control-versus-case differential test using the scripts in `code_pexturb/`.

## Output files

| File in `output/` | Contents |
|---|---|
| `simulated_tpm_index.txt` | Generated exon-pair index |
| `simulated_tpm_psi.txt` | Exon-pair PSI matrix across samples |
| `simulated_tpm_pexturb.txt` | Differential results with feature and gene annotations, `delta_psi`, `pvalue`, and `fdr` |

Column definitions and the current PSI convention are described in [Results interpretation](../README.md#33-results-interpretation). A reference differential result is provided at `expected/group0_vs_group5_pexturb.txt`.

## Auxiliary files

These files are not needed for the command above:

- `data/simulated_index.txt`: precomputed exon-pair index for use with `code_pexturb/calc_psi.py`.
- `data/ground_truth_change_mat.txt`: transcript fold changes used in the simulation.
- `data/simulated_genes.txt`: genes included in the simulation.
- `prepare_subset_gtf.py`: helper used to construct the example GTF subset.
