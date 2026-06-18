# Step-by-Step Usage

The easiest way to run M-SCAN is with a config file:

```bash
mscan run-all --config configs/demo_ccle.yaml
```

For a new dataset, copy `configs/demo_ccle.yaml`, edit the paths and set
`check_reference: false`.

## Run Individual Steps

Prepare RNA only:

```bash
mscan prepare-rna --config configs/demo_ccle.yaml
```

Apply generic mutation factors after gene-matrix RNA preparation:

```bash
mscan apply-mutations --config configs/my_dataset.yaml
```

Build Vmax:

```bash
mscan build-vmax --config configs/demo_ccle.yaml
```

Run pFBA:

```bash
mscan run-pfba --config configs/demo_ccle.yaml
```

Check the bundled Figure 2 reference subset:

```bash
mscan check-reference --config configs/demo_ccle.yaml
```

## Common Overrides

Disable mutation attenuation:

```bash
mscan run-all --config configs/demo_ccle.yaml --skip-mutation
```

Change the mutation residual-capacity factor:

```bash
mscan run-all --config configs/demo_ccle.yaml --mutation-factor 0.10
```

Disable kcat or medium constraints:

```bash
mscan run-all --config configs/demo_ccle.yaml --no-kcat
mscan run-all --config configs/demo_ccle.yaml --no-media
```

Use GLPK instead of Gurobi:

```bash
mscan run-all --config configs/demo_ccle.yaml --solver glpk --skip-reference-check
```

## Manual Commands

The wrapper above calls the same individual scripts shown below.

Prepare CCLE-style RNA:

```bash
python scripts/mscan_prepare_ccle_rna.py \
  --rna data/demo/ccle_rna_depmap_3samples.csv \
  --human-gem-xlsx data/model/Human-GEM.xlsx \
  --protein-coding-genes data/annotation/protein_coding_genes_ensembl.txt \
  --mutation-table data/demo/ccle_mutation_with_prediction_final_3samples.csv \
  --include-lof \
  --mutation-factor 0.05 \
  --preprocessed-out outputs/tutorial_rna_preprocessed.csv \
  --applied-out outputs/tutorial_applied_mutations.csv \
  --out outputs/tutorial_rna_for_vmax.csv
```

Prepare a generic gene-by-sample RNA matrix:

```bash
python scripts/mscan_prepare_gene_matrix.py \
  --rna path/to/gene_by_sample_expression.csv \
  --human-gem-xlsx data/model/Human-GEM.xlsx \
  --protein-coding-genes data/annotation/protein_coding_genes_ensembl.txt \
  --input-scale linear \
  --out outputs/rna_prepared.csv
```

Apply generic mutation factors after generic RNA preparation:

```bash
python scripts/mscan_apply_mutations.py \
  --rna outputs/rna_prepared.csv \
  --mutation-table path/to/mutation_factors.csv \
  --default-factor 0.05 \
  --out outputs/rna_for_vmax.csv \
  --applied-out outputs/applied_mutations.csv
```

Build reaction-level Vmax constraints:

```bash
python scripts/mscan_build_vmax.py \
  --human-gem-xlsx data/model/Human-GEM.xlsx \
  --human-gem-mat data/model/Human-GEM.mat \
  --rna outputs/tutorial_rna_for_vmax.csv \
  --kcat data/kcat/kcat_predictions.csv \
  --fallback-kcat 1.0 \
  --out outputs/tutorial_vmax.csv
```

Run pFBA:

```bash
python scripts/mscan_run_pfba.py \
  --human-gem-mat data/model/Human-GEM.mat \
  --vmax outputs/tutorial_vmax.csv \
  --media data/media/cell_media.csv \
  --solver gurobi \
  --max-samples 3 \
  --out outputs/tutorial_pfba_results.csv
```
