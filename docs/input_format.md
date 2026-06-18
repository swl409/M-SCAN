# Input Formats

M-SCAN has two RNA input modes. Use `rna_mode: ccle` for the bundled DepMap/CCLE
layout and `rna_mode: gene_matrix` for a generic gene-by-sample expression
matrix.

## CCLE Mode

This mode is used by `configs/demo_ccle.yaml` and reproduces the three-sample
Figure 2 smoke test.

CCLE mode accepts either sample-row or gene-row RNA tables. It still interprets
expression values as DepMap-style `log2(TPM+1)`.

Sample-row layout, matching the bundled DepMap/CCLE demo:

```text
cell_line_display_name, ENSG000001..., ENSG000002..., ...
HOP92, 4.21, 3.02, ...
HS578T, 5.11, 0.00, ...
```

Rowname-style CSV files are also accepted. For example, a file saved with cell
line names as row names may be read back with an `Unnamed: 0` sample column:

```text
Unnamed: 0, ENSG000001..., ENSG000002..., ...
HOP92, 4.21, 3.02, ...
HS578T, 5.11, 0.00, ...
```

If your sample column has another name, set `sample_column` in the config or
pass `--sample-column`. The bundled demo keeps the original DepMap-style
`cell_line_display_name` column for reference-check compatibility.

Gene-row layout is also accepted directly:

```text
gene_id,HOP92,HS578T,MDAMB231
ENSG00000141510,4.21,5.11,3.02
ENSG00000171862,3.02,0.00,2.44
```

or as a rowname-style CSV:

```text
Unnamed: 0,HOP92,HS578T,MDAMB231
ENSG00000141510,4.21,5.11,3.02
ENSG00000171862,3.02,0.00,2.44
```

Use `rna_layout: auto`, `sample_rows`, or `gene_rows` in the config if you want
to force a layout. M-SCAN converts DepMap-style `log2(TPM+1)` values to linear
TPM-like values, keeps protein-coding Human-GEM genes, normalizes each sample
to a target sum of 1e6, replaces zero values with the sample-wise minimum
positive value, and writes a standardized `gene_id x sample` table for Vmax
construction.

Optional CCLE mutation table columns:

```text
cell_line_display_name
Gene
ensembl_gene_id
variant_info
Consequence
final_prediction
```

Protein-disrupting VEP-style terms are attenuated by the configured residual
capacity factor. Missense variants are attenuated only when
`final_prediction == LOF` and `include_lof: true`.

M-SCAN does not run LoGoFunc internally. If LoGoFunc missense annotations are
available, provide them in the `final_prediction` column. Rows with
`final_prediction == LOF` are treated as deleterious when `include_lof: true`.

## Generic Gene-Matrix Mode

Use this mode for new datasets that are already in gene-by-sample form.

Required RNA layout:

```text
gene_id,SAMPLE_1,SAMPLE_2,SAMPLE_3
ENSG00000141510,10.2,0.0,4.5
ENSG00000171862,3.1,2.2,0.0
```

This is the recommended general input form: rows are genes and columns are
samples. The first column is treated as the Ensembl gene ID column, regardless
of its name, so rowname-style CSV exports such as `Unnamed: 0,SAMPLE_1,SAMPLE_2`
are accepted. Sample columns must be numeric.

After preparation, M-SCAN writes the standardized Vmax input as:

```text
gene_id,SAMPLE_1,SAMPLE_2,SAMPLE_3
ENSG00000141510,normalized_value,normalized_value,normalized_value
```

The Vmax builder also accepts this standardized table or a rowname-style first
column containing gene IDs.

Supported expression scales:

- `input_scale: log2`: values are converted from log2(TPM+1)-like scale.
- `input_scale: linear`: values are treated as TPM-like non-negative values.
- `input_scale: already_normalized`: values are used without target-sum normalization.

Optional generic mutation table columns:

```text
sample,gene_id,mutation_factor
SAMPLE_1,ENSG00000141510,0.05
SAMPLE_2,ENSG00000171862,0.10
```

`mutation_factor` is optional. If absent, `mutation_factor` from the config is
used for all listed sample-gene pairs.

## kcat Table

Required columns:

```text
reaction_id,gene_id,predicted_kcat_s_inv
MAR00001,ENSG000001...,12.5
```

If `use_kcat: false`, reaction capacities are built from RNA-derived gene
activity only.

## Medium Table

Required column:

```text
HMR
MAR09034
MAR09048
```

Additional columns such as `metabolite`, `component`, and `mapping_confidence`
are allowed and are kept for readability. If `use_media: true`, exchange
reactions listed in `HMR` are opened for uptake; all other exchange uptake
bounds are closed.

## Human-GEM Files

The workflow requires:

- `Human-GEM.xlsx`: reaction and GPR annotation.
- `Human-GEM.mat`: COBRA model loaded by COBRApy.

The bundled demo files are included only to make the tutorial immediately
runnable.
