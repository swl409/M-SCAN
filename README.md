# M-SCAN RNA-to-pFBA Tutorial

M-SCAN builds Human-GEM reaction-capacity constraints from RNA expression,
optional deleterious mutation attenuation, predicted kcat values, and medium
constraints, then runs pFBA. The bundled three-sample CCLE example reproduces
the Figure 2 mini-reference outputs and is also a template for applying the
workflow to new datasets.

## Quick Start

```bash
conda env create -f environment.yml
conda activate M-SCAN
pip install -e .
mscan run-all --config configs/demo_ccle.yaml
```

The same demo can be run with the shell wrapper:

```bash
bash run_tutorial.sh
```

Generated files are written to `outputs/`. The default demo runs the full
RNA-to-pFBA workflow and checks the resulting Vmax and pFBA outputs against the
stored Figure 2 three-sample reference subset.

If the environment already exists:

```bash
conda env update -n M-SCAN -f environment.yml
```

## What the Demo Runs

```text
RNA expression
  -> protein-coding and Human-GEM gene filtering
  -> target-sum normalization and zero-expression floor
  -> optional mutation attenuation
  -> GPR aggregation with or=sum and and=min
  -> optional kcat scaling
  -> medium-constrained Human-GEM pFBA
```

Default inputs:

- `data/demo/ccle_rna_depmap_3samples.csv`
- `data/demo/ccle_mutation_with_prediction_final_3samples.csv`
- `data/kcat/kcat_predictions.csv`
- `data/media/cell_media.csv`
- `data/model/Human-GEM.xlsx`
- `data/model/Human-GEM.mat`

Default outputs:

- `outputs/tutorial_rna_preprocessed.csv`
- `outputs/tutorial_rna_for_vmax.csv`
- `outputs/tutorial_applied_mutations.csv`
- `outputs/tutorial_vmax.csv`
- `outputs/tutorial_pfba_results.csv`

## Apply M-SCAN to Another Dataset

Copy the demo config and edit the paths:

```bash
cp configs/demo_ccle.yaml configs/my_dataset.yaml
```

For a generic gene-by-sample expression matrix, set:

```yaml
rna_mode: gene_matrix
rna: path/to/gene_by_sample_expression.csv
mutation_table: path/to/mutation_factors.csv
check_reference: false
```

Then run:

```bash
mscan run-all --config configs/my_dataset.yaml
```

Useful command-line overrides:

```bash
mscan run-all --config configs/my_dataset.yaml --skip-mutation
mscan run-all --config configs/my_dataset.yaml --mutation-factor 0.10
mscan run-all --config configs/my_dataset.yaml --no-kcat
mscan run-all --config configs/my_dataset.yaml --no-media
mscan run-all --config configs/my_dataset.yaml --solver glpk --skip-reference-check
```

See [docs/input_format.md](docs/input_format.md) for required columns and
[docs/step_by_step.md](docs/step_by_step.md) for manual step-by-step commands.

## Step-Wise CLI

M-SCAN can be run one step at a time:

```bash
mscan prepare-rna --config configs/demo_ccle.yaml
mscan build-vmax --config configs/demo_ccle.yaml
mscan run-pfba --config configs/demo_ccle.yaml
mscan check-reference --config configs/demo_ccle.yaml
```

For generic gene-matrix mode with a separate sample-gene mutation table:

```bash
mscan prepare-rna --config configs/my_dataset.yaml
mscan apply-mutations --config configs/my_dataset.yaml
mscan build-vmax --config configs/my_dataset.yaml
mscan run-pfba --config configs/my_dataset.yaml
```

## Configuration

The config file is a simple `key: value` file read without extra YAML
dependencies. Important keys include:

- `rna_mode`: `ccle` or `gene_matrix`.
- `rna`: RNA input table.
- `rna_layout`: CCLE-mode RNA layout, either `auto`, `sample_rows`, or
  `gene_rows`.
- `sample_column`: CCLE-mode sample/cell-line column. Defaults to
  `cell_line_display_name`; rowname-style `Unnamed: 0` input is also accepted.
- `mutation_table`: optional mutation table.
- `mutation_factor`: residual capacity assigned to mutation-affected genes.
- `use_kcat`: whether to use `data/kcat/kcat_predictions.csv`.
- `use_media`: whether to apply the exchange medium table.
- `solver`: `gurobi`, `glpk`, or `auto`.
- `check_reference`: should be `true` only for the bundled three-sample demo.

M-SCAN does not run LoGoFunc internally. If LoGoFunc missense annotations are
available, provide them in the `final_prediction` column of the mutation table.
Rows with `final_prediction == LOF` are treated as deleterious when
`include_lof: true`.

## GPR Rule

M-SCAN treats isoenzymes connected by `or` as additive and protein complexes
connected by `and` as limited by the minimum subunit value. Mixed rules are
parsed recursively.

Examples:

- `A or B`: `A + B`
- `A and B`: `min(A, B)`
- `(A and B) or C`: `min(A, B) + C`
- `(A or B) and C`: `min(A + B, C)`

## Medium

The default medium is the Ham's F-12-derived Figure 2 benchmark medium mapped to
Human1/Human-GEM exchange reactions. It includes mapped Ham's F-12 components,
core exchange reactions (`H2O`, `O2`, `H+`), minimal feasibility
add-backs. Replace `data/media/cell_media.csv` or set `use_media:
false` to test a different medium definition.

## Notes

The bundled data are intentionally small so the repository can be tested
quickly on another machine. For full-scale runs, increase `max_samples` or set
it to `0` to process all sample columns in the Vmax table.

## License

The tutorial code is released under the MIT License. Included example data and
third-party model/resource files remain subject to their original provider
licenses and terms of use.
