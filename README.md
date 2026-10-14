# LCDM_AnchorItems_Simulation
This is a accompanying repository for simulation code related to "Anchor Items in Cognitive Diagnostic Models: Towards Best Practices for Operational Testing".

`Sim_DataGen.R` generates each dataset once. The `analysis_conditions.csv`
manifest then maps every dataset to the four intended analyses: true or
misspecified Q-matrix crossed with multigroup or pooled calibration. This makes
the analysis comparisons use identical simulated observations.

## Test run

```sh
Rscript Sim_DataGen.R \
  --out-dir=generated_data_test \
  --n-reps=2 \
  --condition-limit=2 \
  --workers=1
```

## Full run

Use a new, empty output directory. Start with one or two workers and increase
only after monitoring memory and disk use.

```sh
Rscript Sim_DataGen.R \
  --out-dir=generated_data_full \
  --n-reps=100 \
  --workers=2
```

The full design contains 90 data conditions and 9,000 generated replications.
`analysis_conditions.csv` contains 360 analysis conditions. Zero-anchor
conditions use `anchor_location = "none"`; single and mixed placement levels
are crossed only with the four- and eight-anchor conditions.

The generator refuses to write into a nonempty output directory unless
`--overwrite` is supplied. Replication seeds are deterministic functions of
the master seed, data-condition index, and replication number, so results do
not depend on worker scheduling. Multiple workers use base-R process
parallelism on macOS and Linux. Generation falls back to one worker on Windows;
this affects runtime but not the generated values.

## Build flexMIRT files

After generation, build all four analysis files for every available
replication:

```sh
Rscript WriteFlexmirtFiles.R --data-dir=generated_data
```

The default settings applied to every file are `MaxE=20000`, `MaxM=100`,
`Etol=1e-4`, `Mtol=1e-7`, and `Processors=1`. Change them with the corresponding
command-line options. The writer also creates `responses_all.csv` for pooled
calibration, using `-99` for structural missing responses, and writes a
validated `flexmirt_manifest.csv`.

The four source flexMIRT templates are stored in this directory with the R
scripts. `WriteFlexmirtFiles.R` uses this directory as its default template
location, so no path back to the earlier workflow is required.

## Fit the same models with GDINA

`RunGDINA.R` fits the LCDM using `GDINA::GDINA(..., model="logitGDINA")`.
It recreates the flexMIRT equality structure by merging grade-item columns only
for the items identified as anchors. Multigroup fits use grade-specific latent
class distributions; pooled fits use one common distribution.

Run all available replications:

```sh
Rscript RunGDINA.R --data-dir=generated_data
```

For a one-fit smoke test:

```sh
Rscript RunGDINA.R --data-dir=generated_data --analysis-limit=1
```

Results are checkpointed after every fit in `gdina_convergence.csv`. Existing
rows are skipped unless `--overwrite` is supplied. Fitted objects are not saved
unless `--save-fits` is requested.

## Run flexMIRT through irtQ

On the Windows computer where flexMIRT is installed, validate the complete
batch without launching the software:

```sh
Rscript RunFlexmirt.R --data-dir=generated_data --dry-run
```

Then run the full manifest:

```sh
Rscript RunFlexmirt.R \
  --data-dir=generated_data \
  --flex-dir="C:/Program Files/flexMIRT 3.6"
```

If flexMIRT is installed in a standard location, `--flex-dir` may be omitted
and `irtQ` will search under `C:/Program Files`. The runner processes only the
validated files in `flexmirt_manifest.csv`, checkpoints after every model, and
resumes by skipping recorded runs. It records first- and second-order status
separately and does not mistake phrases such as "standard error" for fatal
errors.
