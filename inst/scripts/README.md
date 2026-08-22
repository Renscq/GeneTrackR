# GeneTrackR demo-data source model

The `demo_model/` directory is the canonical source for the GeneTrackR example dataset.
Generated user-facing files are written to `inst/extdata/` with the `gtr_demo_` prefix.

Files in `demo_model/` use 1-based closed coordinates. Format-specific coordinate conversion is performed only when output files are generated.

- `chromosomes.tsv`: chromosome names and lengths.
- `genes.tsv`: expected gene-level intervals derived from the transcript model.
- `transcripts.tsv`: canonical transcript, exon, and CDS definitions.
- `samples.tsv`: 36 sample IDs and four balanced GeneA haplotype groups.
- `variants.tsv`: 50 designed variants, functional demo roles, and deterministic genotype patterns.
- `features.tsv`: regulatory and interval feature definitions in canonical coordinates.
- `signal_design.tsv`: transcript-level RNA-seq coverage weights, primary-transcript status, and Ribo-seq translation-density weights.
- `truth.tsv`: stable dimensions and biological/statistical expectations used for documentation and regression tests.

`generate_demo_data.R` builds all user-facing files from these canonical tables. RNA-seq signal is generated over exons as separate plus/minus bedGraph files, whereas Ribo-seq signal is generated as separate plus/minus one-base CDS-density bedGraph files with transcript-oriented three-nucleotide periodicity and initiation/termination peaks.

Regenerate data from a source checkout with:

```sh
Rscript inst/scripts/generate_demo_data.R
Rscript inst/scripts/validate_demo_data.R
```

Do not manually edit generated `gtr_demo_*` files without updating the canonical model and generator.
