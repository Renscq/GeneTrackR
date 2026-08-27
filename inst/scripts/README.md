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
- `signal_design.tsv`: transcript-level RNA-seq coverage weights, primary-transcript status, and Ribo-seq P-site count weights.
- `truth.tsv`: stable dimensions and biological/statistical expectations used for documentation and regression tests.

`generate_demo_data.R` builds all user-facing files from these canonical tables. RNA-seq signal is generated over exons as separate plus/minus bedGraph files, whereas Ribo-seq signal is generated as separate plus/minus moderately dense 1-bp heterogeneous integer P-site-like bedGraph files. Frame 0 is broadly occupied with variable heights, while frame 1 and frame 2 use different codon subsets and irregular lower counts. Total counts are designed around an 80%/10%/10% frame0/frame1/frame2 ratio, and initiation/termination frame-0 counts are approximately two times the internal mean.

Regenerate data from a source checkout with:

```sh
Rscript inst/scripts/generate_demo_data.R
Rscript inst/scripts/validate_demo_data.R
```

Do not manually edit generated `gtr_demo_*` files without updating the canonical model and generator.

Legacy `example_*` files are removed during generation. All package examples, tests, and vignettes use only the `gtr_demo_*` dataset from version 0.5.5 onward.
