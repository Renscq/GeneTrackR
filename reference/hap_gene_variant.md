# Build haplotypes from variants in a gene or transcript region

Extracts variants for a specified gene or transcript, optionally
including upstream and downstream flanking regions, and converts sample
genotype profiles into haplotypes.

## Usage

``` r
hap_gene_variant(
  vcf,
  annotation,
  gene_id = NULL,
  transcript_id = NULL,
  upstream = 0L,
  downstream = 0L,
  strand_aware = TRUE,
  samples = NULL,
  variant_type = NULL,
  genotype_mode = c("code", "string"),
  missing_genotype = NA_character_,
  min_variant_number = NULL
)
```

## Arguments

- vcf:

  A VariantTrack object or VCF file path.

- annotation:

  A gene annotation object used to locate `gene_id` or `transcript_id`.

- gene_id:

  Optional gene ID. Use exactly one of `gene_id` or `transcript_id`.

- transcript_id:

  Optional transcript ID. Use exactly one of `gene_id` or
  `transcript_id`.

- upstream:

  Upstream flanking length in bp. When `strand_aware = TRUE`, upstream
  is interpreted relative to the gene/transcript strand.

- downstream:

  Downstream flanking length in bp. When `strand_aware = TRUE`,
  downstream is interpreted relative to the gene/transcript strand.

- strand_aware:

  Logical. Whether upstream/downstream should follow strand direction.
  Default TRUE.

- samples:

  Optional sample names to keep.

- variant_type:

  Optional variant types to keep.

- genotype_mode:

  Genotype representation. `code` converts genotypes to compact 0/1
  states, where 0 means reference genotype and 1 means any alternate
  allele is present. `string` converts genotypes to a single allele
  label; long InDel alleles are compressed as `iN`, where `N` is the
  displayed REF or ALT allele length, so the haplotype cell uses the
  same label as the corresponding REF/ALT row.

- missing_genotype:

  Missing genotype label. Default is `NA_character_`, which is displayed
  as `NA` in haplotype tables.

- min_variant_number:

  Minimum number of non-missing variants required for a sample. If NULL,
  only samples with complete non-missing genotypes across all retained
  variants are kept.

## Value

A HapVariant object.

## Examples

``` r
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
anno_file <- system.file("extdata", "gtr_demo.genePredExt", package = "GeneTrackR")
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
anno <- read_genepred(anno_file, format = "genePredExt", verbose = FALSE)
hap <- hap_gene_variant(vcf, annotation = anno, gene_id = "GeneA", upstream = 1000, downstream = 500, min_variant_number = 1)
#> [GeneTrackR] Retrieved variants: 13.
hap
#> <HapVariant>
#>   region    : chr1:12339001-12352500
#>   variants  : 13
#>   samples   : 36
#>   haplotypes: 12
hap$haplotypes
#>     hap_id sample_n                         samples varAup01 varA01 varA02
#>     <char>    <int>                          <char>   <char> <char> <char>
#>  1:   Hap1        8 S28;S29;S31;S32;S33;S34;S35;S36        1      1      0
#>  2:   Hap2        7     S10;S11;S13;S15;S16;S17;S18        0      0      1
#>  3:   Hap3        7     S20;S22;S23;S24;S25;S26;S27        1      1      1
#>  4:   Hap4        6         S02;S04;S06;S07;S08;S09        0      0      0
#>  5:  Hap10        1                             S12     <NA>      0      1
#>  6:  Hap11        1                             S30     <NA>      1      0
#>  7:  Hap12        1                             S21     <NA>      1      1
#>  8:   Hap5        1                             S01        0      0      0
#>  9:   Hap6        1                             S05        1      0      0
#> 10:   Hap7        1                             S14        1      0      1
#> 11:   Hap8        1                             S19        1      1      1
#> 12:   Hap9        1                             S03     <NA>      0      0
#>     varA03 varLD01 varLD02 varLD03 varLD04 varLD05 varLD06 varA04 varA05
#>     <char>  <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#>  1:      1       1       1       1       1       1       1      1      1
#>  2:      0       0       0       0       0       0       0      0      0
#>  3:      1       1       1       1       1       1       1      1      0
#>  4:      0       0       0       0       0       0       0      0      0
#>  5:      0       0       0       0       0       0       0      0      0
#>  6:      1       1       1       1       1       1       1      1      1
#>  7:      1       1       1       1       1       1       1      1      0
#>  8:      0       0       0       0       0       0       0      0      0
#>  9:      0       0       0       0       0       0       0      0      0
#> 10:      0       0       0       0       0       0       0      0      0
#> 11:      1       1       1       1       1       1       1      1      0
#> 12:      0       0       0       0       0       0       0      0      0
#>     varLD07
#>      <char>
#>  1:       1
#>  2:       0
#>  3:       1
#>  4:       0
#>  5:       0
#>  6:       1
#>  7:       1
#>  8:       1
#>  9:       0
#> 10:       0
#> 11:       0
#> 12:       0
hap_tx <- hap_gene_variant(vcf, annotation = anno, transcript_id = "TxA1", genotype_mode = "string", min_variant_number = 1)
#> [GeneTrackR] Retrieved variants: 11.
hap_tx$haplotypes
#> Key: <hap_id>
#>    hap_id sample_n                             samples varA01 varA02 varA03
#>    <char>    <int>                              <char> <char> <char> <char>
#> 1:   Hap1        9 S10;S11;S12;S13;S14;S15;S16;S17;S18      A      C      C
#> 2:   Hap2        9 S01;S02;S03;S04;S05;S06;S07;S08;S09      A      T      C
#> 3:   Hap3        9 S19;S20;S21;S22;S23;S24;S25;S26;S27      G      C      G
#> 4:   Hap4        9 S28;S29;S30;S31;S32;S33;S34;S35;S36      G      T      G
#>    varLD01 varLD02 varLD03 varLD04 varLD05 varLD06 varA04 varA05
#>     <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#> 1:       G       A       C       G       T       A      A     i3
#> 2:       G       A       C       G       T       A      A     i3
#> 3:       A       T       T       C       G       C     i3     i3
#> 4:       A       T       T       C       G       C     i3      A
```
