# Build haplotypes from variants in a genomic region

Extracts variants from a user-defined genomic interval and converts
sample genotype profiles into haplotypes.

## Usage

``` r
hap_region_variant(
  vcf,
  chrom,
  start,
  end,
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

- chrom:

  Chromosome name.

- start:

  Region start in 1-based closed coordinates.

- end:

  Region end in 1-based closed coordinates.

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
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
hap <- hap_region_variant(vcf, chrom = "chr1", start = 12339700, end = 12352000, genotype_mode = "code", min_variant_number = 1)
#> [GeneTrackR] Retrieved variants: 12.
hap
#> <HapVariant>
#>   region    : chr1:12339700-12352000
#>   variants  : 12
#>   samples   : 36
#>   haplotypes: 10
hap$haplotypes
#>     hap_id sample_n                         samples varAup01 varA01 varA02
#>     <char>    <int>                          <char>   <char> <char> <char>
#>  1:   Hap1        8 S28;S29;S31;S32;S33;S34;S35;S36        1      1      0
#>  2:   Hap2        8 S19;S20;S22;S23;S24;S25;S26;S27        1      1      1
#>  3:   Hap3        7     S01;S02;S04;S06;S07;S08;S09        0      0      0
#>  4:   Hap4        7     S10;S11;S13;S15;S16;S17;S18        0      0      1
#>  5:  Hap10        1                             S21     <NA>      1      1
#>  6:   Hap5        1                             S05        1      0      0
#>  7:   Hap6        1                             S14        1      0      1
#>  8:   Hap7        1                             S03     <NA>      0      0
#>  9:   Hap8        1                             S12     <NA>      0      1
#> 10:   Hap9        1                             S30     <NA>      1      0
#>     varA03 varLD01 varLD02 varLD03 varLD04 varLD05 varLD06 varA04 varA05
#>     <char>  <char>  <char>  <char>  <char>  <char>  <char> <char> <char>
#>  1:      1       1       1       1       1       1       1      1      1
#>  2:      1       1       1       1       1       1       1      1      0
#>  3:      0       0       0       0       0       0       0      0      0
#>  4:      0       0       0       0       0       0       0      0      0
#>  5:      1       1       1       1       1       1       1      1      0
#>  6:      0       0       0       0       0       0       0      0      0
#>  7:      0       0       0       0       0       0       0      0      0
#>  8:      0       0       0       0       0       0       0      0      0
#>  9:      0       0       0       0       0       0       0      0      0
#> 10:      1       1       1       1       1       1       1      1      1
```
