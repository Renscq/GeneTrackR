# Retrieve variants from a VariantTrack or an indexed VCF file

Retrieves VCF records by optional genomic region, gene/transcript
locator, variant ID, type, or text pattern. With an in-memory
`VariantTrack`, all location arguments may be omitted to query the full
object. For large bgzip-compressed VCF files with a tabix index, a
complete genomic region is queried with
[`Rsamtools::scanTabix()`](https://rdrr.io/pkg/Rsamtools/man/scanTabix.html);
non-regional filters require reading the full file before filtering.

## Usage

``` r
retrieve_vcf(
  object,
  pattern = NULL,
  chrom = NULL,
  start = NULL,
  end = NULL,
  annotation = NULL,
  gene_id = NULL,
  transcript_id = NULL,
  upstream = 0L,
  downstream = 0L,
  strand_aware = TRUE,
  variant_id = NULL,
  variant_type = NULL,
  keep_genotype = TRUE,
  ignore_case = TRUE,
  fixed = FALSE,
  as = c("data.table", "VariantTrack"),
  verbose = TRUE,
  progress = interactive() && isTRUE(verbose)
)
```

## Arguments

- object:

  A VariantTrack object or a path to a VCF/VCF.GZ file.

- pattern:

  Optional text pattern matched against variant ID, REF, ALT, INFO, and
  variant type.

- chrom:

  Optional chromosome name or names. May be used alone as a chromosome
  filter. Required when `start` and `end` are supplied.

- start:

  Optional 1-based region start. Must be supplied together with `end`.

- end:

  Optional 1-based region end. Must be supplied together with `start`.

- annotation:

  Optional GenePred/Feature annotation object used for
  gene/transcript-aware retrieval.

- gene_id:

  Optional gene ID. When supplied, `annotation` is used to resolve the
  gene range.

- transcript_id:

  Optional transcript ID. When supplied, `annotation` is used to resolve
  the transcript range.

- upstream:

  Upstream flanking length in bp for gene/transcript queries.

- downstream:

  Downstream flanking length in bp for gene/transcript queries.

- strand_aware:

  Logical. Whether upstream/downstream should follow gene/transcript
  strand direction.

- variant_id:

  Optional variant ID vector.

- variant_type:

  Optional variant type vector, such as `SNP`, `INS`, `DEL`, or `MNV`.

- keep_genotype:

  Logical. Whether to keep FORMAT and sample genotype columns when
  `object` is a VCF file path.

- ignore_case:

  Logical. Whether pattern matching ignores case.

- fixed:

  Logical. Whether pattern is matched as a fixed string.

- as:

  Output type: `data.table` or `VariantTrack`.

- verbose:

  Logical. Whether to print progress messages.

- progress:

  Logical. Whether to print a compact stage-level progress indicator.

## Value

A data.table or VariantTrack object.

## Examples

``` r
vcf_file <- system.file("extdata", "gtr_demo_variants.vcf", package = "GeneTrackR")
vcf <- read_vcf(vcf_file, mode = "memory", verbose = FALSE)
retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200)
#> [GeneTrackR] Retrieved variants: 10.
#>      chrom      pos variant_id    ref    alt  qual filter
#>     <char>    <int>     <char> <char> <char> <num> <char>
#>  1:   chr1 12339750   varAup01      A      T    60   PASS
#>  2:   chr1 12340250     varA01      A      G    60   PASS
#>  3:   chr1 12340600     varA02      T      C    60   PASS
#>  4:   chr1 12342550     varA03      C      G    60   PASS
#>  5:   chr1 12342620    varLD01      G      A    60   PASS
#>  6:   chr1 12342710    varLD02      A      T    60   PASS
#>  7:   chr1 12342805    varLD03      C      T    60   PASS
#>  8:   chr1 12342920    varLD04      G      C    60   PASS
#>  9:   chr1 12343040    varLD05      T      G    60   PASS
#> 10:   chr1 12343180    varLD06      A      C    60   PASS
#>                                   info FORMAT    S01    S02    S03    S04
#>                                 <char> <char> <char> <char> <char> <char>
#>  1: ROLE=upstream_missing_heterozygous     GT    0/0    0/0    ./.    0/0
#>  2:        ROLE=core_haplotype,high_ld     GT    0/0    0/0    0/0    0/0
#>  3:                ROLE=core_haplotype     GT    0/0    0/0    0/0    0/0
#>  4: ROLE=core_haplotype,protein_effect     GT    0/0    0/0    0/0    0/0
#>  5:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#>  6:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#>  7:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#>  8:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#>  9:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 10:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#>        S05    S06    S07    S08    S09    S10    S11    S12    S13    S14
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    0/1    0/0    0/0    0/0    0/0    0/0    0/0    ./.    0/0    0/1
#>  2:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  3:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#>  4:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  5:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  6:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  7:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  8:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  9:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 10:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>        S15    S16    S17    S18    S19    S20    S21    S22    S23    S24
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    0/0    0/0    0/0    0/0    1/1    1/1    ./.    1/1    0/1    1/1
#>  2:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  3:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  4:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  5:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  6:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  7:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  8:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  9:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 10:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>        S25    S26    S27    S28    S29    S30    S31    S32    S33    S34
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    1/1    1/1    1/1    1/1    1/1    ./.    1/1    0/1    1/1    1/1
#>  2:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  3:    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  4:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  5:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  6:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  7:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  8:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  9:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 10:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>        S35    S36 variant_type    start      end
#>     <char> <char>       <char>    <int>    <int>
#>  1:    1/1    1/1          SNP 12339750 12339750
#>  2:    1/1    1/1          SNP 12340250 12340250
#>  3:    0/0    0/0          SNP 12340600 12340600
#>  4:    1/1    1/1          SNP 12342550 12342550
#>  5:    1/1    1/1          SNP 12342620 12342620
#>  6:    1/1    1/1          SNP 12342710 12342710
#>  7:    1/1    1/1          SNP 12342805 12342805
#>  8:    1/1    1/1          SNP 12342920 12342920
#>  9:    1/1    1/1          SNP 12343040 12343040
#> 10:    1/1    1/1          SNP 12343180 12343180
retrieve_vcf(vcf, variant_type = "SNP")
#> [GeneTrackR] Retrieved variants: 53.
#>      chrom      pos variant_id    ref    alt  qual filter
#>     <char>    <int>     <char> <char> <char> <num> <char>
#>  1:   chr1  1002500     varD01      A      C    60   PASS
#>  2:   chr1  2500000     varX01      A      G    60   PASS
#>  3:   chr1  3004500     varE01      T      G    60   PASS
#>  4:   chr1  5005200     varF01      G      A    60   PASS
#>  5:   chr1  6500000     varX02      A      G    60   PASS
#>  6:   chr1  7003300     varG01      C      T    60   PASS
#>  7:   chr1  9005000     varH01      A      G    60   PASS
#>  8:   chr1 12339750   varAup01      A      T    60   PASS
#>  9:   chr1 12340250     varA01      A      G    60   PASS
#> 10:   chr1 12340600     varA02      T      C    60   PASS
#> 11:   chr1 12342550     varA03      C      G    60   PASS
#> 12:   chr1 12342620    varLD01      G      A    60   PASS
#> 13:   chr1 12342710    varLD02      A      T    60   PASS
#> 14:   chr1 12342805    varLD03      C      T    60   PASS
#> 15:   chr1 12342920    varLD04      G      C    60   PASS
#> 16:   chr1 12343040    varLD05      T      G    60   PASS
#> 17:   chr1 12343180    varLD06      A      C    60   PASS
#> 18:   chr1 12352500    varLD07      G      T    60   PASS
#> 19:   chr1 12353100    varLD08      A      C    60   PASS
#> 20:   chr1 12353700    varLD09      C      A    60   PASS
#> 21:   chr1 12354300    varLD10      T      G    60   PASS
#> 22:   chr1 12354900    varLD11      G      C    60   PASS
#> 23:   chr1 12355500    varLD12      A      G    60   PASS
#> 24:   chr1 12356550     varB01      A      T    60   PASS
#> 25:   chr1 12359600     varB02      G      C    60   PASS
#> 26:   chr1 12362550     varB03      C      T    60   PASS
#> 27:   chr1 12365400     varB04      T      C    60   PASS
#> 28:   chr1 12369500     varC01      G      A    60   PASS
#> 29:   chr1 12372000     varC02      T      G    60   PASS
#> 30:   chr1 12374000     varC03      C      A    60   PASS
#> 31:   chr1 15002500     varI01      T      C    60   PASS
#> 32:   chr1 17500000     varX03      A      G    60   PASS
#> 33:   chr1 20004500     varJ01      G      T    60   PASS
#> 34:   chr2  2001200   varLow01      A      G    60   PASS
#> 35:   chr2  2002500   varLow02      C      T    60   PASS
#> 36:   chr2  2004200   varLow03      G      A    60   PASS
#> 37:   chr2  2006500   varLow04      T      C    60   PASS
#> 38:   chr2  3500000     varX04      A      G    60   PASS
#> 39:   chr2  4004200     varO01      A      T    60   PASS
#> 40:   chr2  6003300     varP01      C      G    60   PASS
#> 41:   chr2  7500000     varX05      A      G    60   PASS
#> 42:   chr2  8500600     varK01      A      G    60   PASS
#> 43:   chr2  8509000    varKL01      G      A    60   PASS
#> 44:   chr2  8512500    varKL02      T      C    60   PASS
#> 45:   chr2  8517000     varL01      C      T    60   PASS
#> 46:   chr2  8525000     varM01      G      C    60   PASS
#> 47:   chr2 10004500     varQ01      A      C    60   PASS
#> 48:   chr2 11000000     varX06      A      G    60   PASS
#> 49:   chr2 12004000     varR01      T      G    60   PASS
#> 50:   chr2 14002500     varS01      C      A    60   PASS
#> 51:   chr2 15500000     varX07      A      G    60   PASS
#> 52:   chr2 16996000  varPair01      G      A    60   PASS
#> 53:   chr2 17004500  varPair02      T      C    60   PASS
#>      chrom      pos variant_id    ref    alt  qual filter
#>     <char>    <int>     <char> <char> <char> <num> <char>
#>                                   info FORMAT    S01    S02    S03    S04
#>                                 <char> <char> <char> <char> <char> <char>
#>  1:            ROLE=background_variant     GT    0/0    1/1    0/0    1/1
#>  2:            ROLE=background_variant     GT    1/1    1/1    0/0    0/0
#>  3:            ROLE=background_variant     GT    0/0    0/0    1/1    0/0
#>  4:            ROLE=background_variant     GT    1/1    0/0    1/1    0/0
#>  5:            ROLE=background_variant     GT    0/0    1/1    0/0    1/1
#>  6:            ROLE=background_variant     GT    1/1    1/1    0/0    0/0
#>  7:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#>  8: ROLE=upstream_missing_heterozygous     GT    0/0    0/0    ./.    0/0
#>  9:        ROLE=core_haplotype,high_ld     GT    0/0    0/0    0/0    0/0
#> 10:                ROLE=core_haplotype     GT    0/0    0/0    0/0    0/0
#> 11: ROLE=core_haplotype,protein_effect     GT    0/0    0/0    0/0    0/0
#> 12:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 13:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 14:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 15:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 16:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 17:           ROLE=high_ld,ld_gradient     GT    0/0    0/0    0/0    0/0
#> 18:                   ROLE=ld_gradient     GT    1/1    0/0    0/0    0/0
#> 19:                   ROLE=ld_gradient     GT    1/1    1/1    0/0    0/0
#> 20:                   ROLE=ld_gradient     GT    1/1    1/1    1/1    0/0
#> 21:                   ROLE=ld_gradient     GT    1/1    1/1    1/1    1/1
#> 22:                   ROLE=ld_gradient     GT    1/1    1/1    1/1    1/1
#> 23:                   ROLE=ld_gradient     GT    1/1    1/1    1/1    1/1
#> 24:                 ROLE=geneB_variant     GT    0/0    0/0    0/0    0/0
#> 25:                 ROLE=geneB_variant     GT    1/1    1/1    1/1    1/1
#> 26:                 ROLE=geneB_variant     GT    0/0    1/1    0/0    1/1
#> 27:            ROLE=geneB_heterozygous     GT    0/0    0/0    0/0    0/0
#> 28:                 ROLE=geneC_variant     GT    1/1    0/0    1/1    0/0
#> 29:                 ROLE=geneC_variant     GT    0/0    0/0    1/1    0/0
#> 30:                 ROLE=geneC_variant     GT    1/1    1/1    0/0    0/0
#> 31:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#> 32:            ROLE=background_variant     GT    1/1    0/0    1/1    0/0
#> 33:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#> 34:                        ROLE=low_ld     GT    0/0    1/1    0/0    1/1
#> 35:                        ROLE=low_ld     GT    0/0    0/0    1/1    0/0
#> 36:                        ROLE=low_ld     GT    1/1    1/1    0/0    0/0
#> 37:                        ROLE=low_ld     GT    0/0    0/0    0/0    1/1
#> 38:            ROLE=background_variant     GT    1/1    1/1    1/1    1/1
#> 39:            ROLE=background_variant     GT    1/1    0/0    1/1    0/0
#> 40:            ROLE=background_variant     GT    0/0    0/0    1/1    0/0
#> 41:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#> 42:                 ROLE=overlap_locus     GT    0/0    0/0    0/0    0/0
#> 43:                 ROLE=overlap_locus     GT    0/0    0/0    0/0    0/0
#> 44:                 ROLE=overlap_locus     GT    0/0    1/1    0/0    1/1
#> 45:                 ROLE=overlap_locus     GT    1/1    1/1    0/0    0/0
#> 46:              ROLE=single_exon_gene     GT    0/0    0/0    0/0    0/0
#> 47:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#> 48:            ROLE=background_variant     GT    0/0    0/0    1/1    0/0
#> 49:            ROLE=background_variant     GT    0/0    0/0    0/0    0/0
#> 50:                ROLE=noncoding_gene     GT    0/0    0/0    1/1    0/0
#> 51:            ROLE=background_variant     GT    1/1    1/1    0/0    0/0
#> 52:                ROLE=two_variant_ld     GT    0/0    0/0    0/0    0/0
#> 53:                ROLE=two_variant_ld     GT    0/0    0/0    0/0    0/0
#>                                   info FORMAT    S01    S02    S03    S04
#>                                 <char> <char> <char> <char> <char> <char>
#>        S05    S06    S07    S08    S09    S10    S11    S12    S13    S14
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  2:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#>  3:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#>  4:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#>  5:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  6:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#>  7:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>  8:    0/1    0/0    0/0    0/0    0/0    0/0    0/0    ./.    0/0    0/1
#>  9:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 10:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 11:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 12:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 13:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 14:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 15:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 16:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 17:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 18:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 19:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 20:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 21:    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 22:    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 23:    1/1    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0
#> 24:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 25:    1/1    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0
#> 26:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 27:    0/1    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    0/1
#> 28:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 29:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#> 30:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 31:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 32:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 33:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 34:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 35:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#> 36:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 37:    1/1    1/1    0/0    0/0    0/0    1/1    1/1    1/1    0/0    0/0
#> 38:    1/1    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0
#> 39:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 40:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#> 41:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 42:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 43:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 44:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 45:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 46:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 47:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 48:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#> 49:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 50:    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0
#> 51:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 52:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 53:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>        S05    S06    S07    S08    S09    S10    S11    S12    S13    S14
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>        S15    S16    S17    S18    S19    S20    S21    S22    S23    S24
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  2:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#>  3:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#>  4:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#>  5:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  6:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#>  7:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>  8:    0/0    0/0    0/0    0/0    1/1    1/1    ./.    1/1    0/1    1/1
#>  9:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 10:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 11:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 12:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 13:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 14:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 15:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 16:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 17:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 18:    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1
#> 19:    0/0    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1
#> 20:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    1/1    1/1    1/1
#> 21:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    1/1    1/1
#> 22:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 23:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 24:    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0
#> 25:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 26:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 27:    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/1    0/0
#> 28:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 29:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#> 30:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#> 31:    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0
#> 32:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 33:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 34:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 35:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#> 36:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#> 37:    0/0    1/1    1/1    1/1    0/0    0/0    0/0    1/1    1/1    1/1
#> 38:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 39:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 40:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#> 41:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 42:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 43:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 44:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 45:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#> 46:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 47:    1/1    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0
#> 48:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#> 49:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 50:    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1
#> 51:    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0
#> 52:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#> 53:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1
#>        S15    S16    S17    S18    S19    S20    S21    S22    S23    S24
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>        S25    S26    S27    S28    S29    S30    S31    S32    S33    S34
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>  1:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  2:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#>  3:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#>  4:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#>  5:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#>  6:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#>  7:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>  8:    1/1    1/1    1/1    1/1    1/1    ./.    1/1    0/1    1/1    1/1
#>  9:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 10:    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 11:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 12:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 13:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 14:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 15:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 16:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 17:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 18:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 19:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 20:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 21:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 22:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 23:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 24:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 25:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 26:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 27:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    0/1    1/1    1/1
#> 28:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 29:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#> 30:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 31:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 32:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 33:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 34:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 35:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#> 36:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 37:    0/0    0/0    0/0    1/1    1/1    1/1    0/0    0/0    0/0    1/1
#> 38:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 39:    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0
#> 40:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#> 41:    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 42:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 43:    1/1    1/1    1/1    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 44:    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1    0/0    1/1
#> 45:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 46:    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 47:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 48:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#> 49:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 50:    0/0    0/0    1/1    0/0    0/0    1/1    0/0    0/0    1/1    0/0
#> 51:    1/1    1/1    0/0    0/0    1/1    1/1    0/0    0/0    1/1    1/1
#> 52:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 53:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#>        S25    S26    S27    S28    S29    S30    S31    S32    S33    S34
#>     <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#>        S35    S36 variant_type    start      end
#>     <char> <char>       <char>    <int>    <int>
#>  1:    0/0    1/1          SNP  1002500  1002500
#>  2:    0/0    0/0          SNP  2500000  2500000
#>  3:    0/0    1/1          SNP  3004500  3004500
#>  4:    1/1    0/0          SNP  5005200  5005200
#>  5:    0/0    1/1          SNP  6500000  6500000
#>  6:    0/0    0/0          SNP  7003300  7003300
#>  7:    1/1    1/1          SNP  9005000  9005000
#>  8:    1/1    1/1          SNP 12339750 12339750
#>  9:    1/1    1/1          SNP 12340250 12340250
#> 10:    0/0    0/0          SNP 12340600 12340600
#> 11:    1/1    1/1          SNP 12342550 12342550
#> 12:    1/1    1/1          SNP 12342620 12342620
#> 13:    1/1    1/1          SNP 12342710 12342710
#> 14:    1/1    1/1          SNP 12342805 12342805
#> 15:    1/1    1/1          SNP 12342920 12342920
#> 16:    1/1    1/1          SNP 12343040 12343040
#> 17:    1/1    1/1          SNP 12343180 12343180
#> 18:    1/1    1/1          SNP 12352500 12352500
#> 19:    1/1    1/1          SNP 12353100 12353100
#> 20:    1/1    1/1          SNP 12353700 12353700
#> 21:    1/1    1/1          SNP 12354300 12354300
#> 22:    1/1    1/1          SNP 12354900 12354900
#> 23:    1/1    1/1          SNP 12355500 12355500
#> 24:    0/0    0/0          SNP 12356550 12356550
#> 25:    1/1    1/1          SNP 12359600 12359600
#> 26:    0/0    1/1          SNP 12362550 12362550
#> 27:    1/1    1/1          SNP 12365400 12365400
#> 28:    1/1    0/0          SNP 12369500 12369500
#> 29:    0/0    1/1          SNP 12372000 12372000
#> 30:    0/0    0/0          SNP 12374000 12374000
#> 31:    0/0    0/0          SNP 15002500 15002500
#> 32:    1/1    0/0          SNP 17500000 17500000
#> 33:    1/1    1/1          SNP 20004500 20004500
#> 34:    0/0    1/1          SNP  2001200  2001200
#> 35:    0/0    1/1          SNP  2002500  2002500
#> 36:    0/0    0/0          SNP  2004200  2004200
#> 37:    1/1    1/1          SNP  2006500  2006500
#> 38:    1/1    1/1          SNP  3500000  3500000
#> 39:    1/1    0/0          SNP  4004200  4004200
#> 40:    0/0    1/1          SNP  6003300  6003300
#> 41:    0/0    0/0          SNP  7500000  7500000
#> 42:    1/1    1/1          SNP  8500600  8500600
#> 43:    0/0    0/0          SNP  8509000  8509000
#> 44:    0/0    1/1          SNP  8512500  8512500
#> 45:    0/0    0/0          SNP  8517000  8517000
#> 46:    1/1    1/1          SNP  8525000  8525000
#> 47:    0/0    0/0          SNP 10004500 10004500
#> 48:    0/0    1/1          SNP 11000000 11000000
#> 49:    1/1    1/1          SNP 12004000 12004000
#> 50:    0/0    1/1          SNP 14002500 14002500
#> 51:    0/0    0/0          SNP 15500000 15500000
#> 52:    1/1    1/1          SNP 16996000 16996000
#> 53:    1/1    1/1          SNP 17004500 17004500
#>        S35    S36 variant_type    start      end
#>     <char> <char>       <char>    <int>    <int>
retrieve_vcf(vcf, pattern = "varA", fixed = TRUE)
#> [GeneTrackR] Retrieved variants: 6.
#>     chrom      pos variant_id    ref    alt  qual filter
#>    <char>    <int>     <char> <char> <char> <num> <char>
#> 1:   chr1 12339750   varAup01      A      T    60   PASS
#> 2:   chr1 12340250     varA01      A      G    60   PASS
#> 3:   chr1 12340600     varA02      T      C    60   PASS
#> 4:   chr1 12342550     varA03      C      G    60   PASS
#> 5:   chr1 12344500     varA04      A    ATG    60   PASS
#> 6:   chr1 12351050     varA05    ATG      A    60   PASS
#>                                  info FORMAT    S01    S02    S03    S04    S05
#>                                <char> <char> <char> <char> <char> <char> <char>
#> 1: ROLE=upstream_missing_heterozygous     GT    0/0    0/0    ./.    0/0    0/1
#> 2:        ROLE=core_haplotype,high_ld     GT    0/0    0/0    0/0    0/0    0/0
#> 3:                ROLE=core_haplotype     GT    0/0    0/0    0/0    0/0    0/0
#> 4: ROLE=core_haplotype,protein_effect     GT    0/0    0/0    0/0    0/0    0/0
#> 5:          ROLE=core_haplotype,indel     GT    0/0    0/0    0/0    0/0    0/0
#> 6:          ROLE=core_haplotype,indel     GT    0/0    0/0    0/0    0/0    0/0
#>       S06    S07    S08    S09    S10    S11    S12    S13    S14    S15    S16
#>    <char> <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#> 1:    0/0    0/0    0/0    0/0    0/0    0/0    ./.    0/0    0/1    0/0    0/0
#> 2:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 3:    0/0    0/0    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 4:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 5:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#> 6:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>       S17    S18    S19    S20    S21    S22    S23    S24    S25    S26    S27
#>    <char> <char> <char> <char> <char> <char> <char> <char> <char> <char> <char>
#> 1:    0/0    0/0    1/1    1/1    ./.    1/1    0/1    1/1    1/1    1/1    1/1
#> 2:    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 3:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 4:    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 5:    0/0    0/0    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1
#> 6:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0
#>       S28    S29    S30    S31    S32    S33    S34    S35    S36 variant_type
#>    <char> <char> <char> <char> <char> <char> <char> <char> <char>       <char>
#> 1:    1/1    1/1    ./.    1/1    0/1    1/1    1/1    1/1    1/1          SNP
#> 2:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1          SNP
#> 3:    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0    0/0          SNP
#> 4:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1          SNP
#> 5:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1          INS
#> 6:    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1    1/1          DEL
#>       start      end
#>       <int>    <int>
#> 1: 12339750 12339750
#> 2: 12340250 12340250
#> 3: 12340600 12340600
#> 4: 12342550 12342550
#> 5: 12344500 12344500
#> 6: 12351050 12351050
vt <- retrieve_vcf(vcf, chrom = "chr1", start = 12339700, end = 12343200, as = "VariantTrack")
#> [GeneTrackR] Retrieved variants: 10.
vt
#> <VariantTrack>
#>   variants  : 10 
#>   format    : VCF 
#>   coordinate: 1-based position 
```
