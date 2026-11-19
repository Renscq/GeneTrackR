test_that("edited factor levels survive feature-level export and subsets", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  names(ranges) <- c("late", "early")
  ranges$group <- factor(c("late", "early"), levels = c("early", "late", "unused"))
  feature <- as_feature(ranges)
  edited <- factor(c("novel", "late"), levels = c("early", "late", "novel", "unused"))
  feature$data$group <- edited

  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  expect_identical(S4Vectors::mcols(observed)$group, edited)
  expect_identical(names(observed), c("early", "late"))
  selected <- retrieve_feature(feature, chrom = "chr1", start = 1L, end = 10L)
  expect_no_warning(selected_ranges <- as_granges(selected, level = "feature"))
  expect_identical(S4Vectors::mcols(selected_ranges)$group, edited[1L])
  empty <- retrieve_feature(feature, chrom = "chr1", start = 50L, end = 60L)
  expect_identical(S4Vectors::mcols(as_granges(empty, level = "feature"))$group,
                   edited[integer()])
})

test_that("IntegerList edits cannot silently truncate or change missing values", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  ranges$members <- IRanges::IntegerList(20L, 5L)
  feature <- as_feature(ranges)
  for (invalid in list(c(1.25, 2.75), 2^40, NaN, Inf, "01")) {
    data.table::set(feature$data, j = "members", value = list(list(invalid, 3L)))
    before <- serialize(feature, NULL)
    expect_no_warning(expect_error(
      as_granges(feature, level = "feature"), "members.*without loss"
    ))
    expect_identical(serialize(feature, NULL), before)
  }
})

test_that("exact numeric IntegerList edits keep annotations in nonempty and empty subsets", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  members <- IRanges::IntegerList(late = 20L, early = 5L)
  S4Vectors::mcols(members)$origin <- c("late_origin", "early_origin")
  S4Vectors::metadata(members) <- list(method = "observed")
  ranges$members <- members
  feature <- as_feature(ranges)
  data.table::set(feature$data, j = "members", value = list(list(c(7, NA_real_), numeric())))

  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  observed_members <- S4Vectors::mcols(observed)$members
  expect_identical(as.list(observed_members), list(early = c(7L, NA_integer_), late = integer()))
  expect_identical(S4Vectors::mcols(observed_members)$origin, c("early_origin", "late_origin"))
  expect_identical(S4Vectors::metadata(observed_members), list(method = "observed"))
  selected <- retrieve_feature(feature, chrom = "chr1", start = 1L, end = 10L)
  selected_members <- S4Vectors::mcols(as_granges(selected, level = "feature"))$members
  expect_identical(as.list(selected_members), list(early = c(7L, NA_integer_)))
  empty <- retrieve_feature(feature, chrom = "chr1", start = 50L, end = 60L)
  empty_members <- S4Vectors::mcols(as_granges(empty, level = "feature"))$members
  expect_length(empty_members, 0L)
  expect_s4_class(empty_members, "IntegerList")
  expect_identical(S4Vectors::metadata(empty_members), list(method = "observed"))
})

test_that("unannotated numeric Rle fields export normalized missing-value types", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  scores <- list(S4Vectors::Rle(c("4.5", NA_character_)),
                 S4Vectors::Rle(factor(c("4.5", NA_character_))))
  for (score in scores) {
    ranges$score <- score
    ranges$exon_number <- S4Vectors::Rle(c("2", NA_character_))
    feature <- as_feature(ranges)
    expect_no_warning(observed <- as_granges(feature, level = "feature"))
    expect_identical(S4Vectors::decode(S4Vectors::mcols(observed)$score), c(NA_real_, 4.5))
    expect_identical(S4Vectors::decode(S4Vectors::mcols(observed)$exon_number), c(NA_integer_, 2L))
  }
})

test_that("standard Rle fields retain containers with normalized value types and annotations", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  score <- S4Vectors::Rle(c("4.5", "2.5"))
  S4Vectors::mcols(score)$origin <- c("late_origin", "early_origin")
  S4Vectors::metadata(score) <- list(scale = "numeric_text")
  ranges$score <- score
  ranges$exon_number <- S4Vectors::Rle(factor(c("2", "1")))
  ranges$feature_id <- S4Vectors::Rle(factor(c("b", "a")))
  feature <- as_feature(ranges)

  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  observed_score <- S4Vectors::mcols(observed)$score
  expect_s4_class(observed_score, "Rle")
  expect_identical(unname(S4Vectors::decode(observed_score)), c(2.5, 4.5))
  expect_null(names(observed_score))
  expect_identical(S4Vectors::mcols(observed_score)$origin, c("early_origin", "late_origin"))
  expect_identical(S4Vectors::metadata(observed_score), list(scale = "numeric_text"))
  expect_identical(S4Vectors::decode(S4Vectors::mcols(observed)$exon_number), c(1L, 2L))
  expect_identical(S4Vectors::decode(S4Vectors::mcols(observed)$feature_id), c("a", "b"))
  empty <- retrieve_feature(feature, chrom = "chr1", start = 50L, end = 60L)
  empty_score <- S4Vectors::mcols(as_granges(empty, level = "feature"))$score
  expect_identical(unname(S4Vectors::decode(empty_score)), numeric())
  expect_identical(S4Vectors::metadata(empty_score), list(scale = "numeric_text"))
})

test_that("additional factor Rle edits retain new levels and aligned annotations", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  category <- S4Vectors::Rle(factor(c("late", "early"), levels = c("early", "late")))
  S4Vectors::mcols(category)$origin <- c("late_origin", "early_origin")
  S4Vectors::metadata(category) <- list(method = "observed")
  ranges$category <- category
  feature <- as_feature(ranges)
  edited <- factor(c("novel", "late"), levels = c("early", "late", "novel", "unused"))
  feature$data$category <- edited

  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  observed_category <- S4Vectors::mcols(observed)$category
  expect_s4_class(observed_category, "Rle")
  expect_identical(S4Vectors::decode(observed_category), edited)
  expect_identical(S4Vectors::mcols(observed_category)$origin, c("early_origin", "late_origin"))
  expect_identical(S4Vectors::metadata(observed_category), list(method = "observed"))
})

test_that("typed logical List edits cannot silently become boolean values", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3L))
  ranges$flags <- IRanges::LogicalList(TRUE)
  feature <- as_feature(ranges)
  feature$data$flags[[1L]] <- c(0.5, 1.5)
  expect_no_warning(expect_error(
    as_granges(feature, level = "feature"), "flags.*without loss"
  ))
})

test_that("exact typed List edits retain element names", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1L, 3L))
  ranges$members <- IRanges::IntegerList(c(original = 1L))
  feature <- as_feature(ranges)
  for (values in list(c(first = 7L, second = 8L), c(first = 7, second = 8))) {
    feature$data$members[[1L]] <- values
    expect_no_warning(observed <- as_granges(feature, level = "feature"))
    expect_identical(as.list(S4Vectors::mcols(observed)$members),
                     list(c(first = 7L, second = 8L)))
  }
})

test_that("typed List edits can replace or remove original element names", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  ranges$members <- IRanges::IntegerList(late = c(z = 1L), early = c(y = 2L))
  feature <- as_feature(ranges)
  edits <- list(list(3L, 4L), list(c(revised = 7L), c(other = 8L)))
  for (values in edits) {
    data.table::set(feature$data, j = "members", value = list(values))
    expect_no_warning(observed <- as_granges(feature, level = "feature"))
    expect_identical(as.list(S4Vectors::mcols(observed)$members),
                     stats::setNames(values, c("early", "late")))
  }
})

test_that("FactorList edits can introduce and reorder shared levels", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  ranges$category <- IRanges::FactorList(factor("late", levels = c("early", "late")),
                                       factor("early", levels = c("early", "late")))
  feature <- as_feature(ranges)
  values <- list(factor("novel", levels = c("novel", "late", "early", "unused")),
                 factor("late", levels = c("novel", "late", "early", "unused")))
  data.table::set(feature$data, j = "category", value = list(values))
  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  expect_identical(class(S4Vectors::mcols(observed)$category), class(ranges$category))
  expect_identical(as.list(S4Vectors::mcols(observed)$category), values)
})

test_that("empty FactorList metadata retains shared levels and annotations", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  category <- IRanges::FactorList(
    factor("late", levels = c("early", "late", "unused")),
    factor("early", levels = c("early", "late", "unused"))
  )
  S4Vectors::mcols(category)$origin <- c("late_origin", "early_origin")
  S4Vectors::metadata(category) <- list(method = "observed")
  ranges$category <- category
  feature <- as_feature(ranges)
  empty_features <- list(
    retrieve_feature(feature, chrom = "chr1", start = 50L, end = 60L),
    as_feature(ranges[integer()])
  )
  for (empty in empty_features) {
    expect_no_warning(observed <- as_granges(empty, level = "feature"))
    observed_category <- S4Vectors::mcols(observed)$category
    expect_identical(levels(unlist(observed_category)), c("early", "late", "unused"))
    expect_identical(observed_category, category[integer()])
  }
})

test_that("classed list metadata round trips without recursive indexing", {
  ranges <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(20L, 5L), c(25L, 9L)))
  values <- list(as.POSIXlt("2020-01-02", tz = "UTC"), as.POSIXlt("2020-01-01", tz = "UTC"))
  ranges$observed_at <- S4Vectors::SimpleList(values)
  feature <- as_feature(ranges)
  expect_no_warning(observed <- as_granges(feature, level = "feature"))
  expect_identical(class(S4Vectors::mcols(observed)$observed_at), class(ranges$observed_at))
  expect_identical(as.list(S4Vectors::mcols(observed)$observed_at), values[c(2L, 1L)])
})
