annotation_parent_fixture <- function() {
    file <- tempfile(fileext = ".gff3")
    writeLines(c("##gff-version 3",
        "chr1\tfixture\tgene\t1\t20\t.\t+\t.\tID=g1;Name=Gene1",
        "chr1\tfixture\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\tfixture\tmRNA\t1\t20\t.\t+\t.\tID=t2;Parent=g1",
        "chr1\tfixture\texon\t1\t5\t.\t+\t.\tID=e1;Parent=t1,t2"), file)
    file
}

test_that("GFF3 preserves complete parents without expanding the primary model", {
    file <- annotation_parent_fixture()
    on.exit(unlink(file), add = TRUE)
    feature <- read_gff(file, verbose = FALSE, progress = FALSE)
    exon <- which(feature$data$feature_id == "e1")
    expect_equal(nrow(feature$data), 4L)
    expect_identical(feature$data$parent_id[exon], "t1")
    expect_identical(feature$data$transcript_id[exon], "t1")
    expect_identical(feature$data$parent_ids[[exon]], c("t1", "t2"))
    expect_identical(feature$data$parent_ids[[which(feature$data$feature_id == "g1")]],
        character())
    expect_identical(feature$exons$transcript_id, "t1")
    expect_match(feature$data$attribute[exon], "Parent=t1,t2", fixed = TRUE)
    ranges <- as_granges(feature, level = "feature")
    restored <- as_feature(ranges)
    expect_identical(restored$data$parent_ids, feature$data$parent_ids)
    selected <- as_feature(ranges[S4Vectors::mcols(ranges)$feature_id == "e1"])
    expect_identical(selected$data$parent_ids[[1L]], c("t1", "t2"))
    expect_error(retrieve_feature(feature, type = "exon"), "parent.*GRanges|GRanges.*parent")
    empty <- feature
    empty$data <- data.table::copy(feature$data[0L])
    expect_length(S4Vectors::mcols(as_granges(empty, level = "feature"))$parent_ids, 0L)
})

test_that("regenerated GFF3 attributes retain parents and reject conflicting edits", {
    file <- annotation_parent_fixture()
    output <- tempfile(fileext = ".gff3")
    on.exit(unlink(c(file, output)), add = TRUE)
    feature <- read_gff(file, verbose = FALSE)
    before <- serialize(feature, NULL)
    write_feature(feature, output, format = "gff", keep_attributes = FALSE)
    restored <- read_gff(output, verbose = FALSE)
    exon <- which(restored$data$feature_id == "e1")
    expect_identical(restored$data$parent_ids[[exon]], c("t1", "t2"))
    expect_identical(serialize(feature, NULL), before)
    edited <- feature
    edited$data <- data.table::copy(feature$data)
    exon <- which(edited$data$feature_id == "e1")
    edited$data$parent_ids[[exon]] <- c("t1", "t2", "t3")
    write_feature(edited, output, format = "gff", keep_attributes = FALSE,
        overwrite = TRUE)
    expect_true(any(grepl("Parent=t1,t2,t3", readLines(output), fixed = TRUE)))
    edited$data$parent_id[exon] <- "t9"
    expect_error(write_feature(edited, tempfile(fileext = ".gff3"),
        keep_attributes = FALSE), "parent_id.*parent_ids|parent_ids.*parent_id")
    legacy <- feature
    legacy$data <- data.table::copy(feature$data)
    legacy$data$parent_ids <- NULL
    write_feature(legacy, output, format = "gff", keep_attributes = FALSE,
        overwrite = TRUE)
    expect_false(any(grepl("Parent=t1,t2", readLines(output), fixed = TRUE)))
})

test_that("GFF parent tokens retain percent escapes and reject malformed list edits", {
    file <- tempfile(fileext = ".gff3")
    on.exit(unlink(file), add = TRUE)
    writeLines(c("##gff-version 3",
        "chr1\ts\tregion\t1\t5\t.\t+\t.\tID=e1;Parent=t%2C1,t2"), file)
    feature <- read_gff(file, verbose = FALSE)
    expect_identical(feature$data$parent_ids[[1L]], c("t%2C1", "t2"))
    output <- tempfile(fileext = ".gff3")
    on.exit(unlink(output), add = TRUE)
    write_feature(feature, output, keep_attributes = FALSE)
    expect_match(readLines(output), "Parent=t%2C1,t2", fixed = TRUE)
    feature$data$parent_ids[[1L]] <- c("t%2C1", NA_character_)
    expect_error(write_feature(feature, tempfile(fileext = ".gff3"),
        keep_attributes = FALSE), "parent_ids")
})

annotation_bed12_fixture <- function() {
    file <- tempfile(fileext = ".bed")
    writeLines(c(
        "chr2\t9\t20\tsame\t1000\t-\t9\t9\t0\t1\t11,\t0,",
        "chr1\t0\t20\tsame\t0\t+\t1\t18\t255,0,128\t2\t5,11,\t0,9,"), file)
    file
}

test_that("BED12 metadata retains blocks thick intervals and RGB in flat Features", {
    file <- annotation_bed12_fixture()
    on.exit(unlink(file), add = TRUE)
    feature <- read_bed(file, verbose = FALSE, progress = FALSE)
    expect_false(inherits(feature, "GenePred"))
    expect_equal(nrow(feature$data), 2L)
    expect_identical(feature$data$chrom, c("chr1", "chr2"))
    expect_identical(feature$data$block_count, c(2L, 1L))
    expect_identical(feature$data$block_sizes, list(c(5L, 11L), 11L))
    expect_identical(feature$data$block_starts, list(c(0L, 9L), 0L))
    expect_identical(feature$data$thick_start, c(2L, 10L))
    expect_identical(feature$data$thick_end, c(18L, 9L))
    expect_identical(feature$data$item_rgb, c("255,0,128", "0"))
    expect_identical(feature$data$name, c("same", "same"))
    restored <- as_feature(as_granges(feature, level = "feature"))
    fields <- c("block_count", "block_sizes", "block_starts", "thick_start",
        "thick_end", "item_rgb")
    expect_equal(as.data.frame(restored$data[, ..fields]),
        as.data.frame(feature$data[, ..fields]))
})

test_that("BED12 read write round trips preserve structural metadata and repeated names", {
    file <- annotation_bed12_fixture()
    output <- tempfile(fileext = ".bed12")
    on.exit(unlink(c(file, output)), add = TRUE)
    feature <- read_bed(file, verbose = FALSE)
    before <- serialize(feature, NULL)
    write_feature(feature, output, format = "bed12")
    observed <- data.table::fread(output, header = FALSE,
        colClasses = list(character = c(9L, 11L, 12L)))
    expect_equal(nrow(observed), 2L)
    expect_identical(observed$V1, c("chr1", "chr2"))
    expect_identical(observed$V4, c("same", "same"))
    expect_equal(observed$V2, c(0L, 9L))
    expect_equal(observed$V7, c(1L, 9L))
    expect_equal(observed$V8, c(18L, 9L))
    expect_identical(observed$V9, c("255,0,128", "0"))
    expect_identical(observed$V11, c("5,11", "11"))
    expect_identical(observed$V12, c("0,9", "0"))
    restored <- read_bed(output, verbose = FALSE)
    fields <- c("chrom", "start", "end", "name", "score", "strand",
        "block_count", "block_sizes", "block_starts", "thick_start", "thick_end", "item_rgb")
    expect_equal(as.data.frame(restored$data[, ..fields]),
        as.data.frame(feature$data[, ..fields]))
    expect_identical(serialize(feature, NULL), before)
    write_feature(feature, output, format = "bed12", coordinate = "granges",
        overwrite = TRUE)
    restored <- read_bed(output, coordinate = "granges", verbose = FALSE)
    expect_equal(as.data.frame(restored$data[, ..fields]),
        as.data.frame(feature$data[, ..fields]))
})

test_that("BED12 metadata stays row aligned after selection and handles empty output", {
    file <- annotation_bed12_fixture()
    output <- tempfile(fileext = ".bed12")
    on.exit(unlink(c(file, output)), add = TRUE)
    feature <- read_bed(file, verbose = FALSE)
    selected <- retrieve_feature(feature, chrom = "chr2", start = 10L, end = 20L)
    write_feature(selected, output, format = "bed12")
    restored <- read_bed(output, verbose = FALSE)
    expect_identical(restored$data$block_sizes, list(11L))
    expect_identical(restored$data$thick_end, 9L)
    selected$data <- data.table::copy(selected$data[0L])
    write_feature(selected, output, format = "bed12", overwrite = TRUE)
    expect_equal(file.info(output)$size, 0)
})

test_that("BED12 import rejects invalid structural fields rather than losing records", {
    valid <- c("chr1", "0", "20", "tx", "0", "+", "1", "18",
        "0", "2", "5,11", "0,9")
    cases <- list(c(10L, "3"), c(11L, "5,0"), c(11L, "5,10.5"),
        c(12L, "1,9"), c(12L, "0,4"), c(12L, "0,10"),
        c(7L, "21"), c(8L, "21"), c(9L, "256,0,0"),
        c(11L, "2147483648,11"), c(11L, "5,,11"),
        c(2L, "0.5"), c(2L, "2147483647"), c(3L, "2147483648"))
    file <- tempfile(fileext = ".bed")
    on.exit(unlink(file), add = TRUE)
    for (case in cases) {
        invalid <- valid
        invalid[as.integer(case[1L])] <- case[2L]
        writeLines(paste(invalid, collapse = "\t"), file)
        expect_error(read_bed(file, verbose = FALSE), "BED12")
    }
})

test_that("BED12 export rejects inconsistent metadata edits before writing", {
    file <- annotation_bed12_fixture()
    on.exit(unlink(file), add = TRUE)
    feature <- read_bed(file, verbose = FALSE)
    feature$data$block_sizes[[1L]] <- c(5, 10.5)
    output <- tempfile(fileext = ".bed12")
    on.exit(unlink(output), add = TRUE)
    expect_error(write_feature(feature, output, format = "bed12"), "BED12")
    expect_false(file.exists(output))
    feature <- read_bed(file, verbose = FALSE)
    feature$data$block_starts <- NULL
    expect_error(write_feature(feature, output, format = "bed12"), "BED12.*missing|BED12.*incomplete")
    expect_false(file.exists(output))
})

test_that("GFF rebuild guards include gene parents and secondary transcript links", {
    file <- tempfile(fileext = ".gff3")
    on.exit(unlink(file), add = TRUE)
    writeLines(c("##gff-version 3",
        "chr1\ts\tbiological_region\t1\t30\t.\t+\t.\tID=r1",
        "chr1\ts\tbiological_region\t1\t30\t.\t+\t.\tID=r2",
        "chr1\ts\tgene\t1\t20\t.\t+\t.\tID=g1;Parent=r1,r2",
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\ts\texon\t1\t5\t.\t+\t.\tID=e1;Parent=t1",
        "chr1\ts\tgene\t40\t50\t.\t+\t.\tID=g2",
        "chr1\ts\tmRNA\t40\t50\t.\t+\t.\tID=t3;Parent=g2",
        "chr1\ts\texon\t40\t50\t.\t+\t.\tID=e3;Parent=t3"), file)
    feature <- read_gff(file, verbose = FALSE)
    expect_error(retrieve_feature(feature, gene_id = "g1"), "parent.*GRanges")
    unrelated <- retrieve_feature(feature, gene_id = "g2")
    expect_identical(unrelated$transcripts$transcript_id, "t3")
    shared_file <- annotation_parent_fixture()
    on.exit(unlink(shared_file), add = TRUE)
    shared <- read_gff(shared_file, verbose = FALSE)
    expect_error(retrieve_feature(shared, transcript_id = "t2"), "parent.*GRanges")
})

test_that("GFF merge renaming updates every parent token and supports legacy rows", {
    file <- annotation_parent_fixture()
    output <- tempfile(fileext = ".gff3")
    on.exit(unlink(c(file, output)), add = TRUE)
    feature <- read_gff(file, verbose = FALSE)
    before <- data.table::copy(feature$data)
    expect_warning(merged <- merge_feature(feature, feature,
        source_names = c("A", "B"), conflict = "rename"), "Duplicated annotation")
    row <- which(merged$data$feature_id == "B_e1")
    expect_length(row, 1L)
    expect_identical(merged$data$parent_id[row], "B_t1")
    expect_identical(merged$data$parent_ids[[row]], c("B_t1", "B_t2"))
    expect_no_error(write_feature(merged, output, format = "gff"))
    restored <- read_gff(output, verbose = FALSE)
    row <- which(restored$data$feature_id == "B_e1")
    expect_identical(restored$data$parent_ids[[row]], c("B_t1", "B_t2"))
    expect_equal(feature$data, before)
    legacy <- feature
    legacy$data <- data.table::copy(feature$data)
    legacy$data$parent_ids <- NULL
    expect_warning(mixed <- merge_feature(legacy, feature,
        source_names = c("A", "B"), conflict = "rename"), "Duplicated annotation")
    expect_no_error(write_feature(mixed, output, format = "gff", overwrite = TRUE))
    row <- which(mixed$data$feature_id == "B_t1")
    expect_identical(mixed$data$parent_ids[[row]], "B_g1")
})

test_that("primary analysis plots and haplotype location preserve legacy behavior", {
    file <- annotation_parent_fixture()
    on.exit(unlink(file), add = TRUE)
    feature <- read_gff(file, verbose = FALSE)
    legacy <- feature
    legacy$data <- data.table::copy(feature$data)
    legacy$data$parent_ids <- NULL
    before <- lapply(feature, function(x) if (is.data.frame(x)) as.data.frame(x) else x)
    for (level in c("feature", "gene", "transcript", "exon")) {
        expect_equal(summary_feature(feature, level = level),
            summary_feature(legacy, level = level))
    }
    expect_equal(get_gene_length_distribution_table(feature, chrom = "chr1",
        start = 1L, end = 20L), get_gene_length_distribution_table(legacy,
            chrom = "chr1", start = 1L, end = 20L))
    plot <- plot_feature_track(feature, "chr1", 1L, 20L)
    expected <- plot_feature_track(legacy, "chr1", 1L, 20L)
    expect_equal(ggplot2::ggplot_build(plot)$data, ggplot2::ggplot_build(expected)$data)
    plot <- plot_region(feature, "chr1", 1L, 20L)
    expected <- plot_region(legacy, "chr1", 1L, 20L)
    expect_equal(ggplot2::ggplot_build(plot)$data, ggplot2::ggplot_build(expected)$data)
    expect_equal(GeneTrackR:::resolve_haplotype_gene_region(feature, gene_id = "g1"),
        GeneTrackR:::resolve_haplotype_gene_region(legacy, gene_id = "g1"))
    expect_equal(GeneTrackR:::resolve_haplotype_gene_region(feature, transcript_id = "t1"),
        GeneTrackR:::resolve_haplotype_gene_region(legacy, transcript_id = "t1"))
    hap <- list(region = list(chrom = "chr1", start = 1L, end = 20L))
    vars <- data.table::data.table(pos = c(2L, 4L))
    mapper <- function(x) as.numeric(x)
    track <- GeneTrackR:::prepare_hap_gene_track(feature, hap, vars, mapper)
    expected <- GeneTrackR:::prepare_hap_gene_track(legacy, hap, vars, mapper)
    expect_false(is.null(track))
    expect_equal(track, expected)
    expect_equal(lapply(feature, function(x) if (is.data.frame(x)) as.data.frame(x) else x),
        before)
})
