signal_ops_fixture <- function(sample_ids = c("A", "B")) {
    signal <- data.table::data.table(
        sample_id = rep(c("A", "B"), each = 2L), chrom = "chr1",
        start = rep(c(1L, 5L), 2L), end = rep(c(2L, 7L), 2L),
        value = c(2, 6, 4, 8), strand = rep(c("+", "-"), each = 2L)
    )
    BwgTrack(
        samples = data.frame(sample_id = sample_ids,
            strand = c("+", "-")[match(sample_ids, c("A", "B"))]),
        data = signal[signal$sample_id %in% sample_ids],
        seqinfo = data.frame(sample_id = sample_ids, chrom = "chr1",
            length = 100L), meta = list(mode = "memory")
    )
}

test_that("signal normalization aligns named factors without mutating inputs", {
    track <- signal_ops_fixture()
    before <- serialize(track, NULL)
    for (method in c("RPM", "CPM")) {
        normalized <- norm_bwg(track, method = method,
            library_size = c(B = 200, A = 100))
        expect_equal(normalized$data$value, c(20000, 60000, 20000, 40000))
        expect_equal(normalized$samples$scale_factor, c(10000, 5000))
        expect_identical(normalized$samples$norm_method, rep(method, 2L))
    }
    automatic <- norm_bwg(track, method = "RPM")
    expect_equal(automatic$samples$scale_factor, c(1e6 / 8, 1e6 / 12))
    scaled <- norm_bwg(track, method = "scale", scale_factor = c(B = 3, A = 2))
    expect_equal(scaled$data$value, c(4, 12, 12, 24))
    custom <- norm_bwg(track, method = "custom", custom_factor = c(B = 0.5))
    expect_equal(custom$data$value, c(2, 6, 2, 4))
    expect_identical(norm_bwg(track, method = "none"), track)
    expect_identical(serialize(track, NULL), before)
    expect_error(norm_bwg(track, method = "scale"), "scale_factor")
    expect_error(norm_bwg(track, method = "custom"), "custom_factor")
    lazy <- track
    lazy$data <- NULL
    lazy$meta$mode <- "lazy"
    expect_error(norm_bwg(lazy, method = "RPM"), "library_size")
    expect_equal(norm_bwg(lazy, method = "CPM",
        library_size = c(A = 100, B = 200))$samples$scale_factor, c(10000, 5000))
})

test_that("signal summaries use inclusive widths and slices retain metadata", {
    track <- signal_ops_fixture()
    snapshot <- function(x) lapply(x,
        function(value) if (is.data.frame(value)) as.data.frame(value) else value)
    before <- snapshot(track)
    summary <- summary_bwg(track, samples = "A")
    expect_identical(summary$sample_id, "A")
    expect_equal(summary$sum_signal, 22)
    expect_equal(summary$mean_signal, 22 / 5)
    expect_equal(summary$max_signal, 6)
    expect_equal(summary$covered_bases, 5)
    regional <- summary_bwg(track, "chr1", 2L, 5L, samples = "A")
    expect_equal(regional$sum_signal, 8)
    expect_equal(regional$covered_bases, 2)
    empty <- summary_bwg(track, "chr2", 1L, 10L)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty), names(summary))
    sliced <- slice_bwg(track, "chr1", 2L, 5L, samples = "B")
    expect_identical(sliced$samples$sample_id, "B")
    expect_identical(sliced$seqinfo$sample_id, "B")
    expect_equal(sliced$seqinfo$length, 100L)
    expect_equal(sliced$data$start, c(2L, 5L))
    expect_equal(sliced$data$end, c(2L, 5L))
    ranges <- slice_bwg(track, "chr1", 2L, 5L, samples = "B", as = "GRanges")
    expect_identical(as.character(GenomicRanges::strand(ranges)), c("-", "-"))
    expect_equal(S4Vectors::mcols(ranges)$score, c(4, 8))
    plain <- slice_bwg(track, "chr1", 2L, 5L, samples = "B", as = "data.frame")
    expect_equal(as.data.frame(plain), as.data.frame(sliced$data))
    expect_identical(snapshot(track), before)
    lazy <- track
    lazy$data <- NULL
    expect_equal(summary_bwg(lazy), lazy$samples)
})

test_that("signal merges resolve aligned duplicates and normalization conflicts", {
    first <- signal_ops_fixture("A")
    second <- signal_ops_fixture("A")
    second$data$value <- c(4, 8)
    before <- lapply(list(first, second), serialize, connection = NULL)
    expect_error(merge_bwg(first, second), "Duplicated sample")
    for (policy in c("sum", "mean")) {
        merged <- merge_bwg(first, second, sample_conflict = policy)
        expect_identical(merged$samples$sample_id, "A")
        expect_equal(merged$data$value,
            if (policy == "sum") c(6, 14) else c(3, 7))
        expect_equal(nrow(merged$seqinfo), 1L)
        expect_equal(merged$seqinfo$length, 100L)
    }
    distinct <- signal_ops_fixture("B")
    distinct$samples$norm_method <- "custom"
    expect_error(merge_bwg(first, distinct), "different normalization")
    merged <- merge_bwg(first, distinct, require_same_norm = FALSE)
    expect_setequal(merged$samples$sample_id, c("A", "B"))
    expect_equal(nrow(merged$data), 4L)
    expect_identical(lapply(list(first, second), serialize, connection = NULL),
        before)
    expect_error(merge_bwg(), "At least one")
    expect_error(merge_bwg(first, list()), "BwgTrack")
})

test_that("signal bins restart by sample and chromosome and preserve empty schemas", {
    input <- signal_ops_fixture()$data
    before <- data.table::copy(input)
    bins <- bin_bwg(input, bin_size = 10L)
    expect_equal(nrow(bins), 2L)
    expect_identical(bins$sample_id, c("A", "B"))
    expect_equal(bins$value, c(4, 6))
    expect_equal(bins$start, c(1L, 1L))
    expect_equal(bins$end, c(7L, 7L))
    expect_equal(input, before)
    empty <- bin_bwg(input[0L])
    invalid <- data.table::copy(input)
    invalid$start <- NA_integer_
    expect_equal(bin_bwg(invalid), empty)
    expect_equal(nrow(empty), 0L)
    expect_identical(names(empty),
        c("sample_id", "chrom", "start", "end", "value", "bin"))
    expect_error(bin_bwg(data.frame(start = 1L)), "signal table")
})

test_that("regional signal tracks retain empty samples and reuse normalized values", {
    plus_file <- tempfile(fileext = ".bedgraph")
    empty_file <- tempfile(fileext = ".bedgraph")
    on.exit(unlink(c(plus_file, empty_file)), add = TRUE)
    writeLines(c("chr1\t0\t4\t2", "chr1\t5\t8\t-1"), plus_file)
    writeLines("track type=bedGraph", empty_file)

    for (mode in c("memory", "lazy")) {
        track <- read_bwg(c(plus_file, empty_file), sample_names = c("plus", "empty"),
            strand = c("+", "*"), mode = mode, use_tabix = FALSE, verbose = FALSE)
        track <- BwgTrack(track$samples, track$data, track$meta,
            seqinfo = data.frame(chrom = c("chr1", "chr2"), length = c(100L, 200L)))
        track <- norm_bwg(track, method = "scale", scale_factor = c(plus = 2))
        region <- retrieve_bwg(track, "chr1", 3L, 6L, as = "BwgTrack",
            keep_empty_samples = TRUE)
        expect_identical(region$meta$mode, "memory")
        expect_identical(as.list(region$data), list(
            sample_id = c("plus", "plus", "empty"), chrom = rep("chr1", 3L),
            start = c(3L, 6L, 3L), end = c(4L, 6L, 6L),
            value = c(4, -2, 0), strand = c("+", "+", "*")
        ))
        expect_identical(as.list(region$seqinfo), list(
            sample_id = c("plus", "empty"), chrom = c("chr1", "chr1"),
            length = c(100L, 100L)
        ))
        narrowed <- retrieve_bwg(region, "chr1", 4L, 4L, samples = "plus", as = "GRanges")
        expect_identical(IRanges::start(narrowed), 4L)
        expect_identical(IRanges::end(narrowed), 4L)
        expect_identical(S4Vectors::mcols(narrowed)$sample_id, "plus")
        expect_identical(S4Vectors::mcols(narrowed)$score, 4)
        gap <- retrieve_bwg(track, "chr1", 5L, 5L, strand = "+",
            strand_policy = "strict", keep_empty_samples = TRUE)
        expect_identical(as.list(gap), list(sample_id = "plus", chrom = "chr1",
            start = 5L, end = 5L, value = 0, strand = "+"))
    }
})

test_that("signal merge aggregation respects interval endpoints and strand", {
    first <- signal_ops_fixture("A")
    second <- signal_ops_fixture()
    second$data <- data.table::rbindlist(list(second$data, data.table::data.table(
        sample_id = "A", chrom = "chr1", start = 2L, end = 3L, value = -3, strand = "+"
    )))
    second$data[1L, value := 4]
    second$data[2L, `:=`(value = 10, strand = "-")]

    merged <- merge_bwg(first, second, sample_conflict = "sum")
    records <- merged$data[sample_id == "A", .(start, end, strand, value)]
    data.table::setorder(records, start, end, strand)
    expect_equal(records, data.table::data.table(
        start = c(1L, 2L, 5L, 5L), end = c(2L, 3L, 7L, 7L),
        strand = c("+", "+", "+", "-"), value = c(6, -3, 6, 10)
    ))
    expect_identical(merged$data$value[merged$data$sample_id == "B"], c(4, 8))
})

for (track_mode in c("memory", "lazy")) {
    test_that(paste("signal merge renaming preserves track identity in", track_mode, "mode"), {
        files <- replicate(5L, tempfile(fileext = ".bedgraph"))
        on.exit(unlink(files), add = TRUE)
        for (i in seq_along(files)) {
            writeLines(sprintf("chr1\t%d\t%d\t%d", i - 1L, i, 10L * i), files[i])
        }
        sample_ids <- c("A", "set1_A", "set8_A", "A", "set8_A")
        tracks <- lapply(list(1:3, 4:5), function(rows) {
            track <- read_bwg(files[rows], sample_names = sample_ids[rows], mode = track_mode,
                use_tabix = FALSE, verbose = FALSE)
            BwgTrack(track$samples, track$data, track$meta, seqinfo = data.frame(
                sample_id = sample_ids[rows], chrom = "chr1", length = 100L + rows
            ))
        })

        merged <- merge_bwg(tracks[[1L]], tracks[[2L]], sample_conflict = "rename")
        merged_ids <- merged$samples$sample_id
        expect_equal(anyDuplicated(merged_ids), 0L)
        expect_match(merged_ids[1L], "^set1_A")
        expect_identical(merged_ids[-1L], c("set1_A", "set1_set8_A", "set2_A", "set2_set8_A"))
        expect_identical(merged$samples$file, normalizePath(files, winslash = "/"))
        expect_identical(merged$seqinfo$sample_id, merged_ids)
        expect_identical(merged$seqinfo$length, 101:105)
        expect_identical(merged$meta$mode, track_mode)
        if (track_mode == "lazy") {
            expect_null(merged$data)
        }
        observed <- retrieve_bwg(merged, "chr1", 1L, 5L)
        expect_equal(nrow(observed), 5L)
        observed_rows <- match(merged_ids, observed$sample_id)
        expect_identical(observed$start[observed_rows], 1:5)
        expect_identical(observed$end[observed_rows], 1:5)
        expect_identical(observed$value[observed_rows], c(10, 20, 30, 40, 50))
    })
}

test_that("duplicate signal renaming rejects mixed storage before a sample is lost", {
    file <- tempfile(fileext = ".bedgraph")
    on.exit(unlink(file), add = TRUE)
    writeLines("chr1\t0\t4\t7", file)
    memory <- read_bwg(file, sample_names = "A", mode = "memory",
        use_tabix = FALSE, verbose = FALSE)
    lazy <- read_bwg(file, sample_names = "A", mode = "lazy",
        use_tabix = FALSE, verbose = FALSE)

    expect_error(
        merge_bwg(memory, lazy, sample_conflict = "rename"),
        "same storage mode"
    )
    expect_error(
        merge_bwg(lazy, memory, sample_conflict = "rename"),
        "same storage mode"
    )
})

test_that("signal bins use start offsets within each sample and chromosome", {
    input <- data.frame(
        sample_id = c(rep("A", 4L), rep("B", 3L)),
        chrom = c(rep("chr1", 3L), "chr2", rep("chr1", 3L)),
        start = c(3L, 7L, 8L, 21L, 101L, 105L, NA_integer_),
        end = c(4L, 7L, 12L, 23L, 102L, 105L, 110L),
        value = c(4, 8, -2, 10, 20, NA_real_, 900)
    )
    expect_equal(as.list(bin_bwg(input, bin_size = 5L)), list(
        sample_id = c("A", "A", "A", "B"), chrom = c("chr1", "chr1", "chr2", "chr1"),
        bin = c(1, 2, 1, 1), start = c(3L, 8L, 21L, 101L), end = c(7L, 12L, 23L, 105L),
        value = c(6, -2, 10, 20)
    ))
})

test_that("tabix diagnostics and benchmarks describe unindexed bedGraph queries", {
    file <- tempfile(fileext = ".bedgraph")
    on.exit(unlink(file), add = TRUE)
    writeLines(c("chr1\t0\t4\t9", "chr1\t5\t10\t-2"), file)
    legacy <- BwgTrack(data.frame(sample_id = "signal", file = file, format = "bedgraph"))
    track <- read_bwg(file, sample_names = "signal", mode = "lazy", use_tabix = FALSE,
        verbose = FALSE)

    diagnostic <- diagnose_tabix(legacy)
    flags <- list(has_tabix = FALSE, use_tabix = FALSE, tabix_backend = NA_character_,
        tabix_empty_fallback = FALSE)
    expect_identical(diagnostic$sample_id, "signal")
    expect_identical(as.list(diagnostic)[names(flags)], flags)
    benchmark <- benchmark_tabix_query(track, "chr1", 4L, 6L, samples = "signal")
    expect_identical(benchmark$sample_id, "signal")
    expect_false(benchmark$use_tabix)
    expect_identical(benchmark$n_records, 2L)
    expect_equal(nrow(benchmark_tabix_query(track, "chr1", 4L, 6L, samples = "absent")), 0L)
})

test_that("indexed bedGraph diagnostics preserve a one-base gap during fallback", {
    skip_if_not_installed("Rsamtools")

    file <- tempfile(fileext = ".bedgraph")
    compressed <- paste0(file, ".gz")
    on.exit(unlink(c(file, compressed, paste0(compressed, ".tbi"))), add = TRUE)
    writeLines(c("chr1\t0\t4\t9", "chr1\t5\t10\t-2"), file)
    Rsamtools::bgzip(file, dest = compressed)
    Rsamtools::indexTabix(compressed, format = "bed")

    track <- read_bwg(compressed, sample_names = "signal", mode = "lazy", verbose = FALSE)
    track$samples$tabix_backend <- "Rsamtools"
    diagnostic <- diagnose_tabix(track)
    expect_true(diagnostic$has_tabix)
    expect_true(diagnostic$use_tabix)
    expect_identical(diagnostic$tabix_backend, "Rsamtools")
    expect_identical(benchmark_tabix_query(track, "chr1", 4L, 6L)$n_records, 2L)
    gap <- retrieve_bwg(track, "chr1", 5L, 5L, tabix_empty_fallback = TRUE)
    expect_identical(as.list(gap), list(sample_id = character(), chrom = character(),
        start = integer(), end = integer(), value = numeric(), strand = character()))
})
