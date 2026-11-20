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
