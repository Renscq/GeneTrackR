test_that("GenePred and BED12 exports preserve inclusive exon boundaries", {
    model <- coverage_gene_model()
    before <- serialize(model, NULL)
    files <- character()
    on.exit(unlink(files), add = TRUE)
    for (coordinate in c("ucsc", "granges")) {
        file <- tempfile(fileext = ".genePredExt")
        files <- c(files, file)
        write_feature(model, file, coordinate = coordinate, format = "auto")
        table <- data.table::fread(file, header = FALSE)
        expect_identical(table$V1, c("plus_tx", "minus_tx", "noncoding_tx"))
        expect_equal(table$V4,
            if (coordinate == "ucsc") c(0L, 20L, 40L) else c(1L, 21L, 41L))
        expect_equal(table$V5, c(12L, 32L, 44L))
        expect_identical(table$V9,
            if (coordinate == "ucsc") c("0,8,", "20,28,", "40,") else
                c("1,9,", "21,29,", "41,"))
    }
    file <- tempfile(fileext = ".bed12")
    files <- c(files, file)
    write_feature(model, file, format = "auto")
    table <- data.table::fread(file, header = FALSE,
        colClasses = list(character = c(11L, 12L)))
    expect_equal(table$V2, c(0L, 20L, 40L))
    expect_equal(table$V7, c(2L, 22L, 40L))
    expect_equal(table$V8, c(10L, 30L, 40L))
    expect_equal(table$V10, c(2L, 2L, 1L))
    expect_identical(table$V11, c("5,4", "5,4", "4"))
    expect_identical(table$V12, c("0,8", "0,8", "0"))
    expect_identical(serialize(model, NULL), before)
})

test_that("BED6 exports normalize scores and strands with explicit chromosome order", {
    input <- data.frame(
        feature_id = c("f10", "f2", "f1"), name = c("ten", "two", "one"),
        chrom = c("chr10", "chr2", "chr1"), start = c(11L, 3L, 1L),
        end = c(12L, 5L, 1L), strand = c("*", "-", "+"),
        score = c(NA, 1001, -2), type = "region"
    )
    feature <- Feature(input)
    before <- serialize(feature, NULL)
    file <- tempfile(fileext = ".bed6")
    fai <- tempfile(fileext = ".fai")
    on.exit(unlink(c(file, fai)), add = TRUE)
    writeLines(c("chr10\t100", "chr1\t100", "chr2\t100"), fai)
    write_feature(feature, file, format = "auto", chrom_order = fai)
    table <- data.table::fread(file, header = FALSE)
    expect_identical(table$V1, c("chr10", "chr1", "chr2"))
    expect_equal(table$V2, c(10L, 0L, 2L))
    expect_equal(table$V5, c(0L, 0L, 1000L))
    expect_identical(table$V6, c(".", "+", "-"))
    write_feature(feature, file, format = "bed6", coordinate = "granges",
        chrom_order = data.frame(chrom = c("chr2", "chr10", "chr1")),
        overwrite = TRUE)
    table <- data.table::fread(file, header = FALSE)
    expect_identical(table$V1, c("chr2", "chr10", "chr1"))
    expect_equal(table$V2, c(3L, 11L, 1L))
    expect_identical(serialize(feature, NULL), before)
    expect_error(write_feature(feature, file, format = "bed6"), "exist")
})

test_that("GFF and GTF exports preserve or rebuild attributes explicitly", {
    file <- tempfile(fileext = ".gff3")
    output <- tempfile(fileext = ".gff3")
    gtf <- tempfile(fileext = ".gtf")
    on.exit(unlink(c(file, output, gtf)), add = TRUE)
    writeLines(c("##gff-version 3",
        "chr1\tfixture\tgene\t1\t20\t.\t+\t.\tID=g1;Name=Gene1;Note=keep",
        "chr1\tfixture\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1",
        "chr1\tfixture\texon\t1\t5\t3\t+\t.\tID=e1;Parent=t1",
        "chr1\tfixture\tCDS\t2\t5\t.\t+\t0\tID=c1;Parent=t1"), file)
    feature <- read_gff(file, verbose = FALSE, progress = FALSE)
    before <- serialize(feature, NULL)
    write_feature(feature, output, format = "auto", sort_output = FALSE)
    observed <- data.table::fread(output, header = FALSE)
    expect_identical(observed$V9, feature$data$attribute)
    expect_true(any(grepl("Note=keep", observed$V9, fixed = TRUE)))
    write_feature(feature, output, format = "gff", keep_attributes = FALSE,
        overwrite = TRUE)
    rebuilt <- data.table::fread(output, header = FALSE)
    expect_false(any(grepl("Note=keep", rebuilt$V9, fixed = TRUE)))
    expect_true(any(grepl("ID=g1", rebuilt$V9, fixed = TRUE)))
    write_feature(feature, gtf, format = "auto", keep_attributes = FALSE)
    observed <- data.table::fread(gtf, header = FALSE, quote = "")
    expect_true(all(grepl('gene_id "g1"', observed$V9, fixed = TRUE)))
    expect_true(any(observed$V3 == "transcript"))
    expect_false(any(observed$V3 == "mRNA"))
    expect_identical(serialize(feature, NULL), before)
    expect_error(write_feature(feature, tempfile(fileext = ".unknown")),
        "Cannot infer")
    expect_error(write_feature(feature, tempfile(), format = "bed6",
        chrom_order = 1L), "chrom_order")
})
