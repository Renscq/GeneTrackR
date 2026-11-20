#!/usr/bin/env Rscript

# Author: Rensc
# Date: 2026-11-20
# Version: dev002
# Function: Compare existing readers with standard readers on shared fixtures.
# Input: A source checkout and optional developer-only reader dependencies.
# Output: Record accounting, field differences, timings, schemas and provenance.

compare_io_tables <- function(native, standard, keys, fields, tolerance = 0) {
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance < 0) {
        stop("Invalid non-negative numeric tolerance.", call. = FALSE)
    }
    native <- as.data.frame(native, stringsAsFactors = FALSE)
    standard <- as.data.frame(standard, stringsAsFactors = FALSE)
    required <- unique(c(keys, fields))
    if (!length(keys) || !all(required %in% names(native)) ||
        !all(required %in% names(standard))) {
        stop("Required comparison columns are missing.", call. = FALSE)
    }
    make_keys <- function(table) {
        columns <- lapply(table[keys], as.character)
        invalid <- vapply(columns, function(x)
            anyNA(x) || any(!nzchar(x)) || any(grepl("\r", x, fixed = TRUE)),
            logical(1))
        if (any(invalid)) stop("Keys contain missing or invalid values.", call. = FALSE)
        values <- do.call(paste, c(columns, sep = "\r"))
        if (anyDuplicated(values)) stop("Duplicate comparison keys.", call. = FALSE)
        values
    }
    native_key <- make_keys(native)
    standard_key <- make_keys(standard)
    index <- match(native_key, standard_key)
    matched <- which(!is.na(index))
    native_only <- native[is.na(index), , drop = FALSE]
    standard_only <- standard[!standard_key %in% native_key, , drop = FALSE]
    differences <- data.frame(key = character(), field = character(),
        native = character(), standard = character())
    for (field in fields) {
        left <- native[[field]][matched]
        right <- standard[[field]][index[matched]]
        same <- is.na(left) & is.na(right)
        present <- !is.na(left) & !is.na(right)
        if (is.numeric(left) && is.numeric(right)) {
            same[present] <- left[present] == right[present] |
                abs(left[present] - right[present]) <= tolerance
        } else {
            same[present] <- as.character(left[present]) == as.character(right[present])
        }
        different <- which(!same)
        if (length(different)) {
            differences <- rbind(differences, data.frame(
                key = gsub("\r", " | ", native_key[matched[different]], fixed = TRUE),
                field = field, native = as.character(left[different]),
                standard = as.character(right[different])))
        }
    }
    equivalent <- !nrow(differences) && !nrow(native_only) && !nrow(standard_only)
    list(summary = data.frame(
        status = if (equivalent) "equivalent" else "different",
        native_rows = nrow(native), standard_rows = nrow(standard),
        matched_rows = length(matched), native_only_rows = nrow(native_only),
        standard_only_rows = nrow(standard_only),
        different_values = nrow(differences)),
        differences = differences, native_only = native_only,
        standard_only = standard_only)
}

measure_io_reader <- function(reader, iterations = 3L) {
    if (length(iterations) != 1L || !is.numeric(iterations) ||
        !is.finite(iterations) || iterations < 1 || iterations != floor(iterations)) {
        stop("`iterations` must be a positive integer.", call. = FALSE)
    }
    warnings <- character()
    timings <- data.frame(iteration = integer(), elapsed_seconds = numeric(),
        object_bytes = numeric())
    value <- NULL
    error <- NULL
    for (iteration in seq_len(iterations)) {
        gc()
        started <- proc.time()[["elapsed"]]
        result <- tryCatch(withCallingHandlers(reader(), warning = function(condition) {
            warnings <<- c(warnings, conditionMessage(condition))
            invokeRestart("muffleWarning")
        }), error = function(condition) {
            error <<- conditionMessage(condition)
            NULL
        })
        if (!is.null(error)) break
        timings <- rbind(timings, data.frame(iteration = iteration,
            elapsed_seconds = proc.time()[["elapsed"]] - started,
            object_bytes = as.numeric(object.size(result))))
        value <- result
    }
    list(status = if (is.null(error)) "ok" else "error", value = value,
        timing = timings, warnings = paste(unique(warnings), collapse = " | "),
        error = if (is.null(error)) "" else error)
}

io_missing <- function(x) {
    x <- as.character(x)
    x[!is.na(x) & x %in% c("", ".")] <- NA_character_
    x
}

io_column <- function(table, name) {
    if (!name %in% names(table)) return(rep(NA_character_, nrow(table)))
    value <- table[[name]]
    if (is.list(value) || methods::is(value, "List")) {
        value <- vapply(as.list(value), function(x)
            if (!length(x) || all(is.na(x))) NA_character_ else
                paste(x, collapse = ","), character(1))
    }
    io_missing(value)
}

io_ranges_frame <- function(object) {
    table <- data.frame(chrom = as.character(GenomicRanges::seqnames(object)),
        start = IRanges::start(object), end = IRanges::end(object),
        strand = as.character(GenomicRanges::strand(object)))
    metadata <- S4Vectors::mcols(object)
    for (name in names(metadata)) table[[name]] <- metadata[[name]]
    table
}

io_with_local_path <- function(file, reader) {
    file <- normalizePath(file, winslash = "/", mustWork = TRUE)
    before <- getwd()
    on.exit(setwd(before), add = TRUE)
    setwd(dirname(file))
    reader(basename(file))
}

io_annotation_table <- function(object, format, backend) {
    if (backend == "native") {
        table <- as.data.frame(object$data)
    } else if (backend == "Bioc.gff") {
        table <- as.data.frame(object)
        names(table)[names(table) == "seqid"] <- "chrom"
    } else {
        table <- io_ranges_frame(object)
    }
    out <- data.frame(chrom = as.character(table$chrom), start = as.integer(table$start),
        end = as.integer(table$end), strand = as.character(table$strand))
    out$strand[out$strand == "."] <- "*"
    out$score <- as.numeric(table$score)
    if (format %in% c("bed6", "bed12")) {
        out$name <- io_column(table, "name")
        if (format == "bed12") {
            out$block_sizes <- NA_character_
            out$block_starts <- NA_character_
            out$thick_start <- NA_integer_
            out$thick_end <- NA_integer_
            if (backend == "native") {
                if (all(c("block_sizes", "block_starts", "thick_start", "thick_end") %in% names(table))) {
                    out$block_sizes <- io_column(table, "block_sizes")
                    out$block_starts <- io_column(table, "block_starts")
                    out$thick_start <- as.integer(table$thick_start) - 1L
                    out$thick_end <- as.integer(table$thick_end)
                }
            } else {
                out$block_sizes <- vapply(as.list(table$blocks), function(x)
                    paste(IRanges::width(x), collapse = ","), character(1))
                out$block_starts <- vapply(as.list(table$blocks), function(x)
                    paste(IRanges::start(x) - 1L, collapse = ","), character(1))
                out$thick_start <- IRanges::start(table$thick) - 1L
                out$thick_end <- IRanges::end(table$thick)
            }
        }
    } else {
        out$type <- as.character(table$type)
        out$source <- as.character(table$source)
        out$phase <- io_column(table, "phase")
        if (format == "gff3") {
            out$id <- io_column(table, if (backend == "native") "feature_id" else "ID")
            out$name <- io_column(table, if (backend == "native") "name" else "Name")
            parent_field <- if (backend == "native") {
                if ("parent_ids" %in% names(table)) "parent_ids" else "parent_id"
            } else "Parent"
            out$parent <- io_column(table, parent_field)
        } else {
            out$gene_id <- io_column(table, "gene_id")
            out$transcript_id <- io_column(table, "transcript_id")
        }
    }
    out
}

io_signal_table <- function(object, backend) {
    table <- if (backend == "native") as.data.frame(object$data) else io_ranges_frame(object)
    data.frame(chrom = as.character(table$chrom),
        start = as.integer(table$start), end = as.integer(table$end),
        value = as.numeric(table[[if (backend == "native") "value" else "score"]]))
}

io_info_values <- function(info, field) {
    vapply(info, function(value) {
        entries <- strsplit(value, ";", fixed = TRUE)[[1L]]
        if (field == "FLAG") return(as.character("FLAG" %in% entries))
        entry <- entries[startsWith(entries, paste0(field, "="))]
        if (!length(entry)) return(NA_character_)
        io_missing(substring(entry[1L], nchar(field) + 2L))
    }, character(1), USE.NAMES = FALSE)
}

io_vcf_table <- function(object, backend, samples) {
    if (backend == "native") {
        table <- object$data
        out <- data.frame(chrom = table$chrom, pos = table$pos,
            ref = table$ref, alt = table$alt, id = table$variant_id,
            end = table$end, qual = table$qual, filter = io_missing(table$filter))
        for (field in c("DP", "AF", "END", "FLAG")) {
            out[[paste0("INFO:", field)]] <- io_info_values(table$info, field)
        }
        for (sample in samples) for (field in c("GT", "AD", "PL")) {
            out[[paste(field, sample, sep = ":")]] <- vapply(seq_len(nrow(table)), function(i) {
                format <- strsplit(table$FORMAT[i], ":", fixed = TRUE)[[1L]]
                values <- strsplit(table[[sample]][i], ":", fixed = TRUE)[[1L]]
                position <- match(field, format)
                if (is.na(position) || position > length(values)) NA_character_ else
                    io_missing(values[position])
            }, character(1))
        }
    } else {
        ranges <- SummarizedExperiment::rowRanges(object)
        out <- data.frame(chrom = as.character(GenomicRanges::seqnames(ranges)),
            pos = IRanges::start(ranges), ref = as.character(VariantAnnotation::ref(object)),
            alt = vapply(as.list(VariantAnnotation::alt(object)), paste,
                character(1), collapse = ","), id = names(ranges),
            end = IRanges::end(ranges), qual = VariantAnnotation::qual(object),
            filter = io_missing(VariantAnnotation::filt(object)))
        for (field in c("DP", "AF", "END", "FLAG")) {
            out[[paste0("INFO:", field)]] <- io_column(VariantAnnotation::info(object), field)
        }
        genotype <- VariantAnnotation::geno(object)
        for (sample in samples) for (field in c("GT", "AD", "PL")) {
            values <- genotype[[field]]
            out[[paste(field, sample, sep = ":")]] <- vapply(seq_len(nrow(object)), function(i) {
                value <- if (length(dim(values)) == 3L) values[i, sample, ] else values[i, sample]
                value <- unlist(value, use.names = FALSE)
                if (!length(value) || all(is.na(value))) return(NA_character_)
                text <- as.character(value)
                text[is.na(value)] <- "."
                paste(text, collapse = ",")
            }, character(1))
        }
    }
    out
}

standard_io_main <- function(args = commandArgs(trailingOnly = TRUE)) {
    help <- paste("Usage: Rscript tools/compare_standard_io.R [options]",
        "  --path DIR          Source checkout (default: current directory)",
        "  --output-dir DIR    Fresh report directory (required)",
        "  --iterations N      Warm-namespace reads per reader (default: 3)",
        "  --help              Show this help",
        "Requires optional Bioc.gff, rtracklayer and VariantAnnotation tools.",
        "Reports differences and reader errors; never changes a runtime backend.",
        sep = "\n")
    if (any(args %in% c("--help", "-h"))) {
        cat(help, "\n")
        return(invisible(NULL))
    }
    options <- list(path = ".", output_dir = NULL, iterations = "3")
    mapping <- c("--path" = "path", "--output-dir" = "output_dir",
        "--iterations" = "iterations")
    index <- 1L
    while (index <= length(args)) {
        if (!args[index] %in% names(mapping) || index == length(args)) {
            stop(help, call. = FALSE)
        }
        options[[mapping[[args[index]]]]] <- args[index + 1L]
        index <- index + 2L
    }
    iterations <- as.numeric(options$iterations)
    if (length(iterations) != 1L || !is.finite(iterations) ||
        iterations < 1 || iterations != floor(iterations)) {
        stop("--iterations must be a positive integer.", call. = FALSE)
    }
    if (is.null(options$output_dir) || file.exists(options$output_dir)) {
        stop("Supply a fresh --output-dir to retain earlier reports.", call. = FALSE)
    }
    packages <- c("GeneTrackR", "Bioc.gff", "rtracklayer", "VariantAnnotation")
    available <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
    if (!all(available)) stop(paste("Missing evaluation packages:",
        paste(packages[!available], collapse = ", ")), call. = FALSE)
    source <- normalizePath(options$path, winslash = "/", mustWork = TRUE)
    source_version <- unname(read.dcf(file.path(source, "DESCRIPTION"))[1L, "Version"])
    if (as.character(utils::packageVersion("GeneTrackR")) != source_version) {
        stop("Installed GeneTrackR version does not match the source checkout.", call. = FALSE)
    }
    fixture_dir <- file.path(source, "tools", "io-evaluation", "fixtures")
    fixtures <- c(gff3 = "annotation.gff3", gtf = "annotation.gtf",
        bed6 = "intervals.bed", bed12 = "transcripts.bed", vcf = "variants.vcf")
    files <- file.path(fixture_dir, fixtures)
    names(files) <- names(fixtures)
    files <- c(files, bigwig = file.path(source, "tests", "testthat", "fixtures",
        "native_reader_test.bigwig"))
    stopifnot(all(file.exists(files)))
    dir.create(options$output_dir, recursive = TRUE)
    output <- normalizePath(options$output_dir, winslash = "/", mustWork = TRUE)
    dir.create(file.path(output, "fixtures"))
    stopifnot(all(file.copy(files, file.path(output, "fixtures", basename(files)))))
    utils::write.csv(data.frame(format = names(files), filename = basename(files),
        md5 = unname(tools::md5sum(files))), file.path(output, "fixtures.csv"), row.names = FALSE)
    description <- utils::packageDescription("GeneTrackR")
    metadata <- c(Package = "GeneTrackR", Version = description$Version,
        R = R.version.string, Platform = R.version$platform,
        SourceRef = "working-tree with recorded source checksums",
        Iterations = as.character(iterations), Timing = "warm namespaces; raw reads; gc outside timing",
        Memory = "returned R object size in bytes, not peak process memory",
        PathHandling = "standard import uses fixture directory and basename; restores working directory",
        Tolerance = "absolute 1e-6 for numeric comparison; NA distinct from zero",
        CreatedUTC = format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
        stats::setNames(vapply(packages[-1L], function(x)
            as.character(utils::packageVersion(x)), character(1)), packages[-1L]))
    write.dcf(as.data.frame(as.list(metadata)), file.path(output, "metadata.dcf"))
    source_files <- c(file.path(source, "DESCRIPTION"),
        list.files(file.path(source, "R"), "[.]R$", full.names = TRUE),
        file.path(source, "tools", "compare_standard_io.R"))
    utils::write.csv(data.frame(filename = substring(source_files, nchar(source) + 2L),
        md5 = unname(tools::md5sum(source_files))),
        file.path(output, "source-checksums.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(output, "session-info.txt"))
    summaries <- list()
    timing <- list()
    statuses <- list()
    schemas <- list()
    for (format in names(files)) {
        file <- files[[format]]
        native_reader <- switch(format,
            gff3 = function() GeneTrackR::read_gff(file, verbose = FALSE, progress = FALSE),
            gtf = function() GeneTrackR::read_gtf(file, verbose = FALSE, progress = FALSE),
            bed6 = function() GeneTrackR::read_bed(file, verbose = FALSE, progress = FALSE),
            bed12 = function() GeneTrackR::read_bed(file, verbose = FALSE, progress = FALSE),
            vcf = function() GeneTrackR::read_vcf(file, mode = "memory", verbose = FALSE, progress = FALSE),
            bigwig = function() GeneTrackR::read_bwg(file, format = "bigwig",
                mode = "memory", verbose = FALSE))
        native <- measure_io_reader(native_reader, iterations)
        candidates <- if (format %in% c("gff3", "gtf")) c("rtracklayer", "Bioc.gff") else
            if (format == "vcf") "VariantAnnotation" else "rtracklayer"
        for (backend in c("native", candidates)) {
            measured <- if (backend == "native") native else measure_io_reader(
                switch(backend, rtracklayer = function() io_with_local_path(file,
                    function(path) rtracklayer::import(path,
                        format = if (format %in% c("bed6", "bed12")) "bed" else format)),
                    Bioc.gff = function() Bioc.gff::readGFF(file,
                        version = if (format == "gff3") 3L else 2L),
                    VariantAnnotation = function() VariantAnnotation::readVcf(file,
                        genome = NA_character_)), iterations)
            label <- paste(format, backend, sep = "-")
            statuses[[label]] <- data.frame(format = format, reader = backend,
                status = measured$status, error = measured$error, warnings = measured$warnings)
            timing[[label]] <- cbind(data.frame(format = rep(format, nrow(measured$timing)),
                reader = rep(backend, nrow(measured$timing))), measured$timing)
            saveRDS(measured$value, file.path(output, paste0(label, ".rds")))
            if (measured$status != "ok") next
            raw_columns <- if (backend == "native") names(measured$value$data) else
                if (format == "vcf") c(names(VariantAnnotation::fixed(measured$value)),
                    names(VariantAnnotation::info(measured$value)),
                    names(VariantAnnotation::geno(measured$value))) else
                        if (backend == "rtracklayer") names(io_ranges_frame(measured$value)) else
                            names(measured$value)
            schemas[[label]] <- data.frame(format = format, reader = backend,
                column = raw_columns, object_class = paste(class(measured$value), collapse = "/"))
            if (backend == "native" || native$status != "ok") next
            if (format == "bigwig") {
                left <- io_signal_table(native$value, "native")
                right <- io_signal_table(measured$value, backend)
                keys <- c("chrom", "start", "end")
            } else if (format == "vcf") {
                samples <- native$value$meta$sample_names
                order <- VariantAnnotation::samples(VariantAnnotation::header(measured$value))
                utils::write.csv(data.frame(position = seq_along(samples),
                    native = samples, standard = order), file.path(output, "vcf-sample-order.csv"),
                    row.names = FALSE)
                stopifnot(identical(samples, order))
                left <- io_vcf_table(native$value, "native", samples)
                right <- io_vcf_table(measured$value, backend, samples)
                keys <- c("chrom", "pos", "ref", "alt")
            } else {
                left <- io_annotation_table(native$value, format, "native")
                right <- io_annotation_table(measured$value, format, backend)
                keys <- c("chrom", "start", "end", "strand",
                    if (format %in% c("gff3", "gtf")) "type" else "name",
                    if (format == "gff3") "id")
            }
            utils::write.csv(left, file.path(output, paste0(label, "-native.csv")), row.names = FALSE)
            utils::write.csv(right, file.path(output, paste0(label, "-standard.csv")), row.names = FALSE)
            fields <- setdiff(names(left), keys)
            result <- compare_io_tables(left, right, keys, fields, tolerance = 1e-6)
            summaries[[label]] <- cbind(data.frame(format = format, reader = backend,
                compared_fields = paste(c(keys, fields), collapse = ";")), result$summary)
            for (component in c("differences", "native_only", "standard_only")) {
                utils::write.csv(result[[component]], file.path(output,
                    paste0(label, "-", component, ".csv")), row.names = FALSE)
            }
        }
    }
    summary <- do.call(rbind, summaries)
    utils::write.csv(summary, file.path(output, "comparison.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, statuses), file.path(output, "reader-status.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, timing), file.path(output, "timing.csv"), row.names = FALSE)
    utils::write.csv(do.call(rbind, schemas), file.path(output, "schemas.csv"), row.names = FALSE)
    print(summary)
    cat("Differences and reader errors are assessment results, not equivalence claims.\n")
    invisible(summary)
}

if (sys.nframe() == 0L) standard_io_main()
