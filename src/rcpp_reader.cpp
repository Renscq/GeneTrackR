#include <Rcpp.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

extern "C" {
#include "bigWig.h"
}

using namespace Rcpp;

namespace {

const size_t kDefaultBufferSize = 128 * 1024;

std::string scalar_string(SEXP x, const char *name) {
  if (!Rf_isString(x) || Rf_length(x) != 1 || STRING_ELT(x, 0) == NA_STRING) {
    stop("`%s` must be a single non-missing character string.", name);
  }
  return Rcpp::as<std::string>(x);
}

int scalar_int(SEXP x, const char *name) {
  if (Rf_length(x) != 1 || Rf_isNull(x)) {
    stop("`%s` must be a single integer value.", name);
  }
  int value = Rcpp::as<int>(x);
  if (value == NA_INTEGER) {
    stop("`%s` must not be NA.", name);
  }
  return value;
}

bool is_remote_path(const std::string &file) {
  return file.rfind("http://", 0) == 0 ||
         file.rfind("https://", 0) == 0 ||
         file.rfind("ftp://", 0) == 0;
}

std::string sanitize_local_path(const std::string &file) {
  if (is_remote_path(file)) {
#ifdef NOCURL
    stop("Remote bigWig URLs are not supported by the current local-only libBigWig backend.");
#else
    return file;
#endif
  }

  const std::string file_prefix = "file://";
  if (file.rfind(file_prefix, 0) == 0) {
    stop("The libBigWig backend expects a local file path, not a file:// URI. Use a path like E:/path/file.bigwig.");
  }

  std::string out = file;
  std::replace(out.begin(), out.end(), '\\', '/');
  return out;
}

class BigWigSession {
 public:
  explicit BigWigSession(const std::string &file) : fp_(nullptr), initialized_(false) {
    if (bwInit(kDefaultBufferSize) != 0) {
      stop("Failed to initialize libBigWig.");
    }
    initialized_ = true;

    std::string path = sanitize_local_path(file);
    fp_ = bwOpen(path.c_str(), NULL, "r");
    if (fp_ == nullptr) {
      cleanup();
      stop("Failed to open bigWig file: %s", path.c_str());
    }
    if (fp_->cl == nullptr || fp_->cl->nKeys <= 0) {
      cleanup();
      stop("Failed to read chromosome information from bigWig file: %s", path.c_str());
    }
  }

  ~BigWigSession() {
    cleanup();
  }

  bigWigFile_t *get() const {
    return fp_;
  }

  chromList_t *chroms() const {
    return fp_->cl;
  }

  int64_t chrom_index(const std::string &chrom) const {
    chromList_t *cl = chroms();
    for (int64_t i = 0; i < cl->nKeys; ++i) {
      if (chrom == std::string(cl->chrom[i])) {
        return i;
      }
    }
    return -1;
  }

 private:
  void cleanup() {
    if (fp_ != nullptr) {
      bwClose(fp_);
      fp_ = nullptr;
    }
    if (initialized_) {
      bwCleanup();
      initialized_ = false;
    }
  }

  bigWigFile_t *fp_;
  bool initialized_;
};

DataFrame seqinfo_impl(const std::string &file) {
  BigWigSession bw(file);
  chromList_t *cl = bw.chroms();

  std::vector<std::string> chrom;
  std::vector<int> length;
  chrom.reserve(static_cast<size_t>(cl->nKeys));
  length.reserve(static_cast<size_t>(cl->nKeys));

  for (int64_t i = 0; i < cl->nKeys; ++i) {
    chrom.push_back(std::string(cl->chrom[i]));
    if (cl->len[i] > static_cast<uint32_t>(std::numeric_limits<int>::max())) {
      length.push_back(NA_INTEGER);
    } else {
      length.push_back(static_cast<int>(cl->len[i]));
    }
  }

  return DataFrame::create(
    _["chrom"] = chrom,
    _["length"] = length
  );
}

DataFrame query_impl(const std::string &file,
                     const std::string &chrom,
                     int start,
                     int end) {
  if (start < 1) {
    stop("`start` must be >= 1 in 1-based closed coordinates.");
  }
  if (end < start) {
    stop("`end` must be >= `start`.");
  }

  BigWigSession bw(file);
  int64_t idx = bw.chrom_index(chrom);
  if (idx < 0) {
    stop("Chromosome `%s` was not found in the bigWig file. Use seqinfo_bwg() to inspect available chromosome names.", chrom.c_str());
  }

  uint32_t chrom_len = bw.chroms()->len[idx];
  uint32_t c_start = static_cast<uint32_t>(start - 1);
  uint32_t c_end = static_cast<uint32_t>(end);

  if (c_start >= chrom_len) {
    return DataFrame::create(
      _["chrom"] = CharacterVector(),
      _["start"] = IntegerVector(),
      _["end"] = IntegerVector(),
      _["value"] = NumericVector()
    );
  }
  if (c_end > chrom_len) {
    c_end = chrom_len;
  }
  if (c_start >= c_end) {
    return DataFrame::create(
      _["chrom"] = CharacterVector(),
      _["start"] = IntegerVector(),
      _["end"] = IntegerVector(),
      _["value"] = NumericVector()
    );
  }

  bwOverlappingIntervals_t *intervals = bwGetOverlappingIntervals(
    bw.get(), chrom.c_str(), c_start, c_end
  );

  if (intervals == nullptr || intervals->l == 0) {
    if (intervals != nullptr) {
      bwDestroyOverlappingIntervals(intervals);
    }
    return DataFrame::create(
      _["chrom"] = CharacterVector(),
      _["start"] = IntegerVector(),
      _["end"] = IntegerVector(),
      _["value"] = NumericVector()
    );
  }

  std::vector<std::string> out_chrom;
  std::vector<int> out_start;
  std::vector<int> out_end;
  std::vector<double> out_value;

  out_chrom.reserve(intervals->l);
  out_start.reserve(intervals->l);
  out_end.reserve(intervals->l);
  out_value.reserve(intervals->l);

  for (uint32_t i = 0; i < intervals->l; ++i) {
    uint32_t s = intervals->start[i] + 1;
    uint32_t e = intervals->end[i];

    if (s < static_cast<uint32_t>(start)) {
      s = static_cast<uint32_t>(start);
    }
    if (e > static_cast<uint32_t>(end)) {
      e = static_cast<uint32_t>(end);
    }
    if (s > e) {
      continue;
    }

    out_chrom.push_back(chrom);
    out_start.push_back(static_cast<int>(s));
    out_end.push_back(static_cast<int>(e));
    out_value.push_back(static_cast<double>(intervals->value[i]));
  }

  bwDestroyOverlappingIntervals(intervals);

  return DataFrame::create(
    _["chrom"] = out_chrom,
    _["start"] = out_start,
    _["end"] = out_end,
    _["value"] = out_value
  );
}

}  // namespace

extern "C" SEXP _GeneTrackR_bw_seqinfo_cpp(SEXP fileSEXP) {
  BEGIN_RCPP
  std::string file = scalar_string(fileSEXP, "file");
  return Rcpp::wrap(seqinfo_impl(file));
  END_RCPP
}

extern "C" SEXP _GeneTrackR_bw_query_cpp(SEXP fileSEXP, SEXP chromSEXP, SEXP startSEXP, SEXP endSEXP) {
  BEGIN_RCPP
  std::string file = scalar_string(fileSEXP, "file");
  std::string chrom = scalar_string(chromSEXP, "chrom");
  int start = scalar_int(startSEXP, "start");
  int end = scalar_int(endSEXP, "end");
  return Rcpp::wrap(query_impl(file, chrom, start, end));
  END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
  {"_GeneTrackR_bw_seqinfo_cpp", (DL_FUNC) &_GeneTrackR_bw_seqinfo_cpp, 1},
  {"_GeneTrackR_bw_query_cpp", (DL_FUNC) &_GeneTrackR_bw_query_cpp, 4},
  {NULL, NULL, 0}
};

extern "C" void R_init_GeneTrackR(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
