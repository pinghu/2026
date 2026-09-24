#!/usr/bin/env Rscript
# Remove columns whose header matches an id in bad.list.
# Usage: Rscript remove_bad_cols.R <data_file> <bad_list> [out_file]
a <- commandArgs(TRUE)
if (length(a) < 2) stop("Usage: remove_bad_cols.R <data_file> <bad_list> [out_file]")
out <- if (length(a) >= 3) a[3] else sub("(\\.[^.]*)?$", "_clean\\1", a[1])

d <- read.delim(a[1], header = TRUE, sep = "\t", check.names = FALSE,
                quote = "", comment.char = "")
bad <- trimws(readLines(a[2], warn = FALSE))
bad <- bad[bad != "" & bad != "bad.id"]           # drop header/blank lines

keep <- !(colnames(d) %in% bad)
message(sprintf("Removed %d of %d columns", sum(!keep), ncol(d)))
write.table(d[, keep, drop = FALSE], out, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
message("Wrote ", out)
