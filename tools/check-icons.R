#!/usr/bin/env Rscript
#### Are the icon names real? ####
#
# Every icon() and submitIcon in R/ is checked against the Font Awesome set that
# shiny actually ships. A name that is not in it renders as nothing and warns once
# at start-up, which is easy to miss in a log -- and the names that catch people
# out are the plausible ones, because Font Awesome Pro has icons the free set does
# not ("chart-diagram" among them).
#
# Needs R and the fontawesome package, so it is not part of rcheck.py, which runs
# where R is not installed.
#
#     Rscript tools/check-icons.R        # from the version directory
#
# Exits 1 if any name is unavailable, so it works as a pre-commit hook.

root <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(root)) root[1] else "shiny-server/R"
if (!dir.exists(dir)) stop("No such directory: ", dir, call. = FALSE)

txt <- unlist(lapply(list.files(dir, pattern = "\\.R$", full.names = TRUE),
                     readLines, warn = FALSE))
pattern <- 'icon\\("[^"]+"\\)|submitIcon = "[^"]+"'
names <- sort(unique(gsub('.*"([^"]+)".*', "\\1",
                          unlist(regmatches(txt, gregexpr(pattern, txt))))))

if (length(names) == 0) {
  cat("No icons found under", dir, "\n")
  quit(status = 0)
}

available <- vapply(names, function(n) {
  isTRUE(tryCatch({ fontawesome::fa(n); TRUE }, error = function(e) FALSE))
}, logical(1))

cat(length(names), "icon name(s) in", dir, "\n")
if (all(available)) {
  cat("OK - every one is in the Font Awesome set shiny ships.\n")
  quit(status = 0)
}
cat("\nNOT AVAILABLE:\n")
for (n in names[!available]) cat("  -", n, "\n")
cat("\nThese render as nothing. Check the free set at https://fontawesome.com/search?o=r&m=free\n")
quit(status = 1)
