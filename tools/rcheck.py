#!/usr/bin/env python3
"""Static checks for the CRAVE R sources.

R is not always available where CRAVE is edited, and there is no R linter in the
deployment image, so this stands in for one. It is deliberately not an R parser: it
strips comments and string literals and then works on the remaining text. That is
enough to catch the classes of mistake that have actually occurred in this codebase,
and it will never catch a type error or anything that depends on runtime values.

Usage:
    python3 tools/rcheck.py [--root shiny-server] [--quiet]

Exit status is 1 if any check marked [HARD] fails, otherwise 0, so it can be wired
into a pre-commit hook or CI step. Checks marked [ADVISORY] print findings for a
human to read and never affect the exit status: they have known false positives,
which are explained where they occur.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ROOT = os.path.normpath(os.path.join(HERE, "..", "shiny-server"))

failures: list[str] = []


# ---------------------------------------------------------------------------
# Lexing helpers
# ---------------------------------------------------------------------------

def strip_all(src: str) -> str:
    """Blank comments and string literals, preserving offsets and line numbers.

    Used by the checks that care about code structure, so that a brace inside a
    string or a `#` inside a URL cannot confuse them.
    """
    out: list[str] = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c == "#":
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'`":
            q = c
            out.append(" ")
            i += 1
            while i < n:
                if src[i] == "\\":
                    out.append("  ")
                    i += 2
                    continue
                if src[i] == q:
                    out.append(" ")
                    i += 1
                    break
                out.append("\n" if src[i] == "\n" else " ")
                i += 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def strip_comments(src: str) -> str:
    """Blank comments but keep string literals, for checks that read id names."""
    out: list[str] = []
    i, n = 0, len(src)
    while i < n:
        c = src[i]
        if c == "#":
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'":
            q = c
            out.append(c)
            i += 1
            while i < n:
                if src[i] == "\\":
                    out.append(src[i:i + 2])
                    i += 2
                    continue
                out.append(src[i])
                if src[i] == q:
                    i += 1
                    break
                i += 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def line_of(src: str, pos: int) -> int:
    return src[:pos].count("\n") + 1


def split_top_level(s: str) -> list[str]:
    """Split an argument list on commas that are not nested inside brackets."""
    parts, depth, cur, q = [], 0, [], None
    for ch in s:
        if q:
            cur.append(ch)
            if ch == q:
                q = None
            continue
        if ch in "\"'`":
            q = ch
            cur.append(ch)
            continue
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


def call_body(src: str, open_paren: int) -> tuple[str, int]:
    """Return the text inside a call whose '(' is at `open_paren`, and the index
    just past its ')'."""
    depth, i = 1, open_paren + 1
    while i < len(src) and depth:
        if src[i] in "([{":
            depth += 1
        elif src[i] in ")]}":
            depth -= 1
        i += 1
    return src[open_paren + 1:i - 1], i


def banner(title: str, kind: str) -> None:
    print()
    print("=" * 78)
    print(f"{title}  [{kind}]")
    print("=" * 78)


# ---------------------------------------------------------------------------
# 1. Delimiter balance
# ---------------------------------------------------------------------------

def check_balance(files, codes) -> None:
    banner("1. DELIMITER BALANCE", "HARD")
    pairs = {")": "(", "]": "[", "}": "{"}
    problems = []
    for f, rel in files:
        stack, line = [], 1
        for ch in codes[f]:
            if ch == "\n":
                line += 1
            elif ch in "([{":
                stack.append((ch, line))
            elif ch in ")]}":
                if not stack:
                    problems.append(f"  {rel}:{line}: unmatched closing '{ch}'")
                elif stack[-1][0] != pairs[ch]:
                    o, ol = stack.pop()
                    problems.append(
                        f"  {rel}:{line}: '{ch}' closes '{o}' opened at line {ol}")
                else:
                    stack.pop()
        for o, ol in stack:
            problems.append(f"  {rel}:{ol}: '{o}' never closed")
    if problems:
        print("\n".join(problems))
        failures.append("delimiter balance")
    else:
        print("OK - every file is balanced.")


# ---------------------------------------------------------------------------
# 2. Undefined function calls
# ---------------------------------------------------------------------------

# Names expected to come from base R or from an attached package. Extend this when
# you start using a new function; an entry appearing in the output only means this
# list has not caught up, not necessarily that anything is wrong.
KNOWN = set("""
if else for while repeat function return switch tryCatch try stop warning message
c list vector character numeric integer logical double unlist lapply sapply vapply
mapply Map Filter Reduce do.call Recall Negate identity
length names setNames seq seq_len seq_along rep rev sort order unique duplicated
anyDuplicated match nchar nzchar rank findInterval which
paste paste0 sprintf format formatC gsub sub grepl grep regmatches gregexpr regexpr
strsplit substr substring trimws toupper tolower startsWith endsWith sQuote dQuote
is.null is.na is.nan is.finite is.function is.list is.numeric is.character
is.logical is.data.frame is.matrix is.environment isTRUE isFALSE identical inherits
as.character as.numeric as.integer as.logical as.list as.matrix as.data.frame
as.vector as.dendrogram as.name as.symbol
abs sign max min sum mean median var sd round signif ceiling floor trunc exp log
log2 log10 sqrt cumsum cumprod pmax pmin range prod
rowMeans colMeans rowSums colSums crossprod outer diag t dim nrow ncol rownames
colnames dimnames cbind rbind matrix array apply sweep scale
head tail intersect setdiff union setequal anyNA complete.cases
data.frame factor levels nlevels droplevels split unsplit table tabulate
file.exists file.path file.info file.access file.remove file.copy dir.exists
dir.create basename dirname normalizePath list.files
readRDS saveRDS readLines writeLines Sys.glob Sys.time Sys.Date Sys.setenv Sys.which
Sys.getenv Sys.getpid system system2 shQuote unlink tempdir tempfile
new.env get get0 exists ls assign rm local on.exit sys.source source library
require requireNamespace loadNamespace attr attributes structure unname
suppressWarnings suppressMessages suppressPackageStartupMessages invisible force
missing nargs all any ifelse stopifnot Vectorize
options getOption conditionMessage conditionCall simpleError simpleWarning
system.file packageVersion R.home nchar merge reshape aggregate by
. .N .SD .I .GRP .BY
set.seed sample runif rnorm p.adjust phyper quantile approx approxfun prcomp cor
difftime as.difftime emptyenv globalenv environment parent.frame match.arg nlevels
.packages utils
tibble as_tibble bind_rows bind_cols mutate transmute select filter arrange
group_by ungroup summarise summarize left_join inner_join full_join anti_join
semi_join distinct relocate rename pull slice slice_head slice_max slice_min
across all_of any_of everything starts_with ends_with contains matches
if_else case_when case_match n row_number first last nth coalesce na_if recode
pick reframe rowwise desc
data.table as.data.table is.data.table rbindlist setorder setorderv setnames
setcolorder fwrite fread CJ setDT copy uniqueN fifelse
ggplot aes geom_point geom_line geom_tile geom_text geom_bar geom_col geom_violin
geom_boxplot geom_abline geom_hline geom_vline geom_label geom_segment geom_smooth
scale_fill_gradient scale_fill_gradient2 scale_fill_gradientn scale_fill_manual
scale_colour_gradient scale_colour_discrete scale_colour_manual scale_color_manual
scale_x_discrete scale_y_discrete scale_x_continuous scale_y_continuous
theme theme_classic theme_minimal theme_bw element_text element_blank element_rect
element_line xlab ylab ggtitle labs annotate coord_equal coord_flip facet_wrap
ggsave rescale trans_new log_breaks
plot_ly layout ggplotly renderPlotly plotlyOutput add_trace add_markers subplot
visNetwork visNetworkOutput renderVisNetwork visOptions visIgraph visNodes visEdges
renderDT DTOutput dataTableProxy selectRows updateSearch datatable formatRound
shinyApp fluidPage fluidRow column navbarPage tabPanel tabsetPanel navlistPanel
sidebarLayout sidebarPanel mainPanel wellPanel conditionalPanel
tags tagList HTML h1 h2 h3 h4 h5 h6 hr br icon div p span img a strong em
actionButton actionLink downloadButton downloadHandler fileInput
textInput textAreaInput numericInput checkboxInput checkboxGroupInput radioButtons
selectizeInput selectInput sliderInput dateInput
updateTextInput updateTextAreaInput updateNumericInput updateCheckboxInput
updateCheckboxGroupInput updateRadioButtons updateSelectizeInput updateSelectInput
updateSliderInput updateSliderTextInput updateTabsetPanel updateNavbarPage
uiOutput renderUI textOutput renderText plotOutput renderPlot verbatimTextOutput
renderPrint tableOutput renderTable
reactive reactiveVal reactiveValues reactiveValuesToList observe observeEvent
eventReactive isolate req validate need bindCache bindEvent
debounce throttle moduleServer NS showModal modalDialog modalButton removeModal
showNotification removeNotification withProgress incProgress setProgress Progress
outputOptions onStop onSessionEnded onFlush getDefaultReactiveDomain session
shinytheme busy_start_up spin_epic add_busy_bar add_busy_spinner
show_modal_spinner remove_modal_spinner materialSwitch updateMaterialSwitch
prettySwitch updatePrettySwitch
colourInput updateColourInput
log_info log_warn log_error log_debug log_trace log_threshold
log_shiny_input_changes
foreach registerDoSEQ
agnes clara pam diana ggdendrogram dendro_data
graph_from_data_frame V E vcount ecount
umap ggVennDiagram scale_x_upset axis_combmatrix AUC complete mice Rtsne
htmlEscape saveWidget
""".split())

# Assignments whose value is callable, so `x()` later is not an undefined function.
CLOSURE_RHS = (
    "function", "reactive", "reactiveVal", "eventReactive", "debounce", "throttle",
    "local", "dataTableProxy", "Negate", "approxfun", "downloadHandler",
)


def check_undefined_calls(files, codes, raws) -> None:
    banner("2. CALLS TO NAMES THAT ARE NEITHER DEFINED HERE NOR KNOWN", "ADVISORY")
    print("A name here is either a package function missing from this script's")
    print("allowlist, or a genuine typo. Check before dismissing.")
    print("Anything bound to a name anywhere in the sources counts as defined,")
    print("since R cannot be asked whether that binding holds a function -- so a")
    print("typo that happens to match a variable name will not be reported.")

    defined = set()
    for f, _ in files:
        code = codes[f]
        # Any binding at all: a name assigned somewhere may legitimately hold a
        # function, e.g. `myGeom <- if (boxplot) geom_boxplot else geom_violin`.
        for m in re.finditer(
                r"(?:^|[;{(,]|<-)\s*([A-Za-z._][A-Za-z0-9._]*)\s*(?:<-|<<-)(?!-)",
                code, re.M):
            defined.add(m.group(1))
        # Infix operators are declared with backticks, which strip_all() blanks.
        for m in re.finditer(r"^\s*`([^`]+)`\s*<-\s*function\s*\(", raws[f], re.M):
            defined.add(m.group(1))
        # Function arguments can be callables supplied by the caller.
        for m in re.finditer(r"function\s*\(([^)]*)\)", code):
            for arg in split_top_level(m.group(1)):
                nm = re.match(r"^([A-Za-z._][A-Za-z0-9._]*)", arg)
                if nm:
                    defined.add(nm.group(1))

    calls = defaultdict(set)
    for f, rel in files:
        for m in re.finditer(r"(?<![$@:.\w])([A-Za-z._][A-Za-z0-9._]*)\s*\(", codes[f]):
            name = m.group(1)
            if name in defined or name in KNOWN:
                continue
            calls[name].add(f"{rel}:{line_of(codes[f], m.start())}")

    if not calls:
        print("OK - every call resolves.")
        return
    for name in sorted(calls):
        where = sorted(calls[name])
        extra = f"  (+{len(where) - 1} more)" if len(where) > 1 else ""
        print(f"  {name:36s} {where[0]}{extra}")


# ---------------------------------------------------------------------------
# 3. Module input/output id cross-reference
# ---------------------------------------------------------------------------

def check_module_ids(root, files, raws) -> None:
    banner("3. SHINY MODULE INPUT/OUTPUT ID CROSS-REFERENCE", "ADVISORY")
    print("Ids built at runtime, e.g. ns(paste0(prefix, '_dist')), look 'missing'")
    print("here. Sub-module ids passed to plotPanelUI or customiseUI look 'unread'.")

    mod_files = [(f, rel) for f, rel in files
                 if re.search(r"R/1[0-8]_mod_", rel.replace(os.sep, "/"))]
    for f, rel in mod_files:
        code = strip_comments(raws[f])
        ui_ids = set(re.findall(r"\bns\(\s*\"([^\"]+)\"\s*\)", code))
        read = set(re.findall(r"\binput\$([A-Za-z._][A-Za-z0-9._]*)", code))
        read |= set(re.findall(r"input\[\[\s*\"([^\"]+)\"", code))
        upd = set(re.findall(
            r"\bupdate[A-Za-z]+\(\s*session\s*,\s*\"([^\"]+)\"", code))
        outs = set(re.findall(r"\boutput\$([A-Za-z._][A-Za-z0-9._]*)", code))
        proxies = set(re.findall(r"\bdataTableProxy\(\s*\"([^\"]+)\"\s*\)", code))

        # DT synthesises these from a table output's id.
        derived = set()
        for t in list(outs) + list(proxies):
            derived |= {t + s for s in ("_rows_selected", "_rows_all",
                                        "_rows_current", "_search_columns",
                                        "_state", "_cell_clicked")}

        referenced = read | upd
        missing = sorted(referenced - ui_ids - derived)
        unread = sorted(ui_ids - referenced - outs)
        print(f"\n-- {rel}")
        print(f"   ns() ids: {len(ui_ids)}   input$: {len(read)}   "
              f"update*: {len(upd)}   output$: {len(outs)}")
        if missing:
            print(f"   referenced, no ns() literal: {missing}")
        if unread:
            print(f"   ns() literal, never read   : {unread}")
        if not missing and not unread:
            print("   OK")


# ---------------------------------------------------------------------------
# 4. plotPanel UI/server pairing
# ---------------------------------------------------------------------------

def check_plotpanel_pairs(files, raws) -> None:
    banner("4. plotPanelUI / plotPanelServer PAIRING", "HARD")
    ui, srv = set(), set()
    for f, _ in files:
        s = raws[f]
        ui |= set(re.findall(r"plotPanelUI\(\s*ns\(\"([A-Za-z0-9_]+)\"\)", s))
        srv |= set(re.findall(r"plotPanelServer\(\s*\"([A-Za-z0-9_]+)\"", s))
    print(f"UI panels    : {len(ui)}")
    print(f"Server panels: {len(srv)}")
    if ui == srv:
        print("OK - every panel has both halves.")
        return
    if ui - srv:
        print(f"  UI without a server: {sorted(ui - srv)}")
    if srv - ui:
        print(f"  server without a UI: {sorted(srv - ui)}")
    failures.append("plotPanel pairing")


# ---------------------------------------------------------------------------
# 5. dplyr shadow-then-read
# ---------------------------------------------------------------------------

DPLYR_VERBS = ("summarise", "summarize", "mutate", "transmute")


def check_dplyr_shadowing(files, codes) -> None:
    banner("5. COLUMNS CREATED AND THEN RE-READ IN THE SAME dplyr CALL", "ADVISORY")
    print("dplyr evaluates these arguments in order and a new column shadows any")
    print("existing one of the same name, so a later argument sees the NEW value.")
    print("Chaining like this is legal and often intended. It is a bug only where")
    print("the later argument meant to read the ORIGINAL column -- which is what")
    print("made Reduce's Others column empty in 5.1.1. Read each one and decide.")

    hits = []
    for f, rel in files:
        code = codes[f]
        for verb in DPLYR_VERBS:
            for m in re.finditer(r"\b" + verb + r"\s*\(", code):
                args, _ = call_body(code, m.end() - 1)
                assigned: list[str] = []
                for a in split_top_level(args):
                    am = re.match(
                        r"^(`[^`]+`|[A-Za-z._][A-Za-z0-9._]*)\s*=(?!=)\s*(.*)$",
                        a, re.S)
                    name = am.group(1).strip("`") if am else None
                    rhs = am.group(2) if am else a
                    for prev in assigned:
                        plain = re.match(r"^[A-Za-z._][A-Za-z0-9._]*$", prev)
                        pat = (r"\b" + re.escape(prev) + r"\b") if plain \
                            else (r"`" + re.escape(prev) + r"`")
                        if re.search(pat, rhs):
                            hits.append(
                                f"  {rel}:{line_of(code, m.start())}  {verb}(): "
                                f"'{prev}' is assigned earlier, then read by "
                                f"'{name or '<unnamed>'}'")
                    if name:
                        assigned.append(name)
    print()
    print("\n".join(dict.fromkeys(hits)) if hits else "  none found")


# ---------------------------------------------------------------------------
# 6. Package manifests vs the Dockerfile
# ---------------------------------------------------------------------------

MANIFESTS = ("CRAVE_PKGS_EAGER", "CRAVE_PKGS_LAZY", "CRAVE_PKGS_INSTALLED_ONLY")


def check_manifests(root) -> None:
    banner("6. R PACKAGE MANIFESTS vs THE DOCKERFILE", "HARD")
    pkg_file = os.path.join(root, "R", "00_packages.R")
    dockerfile = os.path.normpath(os.path.join(root, "..", "Dockerfile"))
    if not (os.path.exists(pkg_file) and os.path.exists(dockerfile)):
        print("SKIP - 00_packages.R or Dockerfile not found.")
        return

    src = open(pkg_file, encoding="utf-8").read()
    docker = open(dockerfile, encoding="utf-8").read()

    def manifest(name):
        m = re.search(name + r"\s*<-\s*c\((.*?)\)\s*\n", src, re.S)
        return sorted(set(re.findall(r'"([^"]+)"', m.group(1)))) if m else None

    blocks = [sorted(set(re.findall(r"'([^']+)'", b)))
              for b in re.findall(r"install\.packages\((.*?)repos", docker, re.S)]

    if len(blocks) != len(MANIFESTS):
        print(f"  Dockerfile has {len(blocks)} install.packages() blocks, "
              f"expected {len(MANIFESTS)} (one per manifest).")
        failures.append("manifest/Dockerfile block count")
        return

    ok = True
    for name, block in zip(MANIFESTS, blocks):
        declared = manifest(name)
        if declared is None:
            print(f"  {name}: not found in 00_packages.R")
            ok = False
            continue
        if declared == block:
            print(f"  {name:26s} MATCH ({len(declared)})")
        else:
            ok = False
            print(f"  {name:26s} MISMATCH")
            print(f"     only in 00_packages.R: {sorted(set(declared) - set(block))}")
            print(f"     only in Dockerfile   : {sorted(set(block) - set(declared))}")
    if not ok:
        failures.append("package manifests")


# ---------------------------------------------------------------------------
# 7. Sizes
# ---------------------------------------------------------------------------

def report_sizes(files, raws) -> None:
    banner("7. FILE SIZES", "INFO")
    total = 0
    for f, rel in files:
        n = raws[f].count("\n") + 1
        total += n
        print(f"  {n:5d}  {rel}")
    print(f"  {total:5d}  TOTAL")


# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT,
                    help="app directory containing app.R and R/ "
                         "(default: ../shiny-server relative to this script)")
    ap.add_argument("--quiet", action="store_true",
                    help="omit the advisory and info sections")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    app = os.path.join(root, "app.R")
    rdir = os.path.join(root, "R")
    if not os.path.isfile(app) or not os.path.isdir(rdir):
        print(f"error: {root} does not look like a CRAVE app directory "
              f"(expected app.R and R/).", file=sys.stderr)
        return 2

    paths = [app] + sorted(
        os.path.join(rdir, p) for p in os.listdir(rdir) if p.endswith(".R"))
    files = [(p, os.path.relpath(p, root)) for p in paths]
    raws, codes = {}, {}
    for p, _ in files:
        raws[p] = open(p, encoding="utf-8").read()
        codes[p] = strip_all(raws[p])

    print(f"CRAVE static checks - {root}")
    print(f"{len(files)} R files")

    check_balance(files, codes)
    check_plotpanel_pairs(files, raws)
    check_manifests(root)
    if not args.quiet:
        check_undefined_calls(files, codes, raws)
        check_module_ids(root, files, raws)
        check_dplyr_shadowing(files, codes)
        report_sizes(files, raws)

    print()
    print("=" * 78)
    if failures:
        print("FAILED: " + ", ".join(failures))
        return 1
    print("All HARD checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
