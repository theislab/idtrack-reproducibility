# The R tools of family G (docs/families/g.md, "The tools"), asked as their documentation says, one call per line.
#
# Usage: Rscript r_tools.R <tool> <organism> <lines.json> <out.json>
#   tool: orgdb | limma | hgnchelper
#   lines.json: [{"qid": ..., "id": ..., "task": "c"|"d"|"e", "kind": "symbol"|"ncbi"|"refseq"|"uniprot"}, ...]
#   out.json:   {"tool": ..., "version": ..., "package_metadata": ..., "calls": [{"qid", "function", "argument",
#                "result": [Ensembl gene ids], "symbols": [...], "error", "seconds"}, ...]}
#
# - orgdb: AnnotationDbi::mapIds(<org package>, keys = <id>, keytype = SYMBOL (task c), ALIAS (task d), ENTREZID,
#   REFSEQ or UNIPROT (task e), column = "ENSEMBL", multiVals = "list").
# - limma: alias2Symbol(<id>, species = <code>) (aliases and symbols to official symbols, from the org package), then
#   the official symbols' Ensembl genes by the org package (mapIds SYMBOL -> ENSEMBL): limma gives symbols, the join is
#   the org package's.
# - hgnchelper (human and mouse): checkGeneSymbols(<id>, species = "human"|"mouse") with the map bundled in the
#   package, then the suggested symbols' Ensembl genes by the org package.
suppressPackageStartupMessages({library(jsonlite); library(AnnotationDbi)})
args <- commandArgs(trailingOnly = TRUE)
tool <- args[1]; organism <- args[2]; lines <- fromJSON(args[3], simplifyVector = FALSE); out <- args[4]
pkg <- switch(organism, homo_sapiens = "org.Hs.eg.db", mus_musculus = "org.Mm.eg.db", sus_scrofa = "org.Ss.eg.db",
              macaca_mulatta = "org.Mmu.eg.db", NA)
code <- switch(organism, homo_sapiens = "Hs", mus_musculus = "Mm", sus_scrofa = "Ss", macaca_mulatta = "Mmu", NA)
if (is.na(pkg)) stop("no org package for ", organism)
suppressPackageStartupMessages(library(pkg, character.only = TRUE))
db <- get(pkg)
keytype_of <- function(l) {
  if (l$task == "c") return("SYMBOL")
  if (l$task == "d") return("ALIAS")
  switch(l$kind, ncbi = "ENTREZID", refseq = "REFSEQ", uniprot = "UNIPROT", NA)
}
ens_of <- function(keys, keytype) {
  if (length(keys) == 0) return(character(0))
  r <- tryCatch(mapIds(db, keys = keys, keytype = keytype, column = "ENSEMBL", multiVals = "list"),
                error = function(e) NULL)
  if (is.null(r)) return(character(0))
  v <- unlist(r, use.names = FALSE)
  sort(unique(v[!is.na(v)]))
}
if (tool == "limma") suppressPackageStartupMessages(library(limma))
if (tool == "hgnchelper") suppressPackageStartupMessages(library(HGNChelper))
calls <- vector("list", length(lines))
for (i in seq_along(lines)) {
  l <- lines[[i]]
  t0 <- Sys.time()
  res <- character(0); syms <- character(0); err <- NULL; fn <- ""
  tryCatch({
    if (tool == "orgdb") {
      kt <- keytype_of(l); fn <- paste0("mapIds(", kt, "->ENSEMBL)")
      v <- if (kt == "REFSEQ") sub("\\.[0-9]+$", "", l$id) else l$id
      res <- ens_of(v, kt)
    } else if (tool == "limma") {
      fn <- "alias2Symbol+mapIds(SYMBOL->ENSEMBL)"
      syms <- alias2Symbol(l$id, species = code, expand.symbols = FALSE)
      res <- ens_of(syms, "SYMBOL")
    } else if (tool == "hgnchelper") {
      fn <- "checkGeneSymbols+mapIds(SYMBOL->ENSEMBL)"
      sp <- if (organism == "homo_sapiens") "human" else "mouse"
      ck <- suppressWarnings(checkGeneSymbols(l$id, species = sp))
      s <- ck$Suggested.Symbol[1]
      syms <- if (is.na(s)) character(0) else trimws(strsplit(s, "///", fixed = TRUE)[[1]])
      res <- ens_of(syms, "SYMBOL")
    }
  }, error = function(e) err <<- conditionMessage(e))
  calls[[i]] <- list(qid = l$qid, `function` = fn, argument = l$id, result = I(res), symbols = I(syms),
                     error = err, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
meta <- tryCatch(as.data.frame(metadata(db)), error = function(e) NULL)
ver <- as.character(packageVersion(pkg))
if (tool == "limma") ver <- paste0("limma ", packageVersion("limma"), "; ", pkg, " ", ver)
if (tool == "hgnchelper") ver <- paste0("HGNChelper ", packageVersion("HGNChelper"), "; ", pkg, " ", ver)
write_json(list(tool = tool, version = ver, package_metadata = meta, calls = calls), out, auto_unbox = TRUE,
           null = "null", digits = NA)
