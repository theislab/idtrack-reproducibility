# The R tools of family G (docs/families/g.md, "The tools"), asked as their documentation says.
#
# Usage: Rscript r_tools.R <tool> <organism> <lines.json> <out.json>
#   tool: orgdb | limma | hgnchelper | seurat | ensdb
#   lines.json: [{"qid": ..., "id": ..., "task": "c"|"d"|"e", "kind": "symbol"|"ncbi"|"refseq"|"uniprot"}, ...]
#   out.json:   {"tool": ..., "version": ..., "package_metadata": ..., "calls": [{"qid", "function", "argument",
#                "result": [Ensembl gene ids], "symbols": [...], "error", "seconds"}, ...]}
#
# - orgdb: AnnotationDbi::mapIds(<org package>, keys = <ids>, keytype = SYMBOL (task c), ALIAS (task d), ENTREZID,
#   REFSEQ or UNIPROT (task e), column = "ENSEMBL", multiVals = "list").
# - limma: alias2SymbolTable(<ids>, species = <code>) (each alias or symbol to its official symbol, from the org
#   package; the one of the lowest Entrez id where several), then the official symbols' Ensembl genes by the org
#   package (mapIds SYMBOL -> ENSEMBL): limma gives symbols, the join is the org package's.
# - hgnchelper (human and mouse): checkGeneSymbols(<ids>, species = "human"|"mouse") with the map bundled in the
#   package, then the suggested symbols' Ensembl genes by the org package.
# - seurat (human): UpdateSymbolList(<symbol>) (HGNC REST through Seurat's own code), one symbol per call with a pause
#   of 0.4 s after it (Seurat asks HGNC two requests per symbol without a pause of its own: about 5 requests per second,
#   half HGNC's documented rate), then the updated symbol's Ensembl genes by org.Hs.eg.db.
# - ensdb: the EnsDb of Ensembl 113 from AnnotationHub (the latest EnsDb it publishes); its genes table
#   (genes(edb, columns = gene_id, gene_name, entrezid)) read for the line's id (task b), name (task c) or NCBI gene
#   id (task e), as ensembldb's GeneIdFilter, GeneNameFilter and EntrezFilter select them.
suppressPackageStartupMessages({library(jsonlite); library(AnnotationDbi)})
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE)
tool <- args[1]; organism <- args[2]; lines <- fromJSON(args[3], simplifyVector = FALSE); out <- args[4]
pkg <- switch(organism, homo_sapiens = "org.Hs.eg.db", mus_musculus = "org.Mm.eg.db", sus_scrofa = "org.Ss.eg.db",
              macaca_mulatta = "org.Mmu.eg.db", NA)
code <- switch(organism, homo_sapiens = "Hs", mus_musculus = "Mm", sus_scrofa = "Ss", macaca_mulatta = "Mmu", NA)
db <- NULL
if (!is.na(pkg)) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  db <- get(pkg)
} else if (tool != "ensdb") stop("no org package for ", organism)
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
edb_genes <- NULL; ensdb_record <- NULL
if (tool == "ensdb") {
  suppressPackageStartupMessages({library(AnnotationHub); library(ensembldb)})
  latin <- switch(organism, homo_sapiens = "Homo sapiens", mus_musculus = "Mus musculus", sus_scrofa = "Sus scrofa",
                  macaca_mulatta = "Macaca mulatta", macaca_fascicularis = "Macaca fascicularis")
  ah <- AnnotationHub(ask = FALSE)
  q <- query(ah, c("EnsDb", latin, "113"))
  # the reference assembly's EnsDb (the hub also holds EnsDbs of a species' other assemblies, e.g. pig's breeds)
  assembly <- switch(organism, homo_sapiens = "GRCh38", mus_musculus = "GRCm39", sus_scrofa = "Sscrofa11.1",
                     macaca_mulatta = "Mmul_10", macaca_fascicularis = "Macaca_fascicularis_6.0")
  q <- q[q$species == latin & grepl("Ensembl 113 EnsDb", q$title, fixed = TRUE) & q$genome == assembly]
  if (length(q) < 1) stop("no EnsDb of Ensembl 113 for ", latin)
  id <- names(q)[1]
  edb <- ah[[id]]
  ensdb_record <- list(ah_id = id, title = q$title[1], genome = q$genome[1], file = unname(AnnotationHub::cache(ah[id])),
                       rdatadateadded = as.character(q$rdatadateadded[1]), metadata = as.data.frame(metadata(edb)))
  edb_genes <- as.data.frame(genes(edb, columns = c("gene_id", "gene_name", "entrezid"), return.type = "data.frame"))
}
if (tool == "seurat") suppressPackageStartupMessages(library(Seurat))
if (tool == "limma") suppressPackageStartupMessages(library(limma))
if (tool == "hgnchelper") suppressPackageStartupMessages(library(HGNChelper))
# The tools are asked as their vectorised calls take a list (one call per keytype or per tool), which is how they are
# documented for many keys; Seurat, which asks HGNC per symbol, one symbol per call with its pause.
ids <- vapply(lines, function(l) l$id, character(1))
t_all <- Sys.time()
lookup <- list(); fn_of <- function(l) ""
map_list <- function(keys, keytype) {
  keys <- unique(keys[!is.na(keys) & keys != ""])
  if (length(keys) == 0) return(list())
  r <- tryCatch(mapIds(db, keys = keys, keytype = keytype, column = "ENSEMBL", multiVals = "list"), error = function(e) NULL)
  if (is.null(r)) return(list())
  lapply(r, function(v) sort(unique(v[!is.na(v)])))
}
if (tool == "orgdb") {
  kts <- vapply(lines, keytype_of, character(1))
  for (kt in unique(kts[!is.na(kts)])) lookup[[kt]] <- map_list(ids[!is.na(kts) & kts == kt], kt)
} else if (tool == "limma") {
  sym <- suppressWarnings(alias2SymbolTable(ids, species = code))
  ens <- map_list(sym, "SYMBOL")
} else if (tool == "hgnchelper") {
  sp <- if (organism == "homo_sapiens") "human" else "mouse"
  ck <- suppressWarnings(checkGeneSymbols(ids, species = sp))
  sug <- lapply(ck$Suggested.Symbol, function(s) if (is.na(s)) character(0) else trimws(strsplit(s, "///", fixed = TRUE)[[1]]))
  ens <- map_list(unlist(sug), "SYMBOL")
} else if (tool == "ensdb") {
  g <- edb_genes
  by_id <- split(g$gene_id, g$gene_id)
  nm <- !is.na(g$gene_name) & g$gene_name != ""
  by_name <- split(g$gene_id[nm], g$gene_name[nm])
  ez <- lapply(g$entrezid, function(x) as.character(unlist(x)))
  by_entrez <- split(rep(g$gene_id, lengths(ez)), unlist(ez))
}
calls <- vector("list", length(lines))
for (i in seq_along(lines)) {
  l <- lines[[i]]
  t0 <- Sys.time()
  res <- character(0); syms <- character(0); err <- NULL; fn <- ""
  tryCatch({
    if (tool == "orgdb") {
      kt <- keytype_of(l); fn <- paste0("mapIds(", kt, "->ENSEMBL)")
      res <- if (is.na(kt)) character(0) else (lookup[[kt]][[l$id]] %||% character(0))
    } else if (tool == "limma") {
      fn <- "alias2SymbolTable+mapIds(SYMBOL->ENSEMBL)"
      syms <- if (is.na(sym[i])) character(0) else sym[i]
      res <- if (length(syms)) (ens[[syms]] %||% character(0)) else character(0)
    } else if (tool == "hgnchelper") {
      fn <- "checkGeneSymbols+mapIds(SYMBOL->ENSEMBL)"
      syms <- sug[[i]]
      res <- sort(unique(unlist(lapply(syms, function(x) ens[[x]]))))
    } else if (tool == "seurat") {
      fn <- "UpdateSymbolList+mapIds(SYMBOL->ENSEMBL)"
      syms <- suppressWarnings(UpdateSymbolList(symbols = l$id, verbose = FALSE))
      res <- unlist(map_list(syms, "SYMBOL"))
      Sys.sleep(0.4)
    } else if (tool == "ensdb") {
      if (l$task == "b") { fn <- "genes(GeneIdFilter)"; res <- by_id[[l$id]] }
      else if (l$task == "c") { fn <- "genes(GeneNameFilter)"; res <- by_name[[l$id]] }
      else { fn <- "genes(EntrezFilter)"; res <- by_entrez[[l$id]] }
      res <- sort(unique(as.character(res)))
    }
  }, error = function(e) err <<- conditionMessage(e))
  calls[[i]] <- list(qid = l$qid, `function` = fn, argument = l$id, result = I(as.character(res)), symbols = I(syms),
                     error = err, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
per_line <- as.numeric(difftime(Sys.time(), t_all, units = "secs")) / max(1, length(lines))
if (tool != "seurat") for (i in seq_along(calls)) calls[[i]]$seconds <- per_line  # a vectorised call's share
meta <- if (is.null(db)) NULL else tryCatch(as.data.frame(metadata(db)), error = function(e) NULL)
ver <- if (is.null(db)) "" else as.character(packageVersion(pkg))
if (tool == "seurat") ver <- paste0("Seurat ", packageVersion("Seurat"), "; ", pkg, " ", ver)
if (tool == "ensdb") { ver <- paste0("ensembldb ", packageVersion("ensembldb"), "; ", ensdb_record$ah_id); meta <- ensdb_record }
if (tool == "limma") ver <- paste0("limma ", packageVersion("limma"), "; ", pkg, " ", ver)
if (tool == "hgnchelper") ver <- paste0("HGNChelper ", packageVersion("HGNChelper"), "; ", pkg, " ", ver)
write_json(list(tool = tool, version = ver, package_metadata = meta, calls = calls), out, auto_unbox = TRUE,
           null = "null", digits = NA)
