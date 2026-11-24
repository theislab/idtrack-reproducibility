# The R tools of family G (docs/families/g.md, "The tools"), asked as their documentation says, one call per line.
#
# Usage: Rscript r_tools.R <tool> <organism> <lines.json> <out.json>
#   tool: orgdb | limma | hgnchelper | seurat | ensdb
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
# - seurat (human): UpdateSymbolList(<symbol>) (HGNC REST through Seurat's own code), one symbol per call with a pause
#   of 0.4 s after it (Seurat asks HGNC two requests per symbol without a pause of its own: about 5 requests per second,
#   half HGNC's documented rate), then the updated symbol's Ensembl genes by org.Hs.eg.db.
# - ensdb: the EnsDb of Ensembl 113 from AnnotationHub (the latest EnsDb it publishes); its genes table
#   (genes(edb, columns = gene_id, gene_name, entrezid)) read for the line's id (task b), name (task c) or NCBI gene
#   id (task e), as ensembldb's GeneIdFilter, GeneNameFilter and EntrezFilter select them.
suppressPackageStartupMessages({library(jsonlite); library(AnnotationDbi)})
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
  q <- q[q$species == latin & grepl("Ensembl 113 EnsDb", q$title, fixed = TRUE)]
  if (length(q) < 1) stop("no EnsDb of Ensembl 113 for ", latin)
  id <- names(q)[1]
  edb <- ah[[id]]
  ensdb_record <- list(ah_id = id, title = q$title[1], file = unname(AnnotationHub::cache(ah[id])),
                       rdatadateadded = as.character(q$rdatadateadded[1]), metadata = as.data.frame(metadata(edb)))
  edb_genes <- as.data.frame(genes(edb, columns = c("gene_id", "gene_name", "entrezid"), return.type = "data.frame"))
}
if (tool == "seurat") suppressPackageStartupMessages(library(Seurat))
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
    } else if (tool == "seurat") {
      fn <- "UpdateSymbolList+mapIds(SYMBOL->ENSEMBL)"
      syms <- suppressWarnings(UpdateSymbolList(symbols = l$id, verbose = FALSE))
      res <- ens_of(syms, "SYMBOL")
      Sys.sleep(0.4)
    } else if (tool == "ensdb") {
      g <- edb_genes
      if (l$task == "b") { fn <- "genes(GeneIdFilter)"; hit <- g$gene_id == l$id }
      else if (l$task == "c") { fn <- "genes(GeneNameFilter)"; hit <- !is.na(g$gene_name) & g$gene_name == l$id }
      else { fn <- "genes(EntrezFilter)"; hit <- vapply(g$entrezid, function(x) l$id %in% as.character(unlist(x)), logical(1)) }
      res <- sort(unique(g$gene_id[hit]))
    }
  }, error = function(e) err <<- conditionMessage(e))
  calls[[i]] <- list(qid = l$qid, `function` = fn, argument = l$id, result = I(res), symbols = I(syms),
                     error = err, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
meta <- if (is.null(db)) NULL else tryCatch(as.data.frame(metadata(db)), error = function(e) NULL)
ver <- if (is.null(db)) "" else as.character(packageVersion(pkg))
if (tool == "seurat") ver <- paste0("Seurat ", packageVersion("Seurat"), "; ", pkg, " ", ver)
if (tool == "ensdb") { ver <- paste0("ensembldb ", packageVersion("ensembldb"), "; ", ensdb_record$ah_id); meta <- ensdb_record }
if (tool == "limma") ver <- paste0("limma ", packageVersion("limma"), "; ", pkg, " ", ver)
if (tool == "hgnchelper") ver <- paste0("HGNChelper ", packageVersion("HGNChelper"), "; ", pkg, " ", ver)
write_json(list(tool = tool, version = ver, package_metadata = meta, calls = calls), out, auto_unbox = TRUE,
           null = "null", digits = NA)
