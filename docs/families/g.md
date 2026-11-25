# G: comparison with other tools

The tools people use for the questions idtrack answers, asked the same lines and judged by the same truth and the same scorer as idtrack ([methods](../methods.md#g-comparison-with-other-tools)): their versions and the dates they were asked, their raw responses kept as they came, a question a tool cannot pose scored `not_supported`, and tables of tool by task.

## State

The comparison sets are published (`g_<family>_v1`, 2026-11-24) and the rule of equal truth is computed for releases 116, 114, 113 and 111; the passes of 2026-11-24 are being asked (one per service, one date per service). The pilot below was asked before the sets were published, at most 100 lines per tool and organism; its responses are kept with the others (endpoints without a form) and are not read by the adapters.

## What it measures

Per tool, task and organism: which lines it can pose (`not_supported` otherwise, with the reason), coverage (lines with an answer), accuracy (the rate of the common scorer, right / (right + wrong + error), with Wilson and cluster-bootstrap intervals) under the views `all` (the tool's answer is exactly the truth's genes) and `best` (a non-empty part of them, or nothing where the truth is nothing), the lines left out by the rule of equal truth and why, the accuracy by distance in releases for old ids, seconds per 1,000 lines with the network, the overlap of the tools' failures, and every disagreement with idtrack sampled and read.

## The tools

Each tool is asked as its own documentation says, with its best settings for the question; its version (or service date) is read from the tool at every pass and recorded with every response. Rates: half the documented rate, one connection per host over all jobs together; one request per second where no rate is documented. Every service below answers without a key or an account; nothing identifying the person who runs the benchmark is sent (a generic client name in the User-Agent, no contact address).

| tool | how it is asked | release it answers at | rate (documented; used) | terms of automated use | organisms | tasks |
| --- | --- | --- | --- | --- | --- | --- |
| idtrack | the tree measured on its graphs of 115 (the families' answers) | 115 | local | – | all five | a, b, c, d, e |
| Ensembl REST | `POST archive/id` (400 ids per call: "limit your requests to up to 400 IDs per batch", the service's wiki, POST requests), `GET xrefs/symbol/:species/:symbol?object_type=gene` (one per line: "a display name for a gene/transcript/translation, a synonym or an externally linked reference"), `POST lookup/symbol/:species` (1,000; display names only, reported apart) | 116 (`info/data`: `{"releases":[116]}`); the only release it serves | "55000 requests over an hour" / "15 requests per second" (github.com/Ensembl/ensembl-rest/wiki/Rate-Limits); 7.5 per second | "Your usage of the data returned by the service is subject to same conditions as laid out on the Ensembl website" (the wiki's home page); a service for programmatic use | all five | a, b, c, d, e |
| BioMart (martservice XML) | one XML query per filter and 500 values (`ensembl_gene_id`; `external_gene_name`; `external_synonym`; `entrezgene_id`; `refseq_mrna`, `refseq_mrna_predicted`, `refseq_ncrna`, `refseq_ncrna_predicted`, `refseq_peptide`, `refseq_peptide_predicted` by the accession's prefix; `uniprotswissprot` and `uniprotsptrembl`), attributes `ensembl_gene_id` and the filter's own, `completionStamp="1"` ("[success]") | the archive of the line's release: 115 = `sep2025.archive.ensembl.org`, archives 105 to 116 answer (104 and older are retired) | none documented; 500 values per query is the advice of biomaRt's maintainer (Bioconductor support, not Ensembl's documents); 1 per second | none found (Ensembl's BioMart pages state no limit) | all five | b, c, d, e |
| g:Profiler g:Convert | `POST api/convert/convert/` (`target = ENSG`, `numeric_ns = ENTREZGENE_ACC`), 1,000 terms per call | current instance `e114_eg62_p19` (`api/util/data_versions`); archived instances e94 to e113 (`api/util/list_archives`) at their own releases; no instance of 115 or 116 | none documented; 1 per second | none found on its API, documentation and archive pages; its API is published for programmatic use | all five | b, c, d, e |
| MyGene.info | `POST v3/query` (1,000 terms per call; "we only take the input query terms up to 5000 maximum"), `species` by taxid, scopes `ensembl.gene` (b), `symbol` (c), `symbol,alias` (d), `entrezgene,retired`, `refseq`, `uniprot` (e) | its build (`v3/metadata`: build 20261117, Ensembl 116) | "limiting the number of massive queries you make to one per second" (mygene.info/blog/update-on-mygene-down-time); 1 per second | mygene.info/terms: "TSRI authorizes you to access and use MyGene.info and the data and information"; nothing on automated use | all five | b, c, d, e |
| HGNC REST | `GET fetch/<field>/<term>`: a symbol is read as HGNC's multi-symbol checker reads it (`symbol`, then `prev_symbol`, then `alias_symbol`, until one finds a record); `entrez_id`, `refseq_accession`, `uniprot_ids` | current (`info`: lastModified 2026-11-23); its Ensembl ids are judged at 116 | "Please only send REST requests at a rate of 10 requests per second" (genenames.org/help/rest); 5 per second | none beyond the rate (the same page: an IP that harms the service is blocked) | human | c, d, e |
| UniProt ID mapping | one job per organism (`POST idmapping/run`, `from = UniProtKB_AC-ID`, `to = Ensembl`, at most 100,000 ids), status polled every 5 s, results streamed | UniProt 2026_04 (`x-uniprot-release`); its Ensembl cross-references are judged at 116 | none documented (rest.uniprot.org/help/id_mapping: 100,000 ids per job); 1 per second | CC BY 4.0 (rest.uniprot.org/help/license); programmatic access documented (help/programmatic_access) | all five (UniProt accessions only: C3) | e |
| BridgeDb web service | `POST <organism>/xrefsBatch/<source>?dataSource=En` (100 ids per call; sources `En`, `L`, `Q`, `S`, `H` for human symbols) | Ensembl 111 (`<organism>/properties`: the Ensembl data source's `DATASOURCEVERSION 111`, build 2024-10-22) | none found; 1 per second | none found | human, mouse, pig, rhesus (no cynomolgus in `contents`) | b, c (human), e |
| pyensembl 2.38.0 | `Genome(..., gtf_path_or_url=<Ensembl's GTF of the release>)`: `gene_by_id` (b), `gene_ids_of_gene_name` (c); the GTF fetched from the address pyensembl builds (`make_gtf_url`), the FASTA files not needed | any release of its species table (pig 90 to 116 on Sscrofa11.1 ...; cynomolgus 103 to 116): 115 exactly, and each line's own release | local (the GTF download: one connection) | – | all five | b, c |
| org.\*.eg.db 3.22.0 (AnnotationDbi) | `mapIds(keys, keytype = SYMBOL / ALIAS / ENTREZID / REFSEQ / UNIPROT, column = "ENSEMBL", multiVals = "list")` | the package's snapshot (Bioconductor 3.22: NCBI's gene2ensembl of its build); judged at 115: the package names no Ensembl release | local | – | human, mouse, rhesus (org.Ss.eg.db has no `ENSEMBL` column; org.Mmu.eg.db no `ALIAS`; no cynomolgus package) | c, d, e |
| AnnotationHub EnsDb 113 (ensembldb 2.34.0) | the EnsDb of Ensembl 113 (`query(AnnotationHub(), c("EnsDb", <species>, "113"))`), its genes table read as `GeneIdFilter` (b), `GeneNameFilter` (c), `EntrezFilter` (e) select it | 113, the latest EnsDb published (AHEnsDbs' metadata stops at 113); judged at 113 by the rule of equal truth. Only 113 is asked: an EnsDb of every older release would add 180 to 540 MB per species and release for lines no other tool poses at those releases | local (the hub's downloads: one connection) | – | all five | b, c, e (NCBI gene ids) |
| Seurat 5.5.1 `UpdateSymbolList` | `UpdateSymbolList(symbol)` one symbol per call, a pause of 0.4 s after each (Seurat asks HGNC REST two requests per symbol without a pause of its own: so about 5 per second, half HGNC's rate), the updated symbol to Ensembl by org.Hs.eg.db; asked after the HGNC pass, never beside it (one host) | current HGNC; judged at 116 | 5 per second | as HGNC | human | c, d |
| gget 0.30.8 `info` | `gget.info(ids, ncbi = False, uniprot = False, pdb = False, json = True)`, 200 ids per call (batches of 1,000 were answered 500 and 504) at Ensembl REST's interval, after the REST pass (one host): its Ensembl REST lookup; gget removes an id's version itself, so its two readings are one | Ensembl's current release (116); judged at 116 | as Ensembl REST | as Ensembl REST | all five | b |
| limma 3.68.5 | `alias2SymbolTable(ids, species)` (its vectorised call: each alias or symbol to one official symbol, the lowest Entrez id where several), the symbols then to Ensembl by the org package | as the org package; judged at 115 | local | – | human, mouse | c, d |
| HGNChelper 0.8.16 | `checkGeneSymbols(ids, species)` with its bundled map (dated 2026-11-01), the suggested symbols then to Ensembl by the org package | the map's date | local | – | human, mouse | c, d |
| NCBI's files (`gene_info`, `gene_history`, `gene2ensembl`, `gene2refseq`) | resolved as NCBI's README describes: a discontinued GeneID to its current one, a symbol or synonym to its GeneIDs, a GeneID or RefSeq accession to Ensembl genes by `gene2ensembl` | the files of the truth store (their date in `SOURCES.tsv`); judged at 116: NCBI's files name no Ensembl release | local | – | all five | c, d, e |

Not run, each with its reason (stated in the tables as `not run here` or `cannot be run from this machine`):

| tool | outcome | why |
| --- | --- | --- |
| Ensembl ID History Converter | cannot be run from this machine | the web tool sits behind a JavaScript challenge of the site's firewall (`/Tools/IDMapper` answers HTTP 202 with `x-amzn-waf-action: challenge` and no body); its Perl route (`IDmapper.pl`, the Perl API) needs Ensembl's MySQL servers, whose ports are closed here; it answers at 116 only |
| pybiomart 0.2.0, biomaRt (R) | scored once, as BioMart | both send the same martservice XML query BioMart is asked here: pybiomart builds it and GETs `martservice?query=<XML>` (`dataset.py`, `self.get(query=ElementTree.tostring(root))`), biomaRt POSTs it (`.submitQueryXML`, `req_body_form(query = query)`), adding a local cache; they add nothing to the answer |
| gget `search`, `ref` | cannot be run from this machine (`search`); not an id mapper (`ref`) | `search` connects to Ensembl's MySQL servers (`mysql-eg-publicsql.ebi.ac.uk` ports 3306, 5306, 4157, 3337, 5316: no connection here, probed 2026-11-24) and searches free text; `ref` lists download addresses |
| BED (Neo4j) | not run here | needs R and a Neo4j server (a 2.2 GB dump of 2026-06-27 for human, mouse, rat, pig and zebrafish, no macaques); the public instance's terms and limits are not documented; the plan's second cut |
| bioDBnet `db2db` | not run here | its web-service documentation answers HTTP 403, so its batch size, rate and terms cannot be read |
| gene-id-resolver | not run here | no tool of this name was found (PyPI, GitHub, the literature) |

## The tasks and the lines

The tasks (`task` of a G line, from its source stratum): **a** an old Ensembl id to its gene at 115 (A1 from an earlier release, A* from its origins, D1's ids from their source release); **b** an Ensembl id at its own release (A1 and A* from 115, D1's ids of 115); **c** a symbol or name (B1, B2, B7, D1's symbols); **d** a previous symbol, alias or synonym (B4, B5); **e** an external identifier (C1 NCBI gene ids, C2 RefSeq, C3 UniProt). D2 (source inference) is claimed by no tool: `not_supported` for all, stated, not asked.

### The comparison sets

Per organism and family, seeded samples of the families' published lines (the latest line of every qid of the base set and its revisions), frozen as G's published sets `g_<family>_v1` with their manifests (`families/g/sets.py`): equal allocation over cells, the first 10,000 lines asked of the web services (`tier = web`), the first 100,000 (the same and more; the whole family where it is smaller) of the local tools (`tier = local`).

| G set | source | cells | lines (web / local): pig, mouse, human, rhesus, cynomolgus |
| --- | --- | --- | --- |
| `g_a1_v1` | A1 | assembly era of the source release x distance to 115 in releases (0, 1-5, 6-15, 16-30, 31+) x expected genes (0, 1, 2+) | 10,000 / 100,000 each |
| `g_astar_v1` | A* `v1` | every hard id from 5 origins (first, last, quartiles of its releases) to 115; the web lines are 2,000 ids' | 10,000 / 44,239; 10,000 / 59,916; 10,000 / 70,804; 10,000 / 22,780; 10,000 / 11,325 |
| `g_b_v1` | B1 (with B2's strata), B4, B5, B7 | the strata | 10,000 / 100,000 each |
| `g_c_v1` | C1, C2, C3 | the strata | 10,000 / 100,000 each |
| `g_d1_v1` | D1 | the strata | pig 10,000 / 58,735; mouse, human 10,000 / 100,000; no D1 for the macaques |

A G line is the source line with `family = g`, `stratum = g_<family>_<source stratum>`, `source_family`, `source_set`, `source_stratum`, `source_cmp`, `tier`, `task`, and its genes compared without versions (`cmp = base`; `expect`, `broad`, `open_within`, `apart` and `support`'s keys without versions): most tools name genes without versions, so every subject, idtrack too, is judged at the gene. Nothing of the truth is computed by G: the truth is the families'.

### Two readings of the input

Every tool is asked each line twice where the two differ, and both tables are published (`families/g/tools.py`, `value`, `in_form`): **best**, the primary table, the reading a careful user gives a tool that takes no versions: an identifier's version removed (an Ensembl stable id, a RefSeq or a UniProt accession followed by `.<n>`; tasks a, b, e; names as they stand), since none of the tools asked documents versioned input and the pilot showed them drop it (Ensembl REST's `archive/id` leaves a versioned id out of its answer; g:Profiler and MyGene.info answer nothing); **stated**, the line exactly as it stands. Lines without a version are asked once and belong to both tables. idtrack takes versioned input and is asked as stated in both. `forms.tsv` gives per tool the posed lines whose readings differ and how many change their outcome.

### Which lines a tool is asked

A line is posed to a tool when the organism and the task are the tool's and the line's release is one it answers at (`families/g/tools.py`, `posable`): a line asked to 115 is asked of the tool's instance at 115 where it has one (BioMart's `sep2025`, pyensembl) and else of its current one (then judged by the rule below); a line asked at another release R (B1's names at R, C's identifiers at R) only of an instance at R (BioMart's archives 105 to 116, g:Profiler's archives, pyensembl). A tool's input limits are reasons too (a value with a comma in a comma-separated batch; a RefSeq prefix BioMart has no filter for; a dataset lacking a filter at that release). Every other line is `not_supported`, with its reason, outside the rate.

## Like for like

The graphs are of release 115 (no graph of 116 is built). A tool that answers at another release T (116: Ensembl REST, MyGene.info, HGNC, UniProt, NCBI's files; 114: g:Profiler's current instance; 111: BridgeDb; 113: EnsDb) is compared, on the lines asked to 115, only on those whose truth is the same at T as at 115; the rule is the same for a later and an earlier release (`families/g/likeforlike.py`). Read from the truth store (Ensembl's GTF and core tables of every release, the core cross-references of 87 on) with the harness's `make_queries.Truth`, nothing of idtrack:

1. every gene the line's truth names (`expect`, `broad`, `open_within`) is in the GTFs of 115 and of T, and no mapping session between them maps it to another id or retires it;
2. the question resolves to the same genes at T as at 115: an Ensembl id by Ensembl's history walked from the line's release to T and to 115 (`walk_series`, strict and wide); a symbol or name by the genes whose GTF name or a core cross-reference's label it is; a previous symbol, alias or synonym also by the genes it is a synonym of (`core_synonyms`, which the truth store must hold at both releases); an external identifier by the genes a core cross-reference of the accession names.

Each table gives, per tool, task and organism, the lines kept and the lines left out with the reason. The synonyms of 116 (`core_synonyms/116.tsv.gz` of each organism: 109 rows for cynomolgus to 73,746 for human, 1.4 KB to 521 KB) were fetched for this rule by family B's fetcher (`families.b.fetch_b synonyms --releases 116`, Ensembl's core dumps over HTTPS, recorded in `SOURCES.tsv`): truth for the tools answering at 116; idtrack stays at 115. The like-for-like rows of idtrack are its answers on the same kept lines; idtrack's rate on every line is reported beside them.

## The store

Everything G fetched lies in the data root, in `<families>/g/raw/`, and **nothing in it is ever removed, overwritten or truncated**: fetching is slow and the services change, so the store is the comparison's evidence and every later scoring, table and manuscript figure is made from it.

| where | what |
| --- | --- |
| `raw/<tool>/<version>/<organism>/` | each response body as it came, named by its sha256 (failed, refused and partial ones too, under their own names); `index.jsonl` (every request ever asked, one row per attempt, appended), `index.parquet` (compiled from it: the last attempt of each request, `records.load_raw`), `headers.jsonl`, `passes.jsonl` (each pass: date, release, lines, requests, the host's interval) |
| `raw/<tool>/_versions/<organism>/` | the version calls of every pass (release, build, data versions, archive lists, BioMart's filters and attributes) |
| `raw/pyensembl/<release>/<organism>/` | Ensembl's GTFs as downloaded (URL, size, sha256, md5, Last-Modified in the index and headers) and pyensembl's calls; `raw/pyensembl/_cache/` its index (rebuilt from the GTF) |
| `raw/_r_work/` | the R tools' inputs and outputs of every pass (their outputs are also stored as responses) |
| `raw/_pilot/`, `raw/_probes/` | the pilot's evaluations and the single probes of 2026-11-24 (g:Profiler's archive list and pages, UniProt's help, the gget and pybiomart packages and biomaRt's sources the design rests on) |
| `raw/STORE.tsv`, `raw/STORE_SUMMARY.tsv` | the store's manifest (`python -m idtrack_reproducibility.families.g.store`): every file with its bytes and sha256, the requests that returned it, their endpoints and dates; per tool, version and organism the files, bytes, requests, usable requests, lines held, dates and the question sets (G's published sets with their sha256); `--check` verifies every file against the index |
| `queries/<organism>/` | the question sets `g_<family>_v1` (published, fixed); `drafts/` the drawn sets of the design and pilot |
| `likeforlike/<organism>/` | the rule of equal truth per release |
| `runs/<subject>/<form>/<organism>/` | the answers made from the store and their scoring (rebuildable from the store) |
| `report/` | the tables |

**The rules of the store.** A response is written once, under its sha256, through a temporary name of its own, never over an existing file; the index and the headers are only appended to; a download streams into a part of its own, which keeps its name when the download is cut. A stored response is usable when no error is recorded and it is not a batch the service refused as a whole (`client.usable`: a status of 400 or more to a request of several lines); an unusable response stays stored, and the next pass asks its lines again under a new request. A pass asks only what the store does not hold (resume by request id and by the lines a usable response holds), so running a pass again fetches nothing that is stored.

**Scoring again without fetching**: `python -m idtrack_reproducibility.families.g.score --subjects all --organisms all` and `python -m idtrack_reproducibility.families.g.report` read the store, the question sets, the truth store and the rule of equal truth only; they open no connection (`tests/test_family_g.py` refuses every socket while the adapters read; the whole scoring was run once with the network blocked, below).

**A service that has not answered every line**: a line the tool can pose and no usable stored response answers is `not_answered` in `tool_task.tsv` (with the dates the service was asked, `raw/STORE_SUMMARY.tsv`), never scored and never dropped; when a later pass of the same tool fetches it (`families.g.fetch`, which asks only what the store lacks), `score` and `report` run again from the store join it to the tables, and nothing else is fetched again.

**Adding a question or a tool**: a new question set is published as `g_<family>_v<n>` (or a new family's set) by `families/g/sets.py`; a pass of each tool then asks only its lines (the stored ones are read, not asked again), and the store grows. A new tool is a class in `families/g/tools.py` (its organisms, tasks, releases, version calls and requests; a local tool in `families/g/local.py`), its host with its documented rate in `client.HOSTS`, and its adapter in `families/g/adapters.py`; it is fetched once into `raw/<tool>/` and scored with the others.

## The client

`families/g/client.py`: one session per host group with a pool of one connection, the interval between requests from the host's rate (above; BioMart's archive hosts are one group), a 429 or 503 waited out for its `Retry-After` (else 60 s doubling to 15 minutes, at most 6 h in all), other 5xx, time-outs and connection errors waited 30 s, 2, 8 and 30 minutes and then recorded with their error and asked again by the next pass; an archive that does not answer its version calls is left for the next pass. One pass of a host is one job, so a host is asked through one connection over all jobs.

**The raw store** (docs/records.md, "Another tool"): `<families>/g/raw/<tool>/<version>/<organism>/`: every response body as it came, in a file named by its sha256; `index.jsonl` (one row per request, `records.RAW_COLUMNS`, appended after its file is written) and `index.parquet` compiled from it (`records.load_raw`); `headers.jsonl` (the response headers: rate-limit and release headers); `passes.jsonl` (each pass: when, the release, the lines, the requests, the host's interval). The version calls of every pass (REST `info/data` and `info/software`; BioMart's `registry`, `version`, `filters`, `attributes` of each archive; g:Profiler's `data_versions` and `list_archives`; MyGene's `metadata`; HGNC's `info`; UniProt's `configure/idmapping/fields` with its release header; BridgeDb's `properties` and `contents`) are stored in `<families>/g/raw/<tool>/_versions/<organism>/`. A request's id is the sha256 of its method, URL, parameters and body: a stored answer is not asked again, and a line held by a stored response of the same endpoint is not asked again, so a pass resumes. pyensembl's GTFs are stored the same way (URL, size, sha256, md5, Last-Modified), its index in `<families>/g/raw/pyensembl/_cache` (rebuilt from the GTF); its calls and the R tools' calls are stored as they returned (`calls_<pass>.jsonl.gz`, indexed like a response). Raw responses and downloaded files are kept forever.

## From responses to answers

One adapter per tool reads its stored responses only (`families/g/adapters.py`) and writes the runner's answers (`<families>/g/runs/<tool>/<version>/<organism>/answers/`, `runner.answer_schema()`): the genes the tool named for the line without versions, the tool's own status (`current`, `replaced`, `retired`, `found`, `not_found`, `other_annotation_only`), `error` where its response failed, the line's share of its request's seconds; a line the tool cannot pose is written with outcome `not_supported` and its reason, which the scorer keeps (`not_supported`, outside the rate). Ensembl REST: a current id is its own answer, a retired one its `possible_replacement` ids. BioMart: the genes of the rows whose filter value is the line's (without case, as the filter matches); a response without `[success]` is an error. g:Profiler: the `converted` ENSG ids of the line's `incoming`. MyGene.info and UniProt: the genes of every hit, without the ids of Ensembl's other annotations of the species (pig's breed assemblies: ids absent from Ensembl's GTF of the tool's release, read from that GTF). HGNC: the `ensembl_gene_id` of the first field of the chain that found a record. BridgeDb: the `En:` ids. The scorer is the common one (`families/scorer.py`), on G's sets, with `run.json` naming the subject (`kind = tool`, `name`, `version`, `asked_at`, `endpoint`). idtrack's rows on the same lines are the families' answers of the tree measured (`families/<family>/runs/<tree>/<organism>/answers/`, read only), scored on G's lines; lines no family run holds are asked through the adapter on the tree and its graphs.

## The pilot

At most 100 web lines per tool and organism (the first 20 of each G set), asked 2026-11-24. Seconds per request (mean, from the raw index) and the like-for-like share of the posed lines:

Per tool and organism: lines posed of the pilot's lines (pig, mouse, human 100; macaques 80, no D1) / of them left out by the rule of equal truth / right : wrong : undetermined under `all` on the lines kept (the common scorer; a pilot, not a result: 20 lines per G set). The rest of the pilot's lines are `not_supported` for the tool (mostly task a, history, which only Ensembl REST claims besides idtrack).

| tool (release) | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | --- | --- | --- |
| Ensembl REST, `archive/id` and `xrefs/symbol` (116) | 94 / 6 / 57:19:12 | 95 / 7 / 65:15:8 | 92 / 4 / 60:19:9 | 73 / 6 / 44:19:4 | 71 / 6 / 46:18:1 |
| Ensembl REST, `lookup/symbol` (116) | 26 / 5 / 13:2:6 | 23 / 7 / 10:3:3 | 21 / 4 / 6:7:4 | 17 / 6 / 10:0:1 | 14 / 3 / 6:4:1 |
| BioMart (115; archives at R) | 51 / 0 / 30:13:7 | 50 / 0 / 36:8:3 | 43 / 0 / 26:11:4 | 36 / 0 / 24:8:1 | 39 / 0 / 23:12:2 |
| g:Profiler (114; archives at R) | 48 / 1 / 27:13:7 | 47 / 2 / 26:16:3 | 41 / 0 / 16:20:5 | 33 / 0 / 18:14:1 | 37 / 1 / 17:17:2 |
| MyGene.info (116) | 48 / 6 / 22:13:7 | 47 / 7 / 27:10:3 | 41 / 4 / 23:10:4 | 33 / 6 / 17:9:1 | 37 / 6 / 14:16:1 |
| HGNC (116) | – | – | 39 / 4 / 16:15:4 | – | – |
| UniProt (116) | 4 / 0 / 2:2:0 | 8 / 0 / 5:3:0 | 9 / 0 / 7:2:0 | 8 / 0 / 5:3:0 | 3 / 0 / 0:3:0 |
| BridgeDb (111) | 22 / 1 / 10:10:1 | 24 / 0 / 17:7:0 | 37 / 3 / 21:10:3 | 16 / 0 / 9:7:0 | – |

BioMart's lines at archives 105 to 114 were still being asked when this was written (1 to 3 per organism). The local tools: pyensembl answered 93 lines at 115 (each GTF fetched in about 1 s, indexed in about a minute); org.\*.eg.db, limma and HGNChelper answered the human, mouse and rhesus pilot lines (33 to 42 per organism) in 5 to 23 s per tool and organism.

Seconds per request in the pilot (the raw index): Ensembl REST `archive/id` 1.55 s per 47 ids, `lookup/symbol` 5.8 s per 20 symbols (so its batches are 200), `xrefs/symbol` 0.39 s; BioMart at 115 0.09 s per query of 5 values; g:Profiler 0.34 s per 12 terms; MyGene.info 0.19 s per 7 terms; HGNC 0.05 s; UniProt 0.03 s per call (a job finishes within 15 s); BridgeDb 0.03 s per call. Every response was HTTP 200 but BioMart's archives that answered 500 (below). Ensembl REST's `archive/id` leaves out an id given with a version (7 of 53 ids of the mouse pilot): the tool's behaviour, measured as it is.

## Expected wall time

The web lines of the full pass (per organism 50,000 lines; the macaques 40,000) and what each tool can pose of them (`families/g/tools.py`), with the pilot's seconds per request and each host's interval; one job per host:

| host | lines posed (pig, mouse, human, rhesus, cynomolgus) | requests | wall time |
| --- | --- | --- | --- |
| Ensembl REST | 46,155; 46,326; 46,344; 36,118; 36,087 | about 95,000 (`xrefs/symbol` one per line) | about 15.5 h |
| HGNC | human 20,341 | about 28,000 | about 1.6 h |
| g:Profiler | 24,202; 25,217; 24,023; 19,965; 21,507 | about 210 | under 1 h |
| MyGene.info | 22,025; 23,098; 21,876; 17,655; 18,339 | about 120 | under 1 h |
| BioMart | 23,702; 24,480; 23,146; 19,199; 19,906 | about 690 | under 1 h where the archives answer; an archive that answers 500 costs its retries (40 minutes per request) and is then left for the next pass |
| UniProt | 3,840; 4,231; 4,317; 3,642; 3,815 | 5 jobs | minutes |
| BridgeDb | 9,540; 10,195; 18,835; 9,808; – | about 490 | minutes |
| pyensembl (local tier) | the b and c lines of 100,000 per family at 115 and at their own releases | one GTF per (release, organism) | hours (downloads and indexing; tens of GTFs per organism) |

## Results

None yet.

## Known limitations

- A tool that answers at another release than the graphs is compared on the lines of equal truth only; the lines left out are counted, not judged.
- The services change on their own rhythm: a pass is one date per service (no second pass).
- NCBI's and the org packages' Ensembl ids come from NCBI's comparison with an Ensembl release the files do not name; they are judged at the release stated above.

## Reproducing from scratch

```
# the environments (Python with pyensembl; R with Bioconductor's packages; gget)
mamba env create -p <prefix>/idtrack_g_py -f environments/g_py.yml
mamba env create -p <prefix>/idtrack_g_r -f environments/g_r.yml
mamba env create -p <prefix>/idtrack_g_gget -f environments/g_gget.yml
# the synonyms of 116 for the rule of equal truth (family B's fetcher), the sets (published once), the rule
python -m idtrack_reproducibility.families.b.fetch_b synonyms --releases 116
python -m idtrack_reproducibility.families.g.sets --organism <o>
python -m idtrack_reproducibility.families.g.likeforlike --organism <o> --releases 116,114,113,111
# one pass per host (one job per host; a pass resumes by request id), the local tools
python -m idtrack_reproducibility.families.g.fetch --tool <ensembl_rest|biomart|gprofiler|mygene|hgnc|uniprot|bridgedb> --organisms all
<idtrack_g_gget>/bin/python -m idtrack_reproducibility.families.g.local --tool gget --organisms all --tier web
python -m idtrack_reproducibility.families.g.local --tool <pyensembl|ncbi_files> --organisms all
python -m idtrack_reproducibility.families.g.local --tool <orgdb|limma|hgnchelper|ensdb|seurat> --organisms all --rscript <idtrack_g_r>/bin/Rscript --workdir <folder>
# the answers from the stored responses, the scoring, the tables
python -m idtrack_reproducibility.families.g.score --subjects all --organisms all --tree main
python -m idtrack_reproducibility.families.g.report
```
