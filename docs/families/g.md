# G: comparison with other tools

The tools people use for the questions idtrack answers, asked the same lines and judged by the same truth and the same scorer as idtrack ([methods](../methods.md#g-comparison-with-other-tools)): their versions and the dates they were asked, their raw responses kept as they came, a question a tool cannot pose scored `not_supported`, and tables of tool by task.

## State

The comparison sets are published (`g_<family>_v1`, 2026-11-24, with their revision sets `g_b_r1`, `g_b_r2`, `g_c_r1`) and the rule of equal truth is computed for releases 116, 114, 113 and 111; every service was asked in one pass (2026-11-24 and 25; BioMart's archives answered only in part, below), every subject scored from the store with the network off, and the disagreements with idtrack read (below). The pilot below was asked before the sets were published, at most 100 lines per tool and organism; its responses are kept with the others (endpoints without a form); a pass that would ask a request the pilot asked reads the pilot's response instead (the same request, the same day).

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

The lines are those of the judge users get (the families' `scores_now`): a decided judge's revision set (`<set>_dec`, decisions not built yet) is not G's source. `g_b_v1` holds B7 lines drawn from `b7_dec`; the revision set `g_b_r1` of each organism (pig 202, mouse 255, human 303, rhesus 201, cynomolgus 172 lines; `sets.py --revise-judge`) holds those qids' lines of the judge users get, and is their truth. A source family's revision set published after G's draw is joined the same way (`sets.py --revise-sources`: the lines whose source line changed its truth, as the source's current line): `g_b_r2` (pig 3,253, mouse 1,807, human 5,496, rhesus 2,829, cynomolgus 1,119 lines: B's revision sets of 2026-11-24 and 25) and `g_c_r1` (mouse 4,754, human 7,794: C's UniProt lines' revisions), so that G judges by the families' current truth.

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

One adapter per tool reads its stored responses only (`families/g/adapters.py`) and writes the runner's answers (`<families>/g/runs/<tool>/<version>/<organism>/answers/`, `runner.answer_schema()`): the genes the tool named for the line without versions, the tool's own status (`current`, `replaced`, `retired`, `found`, `not_found`, `other_annotation_only`), `error` where its response failed, the line's share of its request's seconds; a line the tool cannot pose is never sent to the judge: G writes its row of the record with outcome `not_supported` and the reason (`score.not_supported_rows`, `outcomes/part_not_supported_of_1.parquet`), outside the rate; the common scorer is used unchanged. A request is asked once (its id is the hash of its method, URL, parameters and body): a request of one line (Ensembl REST's `xrefs/symbol`, HGNC's `fetch`) answers every line that asks the same string, and the pilot's response every line of the pass that asks the pilot's request (`adapters.stored_answer`). A line its source family withdrew (`excluded` in a revision set) is scored by no subject, as the families do, and counted apart (`excluded` in `tool_task.tsv`). Ensembl REST: a current id is its own answer, a retired one its `possible_replacement` ids. BioMart: the genes of the rows whose filter value is the line's (without case, as the filter matches); a response without `[success]` is an error. g:Profiler: the `converted` ENSG ids of the line's `incoming`. MyGene.info and UniProt: the genes of every hit (for an Ensembl id looked up in MyGene.info, the id itself where a hit lists it: a hit is an NCBI gene that may list several Ensembl genes), without the ids of Ensembl's other annotations of the species (pig's breed assemblies: ids absent from Ensembl's GTF of the tool's release, read from that GTF). HGNC: the `ensembl_gene_id` of the first field of the chain that found a record. BridgeDb: the `En:` ids. The scorer is the common one (`families/scorer.py`), on G's sets, with `run.json` naming the subject (`kind = tool`, `name`, `version`, `asked_at`, `endpoint`). idtrack's rows on the same lines are the families' answers of the tree measured (`families/<family>/runs/<tree>/<organism>/answers/`, read only), scored on G's lines; lines no family run holds are asked through the adapter on the tree and its graphs.

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

BioMart's pilot lines at archives 105 to 114 (1 to 3 per organism) are not in the table. The local tools: pyensembl answered 93 lines at 115 (each GTF fetched in about 1 s, indexed in about a minute); org.\*.eg.db, limma and HGNChelper answered the human, mouse and rhesus pilot lines (33 to 42 per organism) in 5 to 23 s per tool and organism.

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

The passes of 2026-11-24 and 2026-11-25 (one date per service; the dates, versions and endpoints of every response in `raw/STORE_SUMMARY.tsv`), scored from the store with the network off. The tables are made in `<families>/g/report/`; the summary tables are published in this repository, in [`results/g/main/all/`](../../results/g/main/all/) with their provenance (each names its table in the data root, every registered scored folder's `run.json` and the store's summary: `report --publish`), and the per-line tables (`disagreements.tsv`, about a million lines, and `reading.tsv`, every sampled line read) stay in the data root, listed in the family's data manifest (`manifests/families/g.tsv.gz`), as the other families keep their per-line tables in their runs. The tables: `tool_task.tsv` (per subject, reading of the input, organism and task: the tier's lines, the lines their source family withdrew (`excluded`), `not_supported`, posed, `not_answered`, left out by the rule of equal truth, kept, right, wrong, undetermined, error under `all` and `best`, the rate with its Wilson interval, coverage, and idtrack's rate on the same kept lines), `not_supported.tsv` (the reasons), `forms.tsv` (per tool the posed lines whose two readings differ and how many change their outcome), `by_distance.tsv`, `speed.tsv`, `failures.tsv`, `disagreements.tsv` with `disagreement_classes.tsv` (every disagreement with idtrack classed by the two verdicts of the common scorer: idtrack wrong, the tool wrong, both wrong, the truth cannot decide; a seeded sample of at most 200 per subject, task and organism stratified by the outcome pair is marked for reading) and the reading's `reading.tsv`, `reading_classes.tsv`, `reading_summary.tsv` (below).

The rate of each tool and idtrack's rate on the same kept lines, in per cent, with the kept lines (form `best`, view `all`; a dash where the tool cannot pose the task for the organism):

<!-- table g:rates -->
| tool | task | pig | mouse | human | rhesus | cynomolgus |
| --- | --- | --- | --- | --- | --- | --- |
| BioMart | b | 100 / 100 (1,173) | 100 / 100 (1,721) | 100 / 100 (1,478) | 100 / 100 (1,506) | 100 / 100 (2,228) |
| BioMart | c | 91 / 99 (10,124) | 87 / 98 (9,876) | 68 / 98 (8,952) | 87 / 100 (5,412) | 88 / 100 (6,677) |
| BioMart | d | 73 / 100 (2,951) | 90 / 100 (3,383) | 91 / 99 (3,041) | 76 / 100 (2,478) | 35 / 99 (1,350) |
| BioMart | e | 70 / 100 (8,908) | 70 / 100 (8,879) | 66 / 100 (8,045) | 70 / 100 (8,444) | 72 / 100 (8,476) |
| BridgeDb | b | 93 / 100 (1,172) | 100 / 100 (1,499) | 91 / 100 (1,285) | 99 / 100 (1,506) | – |
| BridgeDb | c | – | – | 70 / 98 (8,607) | – | – |
| BridgeDb | e | 71 / 100 (7,796) | 67 / 100 (7,936) | 68 / 100 (7,626) | 68 / 100 (8,294) | – |
| EnsDb 113 | b | 100 / 100 (8,595) | 100 / 100 (13,216) | 100 / 100 (13,921) | 100 / 100 (9,863) | 100 / 100 (13,124) |
| EnsDb 113 | c | 89 / 100 (80,068) | 74 / 98 (84,766) | 62 / 97 (86,952) | 82 / 100 (68,852) | 94 / 100 (73,354) |
| EnsDb 113 | e | 81 / 100 (11,215) | 60 / 99 (13,580) | 58 / 99 (13,342) | 77 / 100 (11,994) | 76 / 100 (10,442) |
| Ensembl REST | a | 76 / 98 (24,171) | 92 / 99 (23,259) | 94 / 98 (24,503) | 82 / 100 (18,494) | 100 / 100 (17,772) |
| Ensembl REST | b | 100 / 100 (1,173) | 100 / 100 (1,721) | 100 / 100 (1,478) | 100 / 100 (1,506) | 100 / 100 (2,228) |
| Ensembl REST | c | 93 / 98 (9,783) | 90 / 98 (9,677) | 64 / 98 (8,847) | 92 / 100 (5,054) | 93 / 100 (6,320) |
| Ensembl REST | d | 79 / 100 (2,374) | 90 / 100 (3,116) | 85 / 99 (3,011) | 82 / 100 (2,438) | 45 / 99 (1,323) |
| Ensembl REST | e | 68 / 100 (8,200) | 69 / 100 (8,432) | 60 / 100 (7,779) | 69 / 100 (8,298) | 69 / 100 (7,389) |
| Ensembl REST `lookup/symbol` | c | 90 / 98 (9,783) | 83 / 98 (9,677) | 67 / 98 (8,847) | 86 / 100 (5,054) | 81 / 100 (6,320) |
| Ensembl REST `lookup/symbol` | d | 10 / 100 (2,374) | 15 / 100 (3,116) | 6 / 99 (3,011) | 12 / 100 (2,438) | 32 / 99 (1,323) |
| gget | b | 100 / 100 (1,173) | 100 / 100 (1,721) | 100 / 100 (1,478) | 100 / 100 (1,506) | 100 / 100 (2,228) |
| g:Profiler | b | 100 / 100 (1,173) | 100 / 100 (1,721) | 90 / 100 (1,478) | 100 / 100 (1,506) | 100 / 100 (2,228) |
| g:Profiler | c | 92 / 99 (10,274) | 87 / 98 (9,848) | 69 / 98 (9,365) | 88 / 100 (5,659) | 92 / 100 (7,772) |
| g:Profiler | d | 17 / 100 (2,964) | 15 / 100 (3,654) | 8 / 100 (3,530) | 16 / 100 (3,192) | 27 / 99 (1,618) |
| g:Profiler | e | 73 / 100 (9,482) | 73 / 100 (9,642) | 69 / 100 (8,941) | 73 / 100 (9,588) | 74 / 100 (9,084) |
| HGNC | c | – | – | 71 / 98 (8,847) | – | – |
| HGNC | d | – | – | 92 / 99 (3,011) | – | – |
| HGNC | e | – | – | 40 / 100 (7,779) | – | – |
| HGNChelper | c | – | 72 / 99 (88,435) | 61 / 97 (87,412) | – | – |
| HGNChelper | d | – | 32 / 100 (35,349) | 82 / 100 (33,844) | – | – |
| limma | c | – | 64 / 99 (88,435) | 62 / 97 (87,412) | – | – |
| limma | d | – | 80 / 100 (35,349) | 67 / 100 (33,844) | – | – |
| MyGene.info | b | 98 / 100 (1,173) | 100 / 100 (1,721) | 100 / 100 (1,478) | 99 / 100 (1,506) | 100 / 100 (2,228) |
| MyGene.info | c | 90 / 98 (9,811) | 86 / 98 (9,715) | 69 / 98 (8,891) | 83 / 100 (5,081) | 80 / 100 (6,342) |
| MyGene.info | d | 92 / 100 (2,385) | 80 / 100 (3,130) | 73 / 99 (3,024) | 90 / 100 (2,442) | 72 / 99 (1,325) |
| MyGene.info | e | 51 / 100 (8,200) | 69 / 100 (8,432) | 63 / 100 (7,779) | 53 / 100 (8,298) | 55 / 100 (7,389) |
| NCBI's files | c | 23 / 100 (78,091) | 60 / 99 (88,164) | 63 / 97 (87,141) | 74 / 100 (64,161) | 67 / 100 (71,826) |
| NCBI's files | d | 10 / 99 (10,130) | 82 / 100 (35,287) | 71 / 100 (33,715) | 78 / 99 (6,588) | 49 / 99 (3,434) |
| NCBI's files | e | 40 / 100 (43,657) | 67 / 100 (40,655) | 63 / 100 (32,748) | 64 / 100 (44,196) | 57 / 100 (37,098) |
| org.\*.eg.db | c | – | 60 / 99 (88,435) | 59 / 97 (87,412) | 74 / 100 (68,983) | – |
| org.\*.eg.db | d | – | 80 / 100 (35,349) | 68 / 100 (33,844) | – | – |
| org.\*.eg.db | e | – | 68 / 100 (83,657) | 57 / 100 (77,681) | 73 / 100 (82,129) | – |
| pyensembl | b | 100 / 100 (8,595) | 100 / 100 (13,223) | 100 / 100 (13,934) | 100 / 100 (9,863) | 100 / 100 (13,124) |
| pyensembl | c | 90 / 100 (92,841) | 78 / 98 (102,278) | 68 / 98 (103,226) | 86 / 100 (81,506) | 95 / 100 (83,390) |
| Seurat | c | – | – | 59 / 98 (8,891) | – | – |
| Seurat | d | – | – | 38 / 99 (3,025) | – | – |
| UniProt | e | 65 / 100 (3,831) | 48 / 100 (4,226) | 74 / 100 (4,314) | 78 / 100 (3,641) | 75 / 100 (3,812) |
<!-- /table -->

Only Ensembl REST (`archive/id`) and idtrack pose task a (an old Ensembl id to its gene at 115), 44 to 50 % of the web lines. On symbols (task c), aliases (d) and external identifiers (e) no tool's rate reaches idtrack's on the same kept lines in any organism; on Ensembl ids of 115 (task b) most tools and idtrack are right on every line. The rule of equal truth leaves out at most 4.3 % of a tool's posed lines (BridgeDb at 111; NCBI's files 2.7 %, every other tool 2 % or less). What a tool's wrong lines are is read below: most are questions the tool does not claim to answer, and the rates are not a ranking of the tools' quality on the questions they claim.

**BioMart** was asked while some of its archives answered HTTP 500 for hours (2026-11-24 and 25: `dec2021` (105) and `oct2022` to `may2024` (108 to 112)): its lines at 115 and at the archives that answered are scored; the lines of the archives that did not answer are `not_answered` (not answered by the service on these dates, never scored as wrong): pig 34 (c), 263 (d), 247 (e); mouse 40, 275, 306; human 61, 507, 457; rhesus 84, 717, 558; cynomolgus 363, 351, 458. A later pass of the same fetch asks only these and `score --subjects biomart`, `report` and `reading` join them to the tables without anything else being fetched again.

## The reading of the disagreements

Every disagreement with idtrack is listed (`disagreements.tsv`), and a seeded sample of at most 200 per subject, task and organism, stratified by the pair of outcomes, is read from the store, the truth store and the families' runs with the network off (`families/g/reading.py`; `reading.tsv` one row per sampled line and side, `reading_classes.tsv` the classes with three examples each, `reading_summary.tsv` the groups). A sampled line carries the weight of its sampling cell (the cell's lines over its sampled lines), so that the shares below estimate the shares among all of a subject's disagreements of that task, every organism together.

A line where **the tool is wrong** is classed by what the tool did, read at the release it answers at: the service could not read the input (an error answer); the tool named only genes of another Ensembl annotation of the species (NCBI's and MyGene.info's pig genes on breed assemblies); an id that is no gene of the tool's release (a retired id, another annotation's); a gene on a patch or haplotype; a part of the truth's genes (right under the view `best`); the authority's gene (HGNC, MGI, VGNC, NCBI, UniProt name it for the input; the truth reads Ensembl's tables); a reading of another release (the line asks a name, id or accession as an earlier release read it, and the tool reads today's tables); another gene the string names at the tool's release (a shared name, a synonym, a cross-reference); a successor Ensembl's history allows (the wide walk); a synonym a tool that reads display names only does not read (Ensembl REST's `lookup/symbol`); no answer where Ensembl's tables of the tool's release name the truth's gene by the input; no answer where they do not; truly wrong (another current gene by none of these readings). The classes fall in four groups: **truly wrong** (another current gene, or no answer where Ensembl's tables name the truth's gene), **right under another reading** (another annotation, another release's id, a patch copy, a part of the truth, the authority's gene, another gene the string names, a successor of the wide walk), **not a question the tool claims** (a reading of another release, a synonym a display-name lookup does not read, no answer where Ensembl's tables do not name the gene either), **the service failed**.

A line where **idtrack is wrong** is classed by its cause in the families' runs and that cause's row of the known-limitations table (`results/limitations/limitations.tsv`; a cause's row in any family, the source family's map of its wrong answers' classes), else by what idtrack did against Ensembl's tables of the release the line asks to: a string that stands for a missing value, refused (`L-missing-values`); a symbol in another case that is another gene's synonym, or a string no Ensembl table of any release holds, answered (the graph holds Ensembl's strings only: a string held in another case; `L-case-synonym`); else a new cause (`new: ...`), each a row of the maintainers' backlog with its reproducer.

Per tool and task: the sampled lines where the tool is wrong (every organism together) and the weighted shares of the four groups, in per cent; the sampled lines where idtrack is wrong and the weighted share of them that a known limitation states:

<!-- table g:reading -->
| tool | task | tool wrong, read | truly wrong | right under another reading | not a question the tool claims | the service failed | idtrack wrong, read | of them a known limitation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| BioMart | c | 311 | 16 | 18 | 66 | 0 | 110 | 83 |
| BioMart | d | 761 | 21 | 6 | 74 | 0 | 28 | 100 |
| BioMart | e | 578 | 0 | 15 | 85 | 0 | 19 | 100 |
| BridgeDb | b | 215 | 100 | 0 | 0 | 0 | 0 | – |
| BridgeDb | c | 45 | 4 | 22 | 73 | 0 | 32 | 72 |
| BridgeDb | e | 463 | 3 | 13 | 84 | 0 | 19 | 100 |
| EnsDb 113 | c | 196 | 29 | 19 | 53 | 0 | 114 | 80 |
| EnsDb 113 | e | 430 | 0 | 0 | 100 | 0 | 132 | 97 |
| Ensembl REST | a | 207 | 1 | 99 | 0 | 0 | 96 | 100 |
| Ensembl REST | c | 433 | 0 | 42 | 53 | 5 | 76 | 87 |
| Ensembl REST | d | 746 | 0 | 14 | 86 | 0 | 29 | 97 |
| Ensembl REST | e | 529 | 2 | 17 | 82 | 0 | 20 | 95 |
| Ensembl REST `lookup/symbol` | c | 313 | 13 | 43 | 44 | 0 | 111 | 85 |
| Ensembl REST `lookup/symbol` | d | 853 | 5 | 0 | 95 | 0 | 29 | 100 |
| g:Profiler | b | 142 | 100 | 0 | 0 | 0 | 0 | – |
| g:Profiler | c | 363 | 17 | 23 | 60 | 0 | 94 | 85 |
| g:Profiler | d | 857 | 84 | 0 | 16 | 0 | 29 | 100 |
| g:Profiler | e | 573 | 1 | 14 | 85 | 0 | 20 | 95 |
| HGNC | c | 46 | 0 | 28 | 72 | 0 | 38 | 76 |
| HGNC | d | 177 | 3 | 38 | 59 | 0 | 0 | – |
| HGNC | e | 148 | 32 | 1 | 67 | 0 | 6 | 100 |
| HGNChelper | c | 68 | 44 | 18 | 38 | 0 | 66 | 94 |
| HGNChelper | d | 102 | 81 | 7 | 12 | 0 | 26 | 81 |
| limma | c | 68 | 51 | 7 | 41 | 0 | 66 | 98 |
| limma | d | 98 | 70 | 20 | 10 | 0 | 46 | 19 |
| MyGene.info | b | 42 | 100 | 0 | 0 | 0 | 0 | – |
| MyGene.info | c | 330 | 41 | 25 | 33 | 0 | 99 | 87 |
| MyGene.info | d | 810 | 58 | 24 | 19 | 0 | 12 | 33 |
| MyGene.info | e | 513 | 15 | 14 | 71 | 0 | 15 | 100 |
| NCBI's files | c | 226 | 44 | 37 | 19 | 0 | 121 | 75 |
| NCBI's files | d | 421 | 50 | 41 | 9 | 0 | 96 | 94 |
| NCBI's files | e | 354 | 22 | 32 | 46 | 0 | 79 | 99 |
| org.\*.eg.db | c | 112 | 60 | 12 | 29 | 0 | 67 | 85 |
| org.\*.eg.db | d | 89 | 89 | 5 | 7 | 0 | 46 | 18 |
| org.\*.eg.db | e | 170 | 9 | 11 | 80 | 0 | 102 | 98 |
| pyensembl | c | 197 | 34 | 20 | 46 | 0 | 128 | 60 |
| Seurat | c | 37 | 32 | 24 | 43 | 0 | 36 | 83 |
| Seurat | d | 134 | 93 | 3 | 4 | 0 | 0 | – |
| UniProt | e | 771 | 7 | 3 | 90 | 0 | 2 | 100 |
<!-- /table -->

The rules that keep the comparison fair, each with the reading's measured reason, applied from the stored responses alone (nothing fetched again):

- **The lines G judges by** are the families' current lines: a line drawn from a decided judge's revision set (`*_dec`) is replaced by the line of the judge users get (`g_b_r1`), and a line whose source family revised its truth after G's sets were drawn by the source's current line (`g_b_r2`: B's revisions; `g_c_r1`: C's revision of its UniProt lines; the sampled lines where idtrack was wrong on C's UniProt versions were all such lines).
- **MyGene.info's answer to an Ensembl id** (task b) is the id itself where a hit lists it: a hit is an NCBI gene that may list several Ensembl genes.
- **A request is asked once** and answers every line that asks it: Ensembl REST's `xrefs/symbol` and HGNC's `fetch` are one request per string, and a stored request names only the line it was first asked for (pig's symbols: 1,561 of 11,440 web lines ask a string asked for another line, or by the pilot).
- **Lines their source family withdrew** (`excluded`: 6,217, of them 6,193 human RefSeq lines of C2) are scored by no subject and counted apart.
- **A line read at 115** and kept by the rule of equal truth reads alike at the tool's release: it is never "a reading of another release".

What each tool is not built for (the group "not a question the tool claims", never a fault of the tool):

- **Ensembl REST `archive/id`** (task a) names an id's possible replacements at the release that retired it, one step of Ensembl's history; it does not follow them on to the current release. Its sampled wrong lines are almost all such steps: replacements Ensembl retired later, ids that are no gene of 116 ("right under another reading"). A user who asks it again for each replacement walks the history by hand.
- **Ensembl REST `xrefs/symbol`, BioMart, g:Profiler, MyGene.info, HGNC, UniProt, BridgeDb, EnsDb, pyensembl, org.\*.eg.db** read a name, synonym or accession as their one release (or today) reads it. A line that asks a string as an earlier release read it (B1, B4 and B5 names read at their source release, NCBI gene ids and UniProt accessions retired since, RefSeq versions no longer current) asks for history these tools do not claim: such lines are the column "not a question the tool claims" above (from 7 % of org.\*.eg.db's sampled wrong aliases to all of EnsDb's identifiers). Where a tool does claim a kind of history it is asked it and judged on it: NCBI's `gene_history` for discontinued GeneIDs (NCBI's files; MyGene.info's `retired` scope), HGNC's previous symbols and aliases, UniProt's secondary accessions.
- **Ensembl REST `lookup/symbol`** reads display names only, never synonyms: a synonym it does not read is 77 % of its wrong lines on task d (weighted), which is why its rate there is low.
- **Ensembl REST `xrefs/symbol`** cannot read a string with a Greek letter (B7: its database answers HTTP 400 with a collation error), "the service failed".
- **NCBI's files and MyGene.info for pig** link NCBI's genes to Ensembl's annotations of other pig breeds (`ENSSSCG002...`), not to the reference assembly's genes: such a line has no answer of the reference ("another annotation's genes only").
- **g:Profiler, BridgeDb and MyGene.info on Ensembl ids** (task b) know only the genes of their own databases: an Ensembl gene of the release without an entry there has no answer (their sampled wrong lines of task b are all such misses, "truly wrong").

idtrack's wrong lines (the last two columns): most are a known limitation (bridges where Ensembl's history ends, `L-bridges`; Ensembl's history carrying a gene to another place, `L-history-vs-place`; names from other sources, `L-other-name-sources`; a symbol in another case that is another gene's synonym, or a string held only in another case, `L-case-synonym`; values of another database, `L-other-database-values`; a string that stands for a missing value, refused, `L-missing-values`). Every class of them is a row of the known-limitations table (`results/limitations/limitations.tsv`, family `g`: its size in G's sample per organism; `families/g/limits_g.py`). The others had no cause in the families' runs (their lines carry no explained path there); each was asked again with `convert(..., explain=True)` on the released package and the graphs measured, and is a row of that table by G's map (`families/g/limitations_map_g.tsv`: the first row whose strata, organisms and reading's class match the line), and of the maintainers' backlog:

| G's row | lines | what the explained path shows | whose |
| --- | --- | --- | --- |
| `L-alias-copies-in-answer` | 51 | human NCBI aliases of B5: idtrack answers the primary assembly's gene and reports the copies on patches or haplotypes apart (rule `copies_apart`), while B5's line lists them in its answer | the truth's line |
| `L-old-name-bridge` | 74 | a name read at an old release (B1, B2) whose genes Ensembl's history carries to nothing at 115 (the truth: nothing) is answered through a bridge (rule `external_bridge`): 30 at the old gene's place, 35 across an assembly change Ensembl publishes no mapping for (pig Sscrofa10.2 to 11.1 24, rhesus 10, mouse 1; pig `AKR1C3` at 72 also reaches the gene of AKR1C1, another NCBI symbol of the old gene), 8 through a copy outside the reference sequence; one mouse line (`RF00001` at 92) is the RFAM accession of other genes, no bridge | idtrack's bridges; B1 tests no place for a name of several genes, so the 30 at the old place, undetermined by the locus rule for a name of one gene, count wrong |
| `T-copy-display-name` | 17 | a name only copies on patches display at its release (B1 `patchonly`: mouse 15 at 82 to 92, human `CACNA1C-AS2` at 114, `C4B_2` at 108), carried by the gene of the reference sequence at the copies' place in another name source (Vega's clone name, miRBase's label, NCBI's symbol, UniProt's gene name): idtrack answers that gene with the copies apart (rule `copies_apart`), the truth expects the copy | the truth's |
| `L-copy-name-other-gene` | 1 | human `HLA-DRB3` at 70: displayed only by genes of the MHC haplotypes, carried as UniProt's gene name by HLA-DQB1 and HLA-DRB1 on chromosome 6, which idtrack answers (the haplotypes' genes apart) | idtrack's |
| `T-uniprot-genename` | 15 | the best answer is the truth's gene; `all` adds a gene that carries the name as UniProt's gene name, which Ensembl filed as `Uniprot_genename` at human 62 to 72 and pig 63 to 72 and B's truth does not read (`truth_b.NAME_DB` reads `Uniprot_gn`), so it is in neither `expect` nor `broad` (human 13, B1 and B4; pig 2) | the truth's |
| `L-mirna-host-name` | 8 | rhesus and pig at 76: Ensembl displays a miRNA's name for the miRNA gene and its host gene; the graph holds miRBase's names on miRNA genes only (rule `mirbase_genes`), so the miRNA gene alone is answered | idtrack's rule |
| `L-bare-numbers` | 6 | human: spreadsheet serials of dates (B7, 4) and NCBI gene ids NCBI discontinued (C1, 2) read as Havana numbers of releases 51 and 52 (rule `bare_numbers`) | idtrack's rule |
| `L-retired-ncbi-synonyms` | 3 | mouse NCBI gene ids NCBI discontinued or replaced (C1), which Ensembl's records keep as synonyms of a gene's NCBI or UniProt link: the gene of Ensembl's record | Ensembl's records against NCBI's |
| `T-read-release` | 2 | the truth reads `CCDS3` (B5, `held_in_core`) at 87 and `RNU6-283P` (C2) at 89, while at 115 another database holds each (CCDS's NOC2L, HGNC's one gene), where rule `start_release` reads it | the truth's |
| `L-history-vs-place` | 1 | human `TCP10` of D1 at 75: Ensembl's history keeps ENSG00000203690, while ENSG00000290751 at the same place holds more of the old model | a known limitation (D's text) |

## Known limitations

- A tool that answers at another release than the graphs is compared on the lines of equal truth only; the lines left out are counted, not judged.
- The services change on their own rhythm: a pass is one date per service (no second pass).
- NCBI's and the org packages' Ensembl ids come from NCBI's comparison with an Ensembl release the files do not name; they are judged at the release stated above.
- BioMart's archives of 105 and 108 to 112 did not answer on the dates asked: their lines are `not_answered`, counted per organism and task, until a later pass fetches them.
- The reading classes a seeded sample (at most 200 lines per subject, task and organism) by rules over Ensembl's tables of the tool's release; its shares are estimates weighted back to every disagreement, and a class is what the tables show, not a statement of the tool's intent.
- The truth judges every subject by Ensembl's tables of 115 (or of the tool's release by the rule of equal truth): a tool built on another source (NCBI, HGNC, UniProt) is right by that source where the classes say so ("the authority's gene", "another annotation's genes only") and still counted wrong in the rate.

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
# the revision sets: lines drawn from a decided judge's set, and lines whose source family revised its truth later
python -m idtrack_reproducibility.families.g.sets --organism <o> --revise-judge
python -m idtrack_reproducibility.families.g.sets --organism <o> --revise-sources
# the answers from the stored responses, the scoring, the tables, the reading, the store's manifest (no network)
python -m idtrack_reproducibility.families.g.score --subjects all --organisms all --tree main
python -m idtrack_reproducibility.families.g.report
python -m idtrack_reproducibility.families.g.reading
python -m idtrack_reproducibility.families.g.store
# the published tables and this document's tables written from them; G's rows of the known-limitations table;
# the results' provenance
python -m idtrack_reproducibility.families.g.report --publish
idtrack-repro known-limitations --only g
idtrack-repro results --only g,limitations
```
