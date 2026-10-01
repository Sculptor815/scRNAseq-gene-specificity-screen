# HPA API reference and known pitfalls

## Endpoints used

1. **Column discovery**: `GET https://www.proteinatlas.org/search?search=gene_name%3ATP53`
   extracts the column names from `data-abbr="..."` in the HTML. `sc_RNA_*` are the single-cell
   type columns (about 154, excluding `sc_RNA__tau`); `t_RNA_*` are the tissue columns (about 51,
   excluding `t_RNA__tau`).
   Column names change with the HPA version, so discover them dynamically on every run — never
   hard-code them.

2. **Bulk data**: `GET https://www.proteinatlas.org/api/search_download.php`
   Parameters: `search=gene_name:SYM1,SYM2,...`, `format=json`, `columns=<comma-separated column names>`, `compress=no`.
   Returns a list of JSON records; each record contains `Gene`, `Ensembl` and the requested
   columns, whose names look like `Single Cell Type RNA - Platelets [nCPM]` /
   `Tissue RNA - bone marrow [nTPM]`.

3. **Single-gene detail** (functional annotation): `GET https://www.proteinatlas.org/{ENSEMBL}.json`
   Contains `Gene description`, `Protein class` and so on; suitable for a per-gene annotation page.

## Known pitfalls (all encountered in practice)

- **Batch size**: 40 genes per batch, 4 concurrent requests, 120 s timeout, 4 retries (backoff
  3/6/9 s) is the verified, safe configuration.
  The bulk endpoint occasionally does not return for a long time; putting **write failures into
  the retry condition as well** hides real errors — network retries should only cover the request
  itself.
- **Write permission**: the output directory may not be writable (e.g. some external or synced
  drives). Run in a writable working directory first, then copy everything to the target
  directory.
- **Several Ensembl entries with the same name**: for some gene symbols HPA returns multiple
  Ensembl entries (observed in this project for MATR3, POLR2J3). Merging the matrices by
  upper-casing `Gene` and taking the first record picks the wrong one; keep all entries and match
  by Ensembl downstream, and list the duplicated symbols as warnings for manual review. The
  canonical entry is usually the one with the more plausible expression level.
- **Unmatched genes**: about 0.4% of gene symbols cannot be found in HPA (aliases or deprecated
  symbols). Keep the whole row empty and list them under `missing_genes` in the metadata; do not
  silently drop them, and do not backfill from other sources.
- **Missing values**: a missing column in a record = not measured / not included by HPA; write
  None or blank, never 0.
- **Caching**: each batch result is written to disk as atomic JSON keyed by a hash of
  `(dataset, symbols, columns)`; re-running after an interruption resumes automatically, and
  `--refresh` forces a re-download.
