# Designing the specificity criterion

## Core principle: the per-cell-type peak is the only screening criterion

Whether a gene is "specific to a certain cell type" is judged solely by whether its **peak cell
type** in the HPA single-cell-type matrix falls inside the target cell-type group. This is the
conclusion reached after two rounds of iteration in a real project:

- The first version also required "highest expression in bone marrow tissue" for megakaryocytic
  specificity, which wrongly killed a set of genes that are highly specific at the cell level but
  not at the tissue level (a tissue is a mixture of cell populations, so cell-type specificity is
  diluted at the tissue level).
- The correct approach: **the cell-type peak decides inclusion; tissue information is additional
  annotation only**. Most genes that satisfy cell-type specificity are indeed also highest in the
  corresponding tissue (such as bone marrow), but that is a result, not a condition.

## Decision logic (implemented in screen_specificity.py)

For each gene:

1. `top_cell` = the cell type with the highest expression among all cell types; the peak must be > 0.
2. `driver` = the first group in `specificity_groups` that contains `top_cell` (the order of the
   groups in the config is the priority, e.g. put the narrowest lineage group first).
3. **Selected**: `driver` is non-empty and the highest value inside the group ≥ the highest value
   outside the group (i.e. the peak lies inside the target group).
4. **Borderline**: the highest value inside the group ≥ the highest value outside the group ×
   `borderline_fraction` (default 0.85) but the selection criterion is not met; listed separately
   for manual review — typical examples are genes such as RUNX1 whose peak sits at the edge of the
   target group.
5. `specificity_ratio` = highest inside the group / highest outside the group, used for ranking
   and candidate grading; housekeeping genes (such as ACTB, TMSB4X) are expressed at very high
   levels but have a ratio close to 1, so when commenting on candidates judge mainly by ratio and
   secondarily by absolute level.

## Additional tissue-level metrics

With `reference_tissue` configured (e.g. "bone marrow"), the output includes that tissue's
expression value and its rank among all tissues; used for display and cross-validation only, never
for screening.

## Cell-type names

Cell-type names must match the column names returned by HPA exactly (case, hyphens, plurals),
e.g. `Megakaryocyte-Erythroid progenitors`, `Platelets`, `cDC`, `pDCs`.
After running step 1 you can check the names in the config against the header of the output Excel;
the effect of a misspelled name is that the type silently takes part in no group, so it is
advisable to check the `specificity_groups` echo in screening_summary.json.
