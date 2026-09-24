# MATLAB → Python migration notes

This document tracks the port of pleioFDR from MATLAB/Octave to a Python ≥ 3.13
package. It lists where each MATLAB function goes, the MATLAB behaviours that
have to be reproduced exactly, and every intentional difference.

## Proposed package layout

```
pyproject.toml                 requires-python = ">=3.13", console script `pleiofdr`
src/pleiofdr/
    __init__.py
    __main__.py                python -m pleiofdr
    cli.py                     runme.m: argument parsing, config → Options, load inputs, run
    config.py                  TextConfig.m + runme.m defaults; safe MATLAB-literal parser
    options.py                 pleioOpt.m as a dataclass (t1breaks/t2breaks/qqbreaks as properties)
    io.py                      reference/trait .mat loading (v5 via scipy.io, v7.3 via h5py),
                               load_gwas.m, refinfo (.ref) loading, randprune_file, result.mat
    ldmatrix.py                LDMatrix: CSC pattern (int64 indptr, int32 indices), column access
    pruning.py                 FastPrune.m, random_prune_idx_amd.m, random_prune_idx_amd_fb.m
    stats.py                   binofit (Clopper–Pearson), GCcorrect_logpvec.m, fisher_comStats.m,
                               check_sample_overlap.m, correct_sample_overlap.m, logp↔z helpers
    smoothing.py               SparseSmooth2d.m
    lookup.py                  lookup_table.m (ind_look), conj_lookup_table.m, cond_FDR_amd.m,
                               MATLAB-compatible interp2 / histc / hist3
    fdr.py                     pleioFDR_amd.m
    loci.py                    ind_loci_idx.m, locusnumber.m
    analysis.py                pleiotropy_analysis.m (orchestration, logging to <outputdir>/*.log)
    output.py                  save_to_csv.m, save_fdr.m, save_zscore.m
    plotting/
        __init__.py            save_figure.m, filter_points_for_plotting.m, backend selection
        qq.py                  plot_qq_amd.m (QQ and TDR)
        enrichment.py          plot_enrichment_amd.m
        lookup.py              plot_lookup.m
        manhattan.py           plot_Manhattan.m
    fuma/
        combine.py             fuma/*.R ported to pandas (cond / conj-lead / conj-snps)
        novelty.py             fuma/conj_fuma_combined_novelty.py
tests/
    fixtures/make_golden.m     generates MATLAB reference outputs on the demo data
    fixtures/golden_*.mat      (generated) pipeline intermediates per configuration
    fixtures/unit_matlab.mat   (generated) per-function reference outputs
    test_*.py                  unit tests + golden parity tests + end-to-end CLI test
ref4pleioFDR/toolkit/*.py      kept as standalone scripts, modernised for 3.13
run_batch.py                   calls `pleiofdr --config` instead of matlab
MultipleRunUtility.sh          calls `pleiofdr --config` instead of matlab
legacy_matlab/                 original .m files, moved here once parity is proven
```

## Dependencies

| Package    | Why |
|------------|-----|
| numpy      | arrays |
| scipy      | `scipy.io` (v5 .mat), `scipy.sparse` + `spsolve` (smoothing), `scipy.stats`/`scipy.special` (norm, chi2, F, gammaincc) |
| pandas     | `.ref` file, CSV tables, fuma scripts |
| matplotlib | figures (Agg backend unless `onscreen=true`) |
| h5py       | v7.3 .mat files; the full `ref9545380_1kgPhase3eur_LDr2p1.mat` is v7.3 |
| numba      | FastPrune and locusnumber are greedy, inherently sequential walks over up to ~9.5M SNPs × `randprune_n` iterations. They cannot be vectorised, and a pure-Python loop would be orders of magnitude slower than MATLAB's JIT. numba 0.67 supports Python 3.14. |

Optional extras: `ref` (`intervaltree`, for `ref4pleioFDR/toolkit/ld_informed_annot*.py`).
Dev tools: `pytest`, `ruff`.

## Function map

| MATLAB | Python | Notes |
|---|---|---|
| runme.m | `cli.main` | same config keys and defaults |
| TextConfig.m | `config.TextConfig` | `get_cell` uses a literal parser, not `eval` |
| pleioOpt.m | `options.Options` | property validation kept (stattype, exclude regions) |
| pleiotropy_analysis.m | `analysis.run` | |
| load_gwas.m | `io.load_gwas` | first variable starting with `logp` / `z` (case-insensitive), in file order |
| pleioFDR_amd.m | `fdr.pleio_fdr` | |
| cond_FDR_amd.m | `lookup.cond_fdr` | |
| lookup_table.m | `lookup.lookup_table` | |
| conj_lookup_table.m | `lookup.conj_lookup_table` | only reached from plot_lookup when lookups are missing |
| SparseSmooth2d.m | `smoothing.sparse_smooth_2d` | |
| binofit_wrap.m / binofit (toolbox) | `stats.binofit` | Clopper–Pearson via the F distribution, like MATLAB |
| binofit_dale.m | — | Octave-only fallback; dropped |
| GCcorrect_logpvec.m | `stats.gc_correct_logp` | |
| fisher_comStats.m | `stats.fisher_combined` | |
| check/correct_sample_overlap.m | `stats.check_sample_overlap` / `stats.correct_sample_overlap` | |
| FastPrune.m | `pruning.fast_prune` | numba kernel |
| random_prune_idx_amd(_fb).m | `pruning.random_prune_idx` | `repeats` = default / maxout / none |
| ind_loci_idx.m | `loci.independent_loci` | |
| locusnumber.m | `loci.locus_number` | numba kernel |
| save_to_csv / save_fdr / save_zscore.m | `output.*` | same file names, headers, number formats |
| save_figure.m | `plotting.save_figure` | |
| filter_points_for_plotting.m | `plotting.filter_points` | vectorised with `np.unique` on bin ids |
| plot_qq_amd.m, plot_enrichment_amd.m, plot_lookup.m, plot_Manhattan.m | `plotting.*` | |
| is_octave.m | — | dropped |
| suplabel.m, plot_qq_annot.m | — | not called anywhere in the pipeline; not ported |

## MATLAB behaviours that must be reproduced

- **Indexing.** Exclusion regions (`exclude_from_fit`) and locus bounds are 1-based inclusive `[from, to]` in MATLAB. Python stores 0-based half-open ranges internally, and `result.mat` keeps MATLAB's conventions where it exposes indices.
- **Column-major order.** `SparseSmooth2d` vectorises the image with `im(:)`, and `binofit` results are reshaped column-major. Python uses `order="F"`.
- **NaN in min/max.** MATLAB `min`/`max` ignore NaN. In `cond_FDR_amd`, `min(len, max(1, 1 + lp/Δ))` maps NaN to index 1, so `fdrvec0` gets `looktable(1,1)` for undefined SNPs (`fdrvec` is NaN-masked afterwards; `fdrvec0` is not). Use `np.fmin`/`np.fmax`.
- **Sort order.** `FastPrune` uses `sort(abs(x), 'descend')`, which is stable and puts NaN first (NaN and Inf are then dropped). Use a stable argsort on `-abs(x)` over the finite entries.
- **`histc(x, edges)`.** Bin k holds `edges(k) <= x < edges(k+1)`. The last bin holds `x == edges(end)`. Values outside the edges and NaN are ignored.
- **`hist3(..., 'Edges', ...)`.** Same edge rules per dimension. The lp1 edges are shifted down by half a bin, and a trailing `Inf` edge catches `lp2 > 3`. The last row and column are trimmed afterwards.
- **`interp2(V, x, y)`.** Bilinear on the grid `1:ncols × 1:nrows`, NaN outside. The indices are clamped into range before the call. NaN propagation from neighbouring NaN cells is checked against the golden fixtures.
- **`prctile`.** MATLAB's definition uses midpoint interpolation (`method="hazen"` in numpy) and ignores NaN. It is used by the in-house GC.
- **`median`.** Ignores NaN only where the MATLAB code filters explicitly. `median(sig0)` over a vector with NaN returns NaN in MATLAB, and numpy does the same.
- **`norminv` of 0 or 1.** Gives ∓Inf, which matches `scipy.stats.norm.ppf`.
- **`unique(a)`.** Returns the first occurrence after sorting. It is used by `randprune_repeats='none'`.
- **`num2str` in `TextConfig.declare`.** Defaults are stored as strings, and `get_bool` treats any non-zero number as true. `true`/`false` literals come out as 1/0.
- **`%g` / `%e` / `%f` formatting.** File names use `%g` of `fdrthresh` (for example `0.01`, `5e-08`). CSV rows use C-style `%e` / `%f`, and `%d` of NaN prints `NaN`.

## Quirks preserved for parity (candidate bugs; not fixed)

1. `fisher_comStats.m` applies the GC correction only when `~exist('ivec0','var')`, which is never the case. So the Fisher statistic is **never** GC-corrected, and the `pruneidx_fc` argument has no effect. The port keeps this behaviour.
2. `pleiotropy_analysis.m` leaves `pruneidx_fc` undefined when `randprune=false` and `fishercomb=true`, so MATLAB errors. The port uses no pruning in that case instead of crashing (**intentional difference**).
3. In `plot_qq_amd.m` / `plot_enrichment_amd.m`, `qqci_max`/`qqci_min` are built from the already-updated `qqmat`. They only affect the dashed CI lines, which are drawn only when `show_ci=true`. Kept.
4. `mat_qq*` and `mat_enrich*` in `result.mat` hold the matrices for the **last** conditioning trait only, because each trait overwrites them. Kept.
5. `SparseSmooth2d` uses `lamvec(1)` for both directions. Kept; both entries are `smf` anyway.
6. `plot_Manhattan.m` has `showgenes = false` hard-coded, so gene labels are never drawn. The gene-label code is not ported.

## Intentional differences

- **Figure formats.** MATLAB writes `.fig`, `.png` and `.svg`. `.fig` is a MATLAB-only format and is not produced. PNG and SVG keep the same base names.
- **Random numbers.** Prune indices are drawn with `numpy.random.Generator` (seedable with `--seed`). They differ from MATLAB's `unifrnd` stream, so exact comparisons with MATLAB go through `randprune_file` (prune indices exported from MATLAB), which is already a config option.
- **Octave branches** (`is_octave`, `binofit_dale`, gnuplot toolkit) are removed.
- **`exit_matlab_upon_completion`** is accepted and ignored. `mlibrary` is accepted and ignored.
- **`randprune=false`.** Works instead of erroring (see quirk 2).

## Verification status

- MATLAB R2022b is installed on the development machine but cannot run: it is an x86_64 build, and Rosetta 2 is not installed. Octave is not installed. The golden fixtures (`tests/fixtures/make_golden.m`) still need to be generated.
