# stormbreaker 0.3.47

* English guides, runnable tutorials and main workflow help examples have been updated. Core CI now includes tutorial behavior assertions.
* `sb_complete_run()` retains log-rank status and errors in survival risk-group tables and exported issues; usable predictions and plots survive a failed test.
* `sb_model_shap()` and `sb_tree_shap()` validate a saved fit once, reusing its prediction state inside permutations while retaining output checks. The estimator, seeds and score scales are unchanged.
* `surv_risk_plots()` distinguishes successful, single-group, no-event and failed log-rank tests. Computation errors warn and retain their original message instead of becoming an unexplained NA.

# stormbreaker 0.3.46

* Every public diagnostic and plotting entrance now has a concise recommended name. New aliases include `forest_diag()`, `path_diag()`, `pamr_diag()`, `mclust_diag()`, `native_diag()` and `plot_diag()`.
* The 85 additional names also cover remaining plot/data helpers and further shorten SHAP and generic drawing calls: `plot_shap_imp()`, `plot_shap_dep()`, `plot_shap_one()`, `plot_cm()` and `plot_surv()`. All 221 names from 0.3.45 remain supported.
* Shipped current guides and examples use the preferred short names. The complete map and a native Lasso example ship with the package. Model implementations, arguments, result fields, coordinates and export behavior are unchanged.

# stormbreaker 0.3.45

* Six executable scenarios ship with five synthetic CSVs and an offline illustrated vignette. Each scenario documents input/output objects and plotting.
* `models()`, `fit_final()`, `predict_final()`, `finish_run()`, `clinical_plots()` and common plotting functions have shorter names. The 45 short names are exact aliases, not new model implementations; all 176 previous exports remain supported without deprecation warnings.
* `plot_calib()`, `plot_perf()`, `plot_prob()` and `plot_heatmap()` replace long plotting names in the short-name guide. Arguments, modeling defaults, numeric results and concise figure styling are unchanged.

# stormbreaker 0.3.44

* Figures use concise default styling: explanatory subtitles and captions remain available as metadata or report notes rather than drawing inside the plot. Font sizes and numeric results are unchanged.
* Clinical display labels are shorter; ROC AUC annotations use curve colours instead of repeating model names. Direct performance figures also use readable model names while preserving original IDs.
* SHAP and coefficient axes use readable score units. Existing plotting calls and modeling defaults are unchanged.

# stormbreaker 0.3.43

* Packages the revised 0.3.41 source with doctor-facing plotting repairs and algorithm-view guidance. No production R code is changed during this rebuild.
* The full model menu and original 0.3.41 modeling defaults are retained.

# stormbreaker 0.3.41

* `bi_run()` and `surv_run()` accept optional decision thresholds; `surv_run()` also accepts a prediction horizon. Omitted settings retain the previous behavior.
* Doctor-facing comparison figures use readable model names instead of internal parameter IDs; saved IDs and numeric results are retained.
* `sb_clinician_views()` accepts doctor run and completed-run results directly. Named-cohort calibration and decision curves now retain their saved coordinates alongside ROC views, without refitting.
* `sb_plot_nomogram()` keeps instructions in plot metadata and report captions instead of a figure footer; point conversion and probabilities are unchanged.
* Unifies the formal-run 0.3.40 implementation with the later binary modeling entrypoints. Flexible-survival numerical settings, canonical/legacy name handling and retired-learner rules are retained; historical research results are not regenerated.
* Binary brglm2, LDA, HDA, FDA, MDA and local LDA use the formal-run training-only constant-sum representation repair when applicable. New states retain original feature names and the fitted transform; historical project representation states and ordinary native states remain readable. Linear brglm2/LDA views map back to original gene units.
* `bi_pipeline()` accepts ordered binary method names; `bi_run()` enumerates supported routes within one chosen pool. The binary whole-menu default is stages 1 and 2; three-stage expansion is explicit. `bi_selector()` is internal compatibility code rather than an exported entrypoint.
* `sb_start()` presents modeling and validation without a separate binary feature-selection task. Advanced feature analysis, survival syntax and multiclass interfaces are retained.

# stormbreaker 0.3.40

* The Gompertz time-unit regression test now supplies the fitted training-coordinate transform to its native-engine comparison. Production R code, prediction states and model defaults are unchanged from 0.3.39.

# stormbreaker 0.3.39

* Flexible parametric survival fitting now uses training-only centered/unit-SD numerical coordinates and a BFGS iteration ceiling of 1000. Spline fitting starts from a converged training-fold intercept-only model. Convergence and Hessian checks remain enforced.
* Public flexible-survival coefficients, risk scores, survival probabilities and compact states remain on the original input scale. Model and state metadata record numerical fitting settings; old compact prediction states remain readable.

# stormbreaker 0.3.38

* Removed the binary and survival REML GAM learners and binary SLOPE, including all menu combinations, because of excessive runtime and computational cost. GAM boosting remains available. Historical releases/results and the separate multiclass implementation are unchanged.

## Earlier unversioned workflow notes

- `bi_run()` and `surv_run()` print a compact summary while retaining all
  result fields. Optional plotting failures are recorded without losing the
  model comparison or fitted models.
- `sb_complete_run()` reuses saved algorithm plots when no new explanation
  is requested and prints a summary rather than patient-level predictions.
- `sb_export_run()` keeps column headers in empty cohort and figure indexes.

# stormbreaker 0.3.37

- Algorithm identifiers no longer end in `_fixed`. This naming change does
  not change numerical defaults, native parameter selection, or fitted values.
  The `fixed_features` selector retains its existing meaning and name.
- `sb_pipeline()` accepts former learner identifiers and constructs pipelines
  with canonical names. Registered menus and display labels use the new names.
- `sb_fit_final()` and `sb_predict_final()` normalize legacy model metadata
  without changing expression values, gene names, or saved native model payloads.
- Stack training and repeated comparisons preserve the stack structure while
  normalizing legacy member identifiers.

# stormbreaker 0.3.36

- `sb_binary_menu_pipelines()` constructs the full registered binary menu
  from one declared feature panel, with explicit Top-K and glmnet tuning
  settings. Seven MIME family references are labelled as family-level only;
  fixed adapters are not presented as reproducing MIME's native tuning grids.

# stormbreaker 0.3.35

- Survival CForest predictions split assessment sets larger than the training
  set into native prediction batches, avoiding the reproduced process crash
  without changing the forest or per-patient probability calculation.

- ahaz learners resolve tied deaths with reproducible, bounded fit-only time
  offsets. Original endpoints remain unchanged for calibration and evaluation;
  fitted audit metadata records the policy, seed, affected counts and size.

- `sb_survival_menu_pipelines()` saves the full survival registry and all 117
  MIME variants in a source-backed, deduplicated menu. MIME-common tuning
  settings use 10 inner folds (PLS-Cox: 5), seed 1, and up to 10000 GBM trees.
  Fixed-feature routes use the caller's same declared feature vector;
  extensions retain native defaults and low-level APIs remain unchanged.

# stormbreaker 0.3.34

- StepCox and CV-tuned PLS-Cox separate likelihood convergence precision from
  the native infinite-coefficient warning threshold. StepCox uses tighter
  likelihood convergence and saves its effective controls; PLS-Cox no longer
  implicitly tightens the warning threshold when tightening convergence.
  Genuine fitting warnings still stop the affected route.
- CV-tuned PLS-Cox uses tighter native Cox likelihood convergence in both
  component search and final fitting. Component limits, selection criterion
  and strict warning handling are unchanged; effective controls are saved.
- `sb_predict_final()` and comparison tables label exactly constant
  predictions without changing fitted states or metrics. `sb_select_model()`
  excludes constant validation predictions and models constant within every
  development assessment fold, recording the reason without inspecting test
  or external cohorts.
- Gamboost survival fits use native dense spline calculations to avoid sparse
  QR augmentation warnings on larger analysis sets. Model definitions and
  strict handling of other fitting warnings are unchanged.
- Gompertz fits use time-scale-aware native starting values and optimizer
  scaling, preventing nonfinite initial likelihoods from long follow-up.
  Input time units, risk scores, native model family and convergence checks
  are preserved.
- LPC reports the native automatic threshold search's 50-predictor minimum
  when smaller inputs produce an empty threshold, without padding predictors
  or substituting a threshold.

- Ranger C-index survival fits support the opt-in execution option
  `stormbreaker.ranger_cindex_threads` (1--4, default 1). Parallel repeated-CV
  workers force one engine thread. Forest predictions are unchanged in paired
  tests; native importance reduction may differ at floating-point precision.

- Akritas survival prediction computes native curves once per exactly distinct
  query CDF value within each call and restores the original patient order.
  Training, risk and horizon probability definitions remain unchanged.

- Bagged-kNN survival prediction reuses bootstrap-specific quadratic products
  to screen neighbors, then computes all potentially selected distances with
  the original centered expression. Native neighbor ordering, KM aggregation,
  training and bootstrap budgets are unchanged; near-ties are recomputed.

- CV-tuned GBM-Cox reuses warning-free search diagnostics when refitting the
  selected tree count. Serial refitting retains the final native CV fold to
  reproduce the original full-model RNG stream; an upper-bound optimum reuses
  the already identical complete search fit. Warned searches and one-tree
  fits retain their original execution path. Tree budgets and folds are unchanged.

- Removed binary and survival VSOLassoBag PST selectors and their routes.
  Existing PST pipeline specifications are no longer supported; no substitute
  selector is chosen automatically. Historical analysis records are unchanged.

- Bounded research execution can retain completed Boruta selection
  across downstream worker timeouts, without substituting fixed
  features or sharing across different analysis samples. Ordinary comparisons
  retain their existing automatic in-memory reuse.

- Repeated CV accepts optional `workers = 2` for local repeat-level parallelism,
  preserving within-repeat sharing and result order. The serial default and old
  positional arguments remain unchanged. Sequential caller-RNG learners require
  serial execution rather than silently changing their random streams.

- Lasso/Elastic Net Cox feature selection reuses the identical full-analysis
  lambda path already fitted by native CV rather than fitting it again.
  RSF minimal-depth selection omits the unused proximity matrix and uses private
  reference containers during native tree traversal to avoid recursive copying;
  native depth, importance, selection, error-rate and RNG semantics are retained.

- Expression-local learner reuse now covers the registered learner families in
  both CV and final fitting, including complete analysis-local native searches.
  Matching includes relevant tuning context and native RNG state where required;
  failed or warning-bearing fits remain separate. Nonconsecutive identical final
  inputs can reuse a fit. A shared 256 MiB retained-learner-entry budget releases
  least-recently-used entries; oversized entries are not retained. This is not
  a process-memory limit, and released fits are recomputed normally if needed.

- RSF keeps its native 0.632 SWOR sampling rule but uses a namespace-scoped
  sample-size callback, avoiding capture of the full native training frame in
  fitted and portable models. Native risk and survival predictions are unchanged.

- BART and cforest survival predictions now reuse one warning-free native
  distribution for risk and horizon probability within the same CV assessment
  or saved-state prediction call. Native values, settings and per-call failure
  handling are unchanged; randomForestSRC predictions retain their original path.

- Identical `none`/`fixed_features` final penAFT, PLS-Cox-CV and CoxBoost-CV
  fits now share their complete training-local search and model in the existing
  expression-local reuse context. Per-route diagnostics, warning metadata and
  horizon processing remain separate. Fits with stored component-search warnings
  are not reused, just like fits with emitted warnings.

- Expression-local final-fit reuse now includes identical `none`/`fixed_features`
  RLT survival and CV-tuned GBM-Cox fits. GBM shares its complete native search
  and selected model only on identical training inputs and settings; it does not
  transfer tuning between folds. Each route retains its own diagnostics, state
  and horizon extraction. Warnings and failed fits are never reused.

- Survival CV reuses the same assessment-fold bagged-kNN curve for
  risk and fixed-horizon survival probability. Horizon-only errors retain
  successful risks, warned predictions are not reused, and no curve is shared
  across folds or stored in the CV result.

- Bagged-kNN prediction reuses each bootstrap covariance inverse across query
  patients while retaining native Mahalanobis distances, ordering, KM and
  aggregation. Akritas uses the same native multivariate empirical-CDF routine
  without reconstructing R6 objects, and evaluates its risk grid and requested
  horizon together. Model budgets, saved states and prediction scales are unchanged.

- Within expression-local reuse, identical `none`/`fixed_features` final fits
  of fixed survival-SVM learners share native training. Patient IDs, ordered
  inputs, outcomes and learner settings must match. Every route retains its own
  state extraction, horizon baseline and diagnostics; CV and tuning are unchanged.
  No warned or failed training result is reused.

- Selector reuse now includes survival Ridge, PLS-Cox and survival-SVM weight
  screening and Elastic Net selection without changing their native budgets.
  Captured diagnostic histories are distinguished from compact fits, including
  Lasso and CoxBoost paths, so reuse cannot discard a requested path or loading.

- Boruta accepts explicitly declared smaller selector budgets
  through `selector_options`; defaults and default pipeline IDs are unchanged.
  Changed budgets are recorded in IDs and fitted selection metadata.
- Survival boosting selectors also
  reuse identical successful, warning-free analysis-fold fits.

- Added native Akritas conditional survival prediction with fixed empirical-CDF
  bandwidth 0.5, direct survival probabilities, and negative step-RMST risk.
  Portable prediction requires sensitive training predictors and paired outcomes;
  patient IDs and row names are removed.

- Added binary PLSR screening using native logistic PLS coefficients times
  analysis-fold SD, retaining the top half (minimum two positive weights when
  available); no threshold search or gene padding.

- Added binary and survival XGBoost selectors based on positive native tree
  gain, using the existing fixed native engines and fold-local refitting.

- Binary GamBoost now also supports fold-local active-spline gene selection,
  with native Binomial-logit loss; it is not an alias of linear glmBoost.

- Survival `GamBoost` selection now uses native Cox bbs spline terms (dfbase 4,
  100 steps, shrinkage 0.1), independently of the linear `glmBoost` selector.
  Only active original features are passed downstream; displayed scores are
  spline coefficient-vector norms. Selection is repeated within training folds.

- `bi_run()` and `surv_run()` accept a named `pipelines` list, comparing single learners and different ordered selector chains in one call with shared folds.
- `sb_complete_run()` reports external survival risk-group counts at the locked development threshold even without plots; single-group limitations appear in exported issues and `risk-groups.csv` without changing predictions or cutpoints.

- `sb_run_issues()` handles native plot-unavailability tables row by row, without interrupting doctor-facing runs or duplicating the same issue.

- New multilayer binary networks use two smaller hidden layers (3, 2), probability squared-error loss rather than cross-entropy, and up to 500000 optimization steps; the convergence threshold is unchanged. New FDA fits use linear polynomial bases. Local LDA uses a training-size-capped neighbourhood of at least twice the predictor-plus-class count (minimum 30). Multiclass FDA/local LDA share these defaults.
- New RBF survival-SVM fits use the median squared distance of standardized training predictors for their native kernel width. Cox calibration remains required and unchanged. Previously saved pipeline parameters retain their old loss, basis, neighbourhood and kernel-width rules when refitted.

- Binary `nnet` now allows 2000 optimization iterations rather than 200; convergence is still required and the architecture, penalty, seed and tolerance are unchanged. Previously saved 200-iteration pipelines retain their recorded budget instead of being silently upgraded.

- Reloaded GAM fits use the native mgcv predictor for smooth plots even in a fresh R session. Empty plot collections return an empty page list instead of failing during pagination.

- Shared plot styling uses clearer sans-serif text, lighter grids and consistent spacing. Mixture-correlation heatmaps have square cells and angled labels; selector OOB curves use a restrained blue/orange/grey palette. Plotted values are unchanged.
- `sb_plot_pages()` also paginates multi-column feature facets, keeping all observations and uncertainty/support layers rather than compressing many GAM panels onto one page.

- Univariable AUC and Cox screening expose their saved native ranking statistics through the unified plotting entrypoint without refitting.

- `bi_run()` and `surv_run()` default to automatic plots: comparison-only for batches, algorithm-specific plots for one candidate. `sb_complete_run()` continues the selected candidate. Explicit `plots = TRUE/FALSE` retains eager/disabled behaviour. `multi_run()` similarly defers batch method plots.
- `sb_compare_models()` and `sb_compare_cohorts()` no longer open a browser by default; explicit `view = TRUE` and report exports remain available.

- Windows multinomial ABESS fits run in a fresh R worker to contain a reproduced native access violation in mixed-engine sessions; data, native algorithm and fitting parameters are unchanged, and the original fitted object is returned for prediction and plotting.
- LightGBM uses native short logloss aliases at the training boundary so Bioconductor list formatting cannot truncate metric names; recorded canonical metrics and native predictions are preserved. Multiclass RDA explicitly calls klaR's predictor to avoid unrelated RDA method registrations.
- Multiclass GAM, prediction-rule ensembles, PAMR, Gaussian processes and neighbourhood methods expose their native smooths, rule/CV results, centroid contrasts, basis weights and query neighbourhoods respectively. Neural connection colours use a symmetric pseudo-log scale without modifying weight values.

- Multiclass conditional forests, bagging, RLT, oblique forests and weighted-subspace forests include native tree-size diagnostics, available split counts and individual-tree OOB errors without surrogate fitting.

- Multiclass EDDA/MclustDA include native component means and available BIC search tables; HDDA and HDA include original eigenvalue spectra and projection loadings without refitting.

- Multiclass LIBLINEAR L1/L2 logistic and XGBoost linear expose original class-specific weights; XGBoost uses full-precision saved weights rather than rounded text dumps.

- The multiclass method entrypoint now returns compact plots with per-figure explanations retained separately in `notes`, including native tree explanations.

- Multiclass conditional inference and evolutionary trees expose every original branch and training-node count through the existing tree renderer, without surrogate fitting or pruning.

- Multiclass FDA/MDA/PDA, stabilized LDA and shrinkage discriminant analysis include native projected centroids, loadings, mixture weights or shrinkage diagnostics as applicable.

- Multiclass kernel naive Bayes includes saved class-conditional density curves; RDA includes native class means and regularized within-class correlation heatmaps.

- Multiclass C5.0 trees/rules expose native feature-usage plots; regularized forests expose native impurity importance and class-specific out-of-bag error histories.

- Multiclass logistic, stepwise logistic and bias-reduced models include class-specific coefficient plots; stepwise models also expose the original accepted AIC history.

- Multiclass neural networks expose the original input-to-hidden and hidden-to-output connection weights, including biases, as separate heatmap panels.

- `sb_plot_pages()` paginates sample-level algorithm figures without dropping observations; `sb_plots()` now returns ready-to-render pages and separate explanatory notes.
- Selector coefficient paths retain gene legends; SuperPC CV views separate component counts and identify the component actually used for threshold selection.
- Multiclass decision trees include their native node structure; multiclass SVMs include the original pairwise decision-score distributions, not only support-vector counts.
- `sb_export_run()` paginates algorithm figures, respects requested figure dimensions and exports a separate figure-notes table.

- `sb_complete_run()` continues one explicitly chosen doctor-run pipeline without repeating comparison, reusing its final fit and evaluating named external cohorts.
- `sb_complete_run()` connects survival risk panels and optional cohort explanations using development backgrounds; failing cohorts remain separate from successful results.
- `sb_export_run()` writes comparison, stage genes, cohort predictions/metrics, PDF figures and readable issues to a new directory without overwriting existing output.

# stormbreaker 0.3.33

- `bi_selector("SVM_RFE")` adds linear recursive elimination with grouped inner-CV feature-count selection; `sb_plots()` supplies accuracy and error curves.

- `sb_threshold_metrics()` includes binary Cohen Kappa; `sb_plot_thresholds(metrics = "kappa")` preserves negative agreement instead of clipping it.

# stormbreaker 0.3.32

- `sb_selector_diagnostics(type = "boruta_boxplot")` now draws true predictor/shadow boxplots from saved Boruta histories, also routed through `sb_plots()`.
- `sb_pipeline()` extends Lasso, Enet, CoxBoost and RSF screens to compatible fixed-parameter survival learners while preserving native-internal-CV restrictions.
- `surv_selector()` adds Ridge, PLSR and linear SVM screens with fixed defaults: absolute native coefficient times training SD, retaining the top half (at least two positive-weight genes where available). These are explicit ranking heuristics, not native significance tests, VIP or SVM-RFE. Ridge inner CV preserves patient groups.
- `sb_plots()` includes the new screening weights, captured Ridge coefficient/CV paths and PLS loadings. Existing plotting helpers are reused.
- `surv_pipeline()` accepts PLSR as a short learner name and keeps the old plsrcox name working. The survival Mboost alias explicitly means linear glmBoost, not every mboost engine.

# stormbreaker 0.3.31

- `sb_model_shap()` reuses the approximate permutation SHAP estimator for saved non-tree models. `sb_plots()` now connects LogitBoost and RDA explanations, including grouped importance, beeswarm, waterfall and dependence views.
- `sb_gene_overlap()` returns exact exclusive intersection members and common genes. `sb_plot_gene_overlap()` draws two-to-four-set Venn diagrams (optional ggvenn) or intersection bars without a new plotting dependency.
- `surv_selector()` adds GBM positive native influence linear Cox glmBoost nonzero coefficients, and SuperPC native survival-score screening, with fixed defaults, ordered-chain support and score/training-loss plots. These are explicit package rules, not claims of reproducing unspecified paper methods. Native internal-CV combinations that cannot refit selectors remain excluded.

# stormbreaker 0.3.30

- `sb_plot_calibration()` removes process/status text and explanatory footnotes
  from default plots while retaining full context in returned attributes.
- `bi_run()` and `surv_run()` return an issues table with readable Chinese
  guidance while retaining original failure reasons. `sb_run_issues()` also
  explains previously saved doctor-run results without refitting.
- `surv_run()` no longer attempts an ordinary binary ROC. Time-dependent
  survival ROC remains available through `surv_plot_roc()` at explicit times.

# stormbreaker 0.3.29

- `sb_plot_selector_scores()` and `sb_plots()` now expose saved RSF minimal
  depth/importance, SIS ranks, information-selection scores, LPC scores and
  stability-selection probabilities without refitting.

# stormbreaker 0.3.28

- `sb_algorithm_diagnostics()` adds the native ABESS coefficient/support path,
  automatically included by `sb_plots()` alongside the GIC path.
- `sb_plots()` now draws captured Lasso/Enet selector coefficient paths and
  RF selector OOB histories, including ordered chains and legacy Lasso-StepCox.
  ABESS and CoxBoost selector paths, StepGLM AIC steps and glmBoost training
  loss are also connected to the common plotting entrypoint.
  Capture reuses the original fit; plotting never refits the model.
- `sb_plot_cohort_metrics()` plots C-index or AUC across cohorts.
- `surv_meta_risk()` provides cohort Cox estimates and fixed/random-effect
  meta-analysis for high/low groups or continuous saved-model scores.
- `surv_risk_plots()` adds KM confidence bands and returns group Cox HR/CI/P.

# stormbreaker 0.3.27

- `sb_plot_prediction_confusion()` directly plots aligned binary labels and
  probabilities for apparent, external or OOF predictions, with fixed threshold
  and optional class labels. It is not added to automatic reports.

# stormbreaker 0.3.26

- `surv_selector()` supports fixed-alpha Cox Enet and ordered chains.
- `sb_plot_model_ranking(metric = "c_index")` ranks consistently defined
  C-index estimates without labelling them as AUC.
- `surv_risk_reference()` locks the training median risk threshold and gene
  scaling; `surv_risk_plots()` shares patient order across expression, risk and
  outcome plots and provides KM curves, log-rank P and at-risk counts.
- `surv_plot_roc()` draws cumulative/dynamic IPCW ROC at multiple times,
  preserving unavailable horizons instead of inventing curves.

# stormbreaker 0.3.25

- `bi_selector()` supports Enet (ElasticNet) and GBM in ordered binary chains,
  reusing fold-local selection and fixed defaults. Enet keeps nonzero
  coefficients; GBM keeps positive relative influence, potentially all genes.
- `sb_plots()` includes Enet CV/coefficients and GBM selector importance and
  saved training loss, without refitting.

# stormbreaker 0.3.24

- `sb_plot_shap_importance(groupby = ...)` compares groups using existing SHAP
  values and a common global feature order, without fitting new models.

# stormbreaker 0.3.23

- `bi_run()`, `surv_run()` and `multi_run()` provide task-specific short-name
  entrypoints with fixed algorithm defaults and no automatic browser.
- `bi_selector()` and `surv_selector()` create ordered selector lists for
  `bi_pipeline()` and `surv_pipeline()`. Stage evidence is saved and the whole
  chain refits within supported inner/outer training folds. Existing admission
  restrictions remain; only a small explicit set of three-stage presets is used.
- `sb_models()` lists all task-specific short names; no hyperparameter grid.
- `sb_plots()` gathers saved native and selector plots, including nearest
  centroid profiles and Bernoulli probabilities. Unavailable evidence is visible.
- `sb_tree_shap()` adds seeded permutation interventional SHAP for saved tree
  predictions, explicitly approximate rather than exact TreeSHAP. Existing SHAP
  plots accept these results without relabelling their score scale or method.

# stormbreaker 0.3.22

- `sb_cox_forest()` fits separate gene-level Cox models and returns HRs,
  pointwise Wald intervals, raw/BH P values and a log-scale forest plot.
  All genes and failed fits remain in the table; no significance filtering.

# stormbreaker 0.3.21

- `sb_pipeline()` adds Euclidean nearest-centroid and Bernoulli naive Bayes
  binary learners, with no additional dependencies. Nearest-centroid scores
  use squared-distance softmax (not calibrated risks); Bernoulli NB uses
  strict x > 0 binarization, Laplace smoothing alpha 1 and empirical priors.

# stormbreaker 0.3.20

- `sb_pipeline()` adds binary LASSO, Stepglm, glmBoost and RF importance
  selectors through `selector_options`. Feature selection is refitted inside
  training folds, including glmnet tuning and SVM probability calibration;
  LASSO inner folds respect patient groups. Final fits retain selected genes
  and native coefficient/importance tables without fallback genes.

# stormbreaker 0.3.19

- `sb_bootstrap_calibration()` refits a declared binary pipeline in ordinary
  patient bootstrap samples and saves apparent/optimism-corrected LOWESS curves.
- `sb_plot_dca()` optionally standardizes binary net benefit and adds a
  cost:benefit axis without changing thresholds or predictions.
- `sb_plot_model_ranking()` combines an explicitly ranked cohort heatmap with
  mean-AUC bars; incomplete ranking-cohort data are never silently averaged.
- `sb_plot_nomogram()` optionally overlays reference densities and one patient's
  exact points and probability, without fitting a surrogate model.
- `sb_plot_prediction_roc()` draws apparent-training, external or OOF score ROC
  curves with AUC and optional patient-bootstrap confidence-interval labels.
- `sb_plot_roc()` now pools OOF predictions within each resampling context by
  default; `by_fold = TRUE` retains per-fold plots. Saved aggregate metrics stay
  unchanged. Optional AUC/interval labels describe the displayed curve.

# stormbreaker 0.3.18

- `sb_start()` opens a local Chinese guide with three choices: feature screening,
  model comparison and validation of a saved model. It supplies executable
  examples and workflow-specific dependency instructions without running them.
- `sb_compare_models()` now defaults to three complementary models, rather than
  the entire catalog. Explicit `models = NULL` retains the previous all-model
  behavior. `sb_compare_multiclass()` also has a three-model default. Models are
  not selected by observed performance, and missing engines are not replaced.
- `sb_clinician_views()` defaults to key ROC, calibration and saved decision
  curves. Use `detail = "full"` for precision-recall and probability distributions;
  explicitly requested diagnostics and explanations remain available in both.

# stormbreaker 0.3.17

- `sb_analyze_methods()` additionally supports multiclass mRMR, JMI, CMIM,
  DISR and JMIM, with saved selected genes and native sequential-score plots.
  Native early stopping is retained; missing genes are never filled in to
  reach the requested maximum subset size.
- Multiclass adapters have been added for EDDA/MclustDA, HDDA/HDA, RLT,
  aorsf/wsrf, multinomial logistic/stepwise/bias-reduced regression, GAM
  and PRE. Experimental SLOPE is explicitly unsupported while its multiclass
  validation remains incomplete; it is not part of this release's usable set.

- Multiclass prediction now also supports FDA/MDA/PDA, stabilized/shrinkage/
  local discriminants, PAMR, Gaussian processes, kernel-weighted KNN and two
  LiblineaR logistic solvers. FDA's single-row native polynomial-basis issue
  and local LDA's missing single-row class names are handled without refitting.
- Native multiclass objectives are available for XGBoost tree/DART/linear,
  LightGBM GBDT/DART/GOSS/RF and CatBoost. Boosters retain portable native model
  bytes or strings rather than live external pointers. Saved native importance
  is available for tree variants; gain definitions are kept engine-specific.

- `sb_analyze_methods()` and `sb_compare_multiclass()` additionally support
  ridge, five e1071 SVM variants, kernel naive Bayes, RDA, softmax nnet,
  deterministic KNN, ranger/ExtraTrees, C5.0 trees/rules, conditional trees/
  forests, bagging, regularized random forest, tree and evtree. SVM probabilities
  use the native one-versus-one coupling; KNN returns all class vote fractions.
  These methods do not fabricate gene-selection sets. SVM support-vector counts
  and ranger importance plots are available; other native diagnostics are not
  implicitly claimed available merely because prediction is supported.

- `sb_analyze_methods()` supports native multiclass Boruta, varSelRF and VSURF
  selection with saved gene sets, importance/history plots and no hidden
  classifier. Empty multiclass Boruta selections remain empty, not replaced.
- `sb_analyze_methods()` and `sb_compare_multiclass()` support ABESS multinomial
  GIC-selected subsets, class coefficients, native probabilities and a subset
  selection curve. Class ordering follows the input factor rather than the
  incorrect coefficient column labels in ABESS 0.4.11.

# stormbreaker 0.3.16

- `sb_analyze_methods()` and `sb_compare_multiclass()` now support native LDA,
  QDA, Gaussian naive Bayes and rpart classification trees. Class statistics
  are not relabelled as gene selection; trees expose split usage separately
  from importance. Singular fits remain failures without fallback learners.

# stormbreaker 0.3.15

- `sb_analyze_methods()` supports native multiclass random forest, multinomial
  LASSO and Elastic Net. Saved class-specific gene results are available through
  the optional `class` argument in feature-query functions.
- `sb_compare_multiclass()` evaluates the three native methods on shared outer
  folds with inner glmnet tuning and patient-group separation when declared.
- `sb_multiclass_task()`, `sb_predict_multiclass()` and `sb_multiclass_metrics()`
  preserve explicit class mappings and provide class probabilities, per-class
  metrics and macro summaries. Confusion plotting is explicit only through
  `sb_plot_multiclass_confusion()`; existing report defaults are unchanged.

# stormbreaker 0.3.14

- `sb_analyze_methods()` independently runs explicitly named Boruta, varSelRF,
  VSURF, LASSO, Elastic Net, MCP, SCAD, ABESS and random-forest analyses, returning
  saved plots, native evidence and gene results without changing clinical-report
  defaults. Failures remain method-specific; no files are exported automatically.
- `sb_method_catalog()` shows task support and selection versus ranking.
- `sb_method_features()` and `sb_selected_features()` retrieve saved gene results
  without refitting. Random forest provides ranking, not a fabricated subset.

# stormbreaker 0.3.13

- Clinical figures use readable model names and short context labels without
  internal IDs, successful-run labels or unnamed-cohort placeholders. Model
  settings stay in the report body; normal availability rows are hidden.

- Confusion matrices separate model/context cells instead of overlaying counts,
  and scale report height with panel count. Threshold panels use shorter labels.
- Visual reports place explanatory text in a separate body section rather than
  small-print captions under each figure.

# stormbreaker 0.3.12

- `sb_clinician_views()` no longer accepts the candidate-only `fit`, `newdata`
  and `complete_native` shortcuts. Technical diagnostics must be extracted
  explicitly and passed through `native` or `selectors`; default and expanded
  reports do not compute them. Individual diagnostic functions remain available.

# stormbreaker 0.3.11

- `sb_plot_clinician_panel()` combines explicitly selected clinical figures with
  panel letters and separate captions, retaining original data and legends.

- `sb_clinician_views(fit = ...)` discovers native diagnostics from a saved model
  alongside saved evaluation results. Detailed extraction and query rows are
  explicit options; report assembly does not retrain the model.

- `sb_clinician_views()` keeps explanations in separate figure captions for all
  topics, so its plots stay uncluttered in both direct exports and HTML reports.

- Cancerclass exposes original centroids and its existing score-to-probability
  mapping, explicitly distinct from a calibration-reliability assessment.

- Local LDA exposes explicit-query native neighborhoods, distance weights,
  tied neighbors and native class-count exclusions without retuning.

- HDA and HDDA expose retained native subspaces, HDA coordinate densities,
  HDDA class spectra and original-unit centers without recomputing PCA.

- EDDA and MclustDA expose original component means, fitted correlations and
  within-class mixing weights without reselecting or refitting components.

- SDA exposes original shrinkage intensities and pre-rounding log-odds
  coefficients; SLDA exposes original subspace loadings and class centers.

- FDA, PDA and MDA expose original discriminant centers and explicit-query
  variates; MDA also exposes within-class latent mixture weights.
  Singleton FDA interaction queries retain native polynomial dimensions,
  producing the same result as the corresponding row in a multi-row query.

- Gaussian-process classifiers expose original basis weights and explicit-query
  kernel contributions, not feature SHAP or unsupported posterior intervals.

- Random ferns expose original split-use counts and real-valued threshold
  distributions, explicitly separated from importance and clinical cutoffs.

- Rotation forests expose all original rotation matrices, tree sizes and
  explicit-query per-tree probabilities. Rotated axes are labeled separately
  from original variables; loading signs are not feature effects or importance.

- C5.0 tree and rule ensembles expose captured native usage/split summaries
  and trial sizes. Unparseable importance in unsplit ensembles is explicitly
  unavailable; compact models never fabricate missing summaries.

- The tree-package classifier exposes every original node/split and terminal
  class probability. Native split precision, internal-to-original feature
  mapping and training-node counts remain intact; no pruning or refit.

- Evolutionary classification trees expose their original structure and node
  memberships. Native split indices map through the retained predictor schema,
  not fitted membership columns; no search, pruning or surrogate fitting occurs.

- Binary, discrete-time survival and AFT BART expose retained posterior
  prediction draws for explicit queries, with means and pointwise equal-tail
  credible limits. AFT retains latent log-time location, not outcome prediction
  intervals. Native discrete survival grids render as steps; no new MCMC.

- Bagged survival KNN exposes original bootstrap multiplicities/coverage and,
  for explicit queries, every native bag's survival curve with the ensemble
  mean. Bag variation is not a confidence band and bootstrap coverage is not
  OOB performance. No new bootstrap, fitting or calibration is performed.

- Ordinary and kernel-weighted binary KNN expose explicit-query neighborhoods:
  original training row, distance, class and native/normalized voting weights.
  Native bundles accept optional query data; without it KNN reports unavailable
  instead of silently evaluating training rows. No selection, tuning, prediction
  mapping or self-neighbor behavior changes.

- Robust Mqle logistic regression supports its original coefficients,
  faithful nomograms and exact additive log-odds SHAP. Robust estimation is
  preserved; no ordinary logistic surrogate or missing robust weights are
  reconstructed.

- LIBLINEAR L1 and L2 logistic models support original-unit coefficients,
  faithful nomograms and exact additive log-odds SHAP. Native class orientation,
  bias and training standardization are retained; no surrogate fit or ordinary
  unpenalized confidence interval is introduced.

- Bagged classification trees expose primary split counts, complete tree sizes
  and native rpart importance for every bag/feature. Importance is averaged
  across all bags only for display and may credit surrogate splits; disabled
  OOB evaluation is not reconstructed.

- Binary and survival conditional inference forests expose complete native
  split usage and tree sizes. Original model-frame indices correctly distinguish
  predictors from response columns; structural counts are not conditional
  permutation importance or newly computed statistical tests.

- Binary and survival accelerated oblique forests expose original leaf-count
  distributions and recorded OOB scores. AUC-ROC and Harrell's C-index retain
  their native meanings; final-only evaluations stay single points rather than
  invented learning curves. No R6 model updates or importance computation occur.

- Binary and survival reinforcement-learning forests expose complete native
  split-variable counts and tree-size distributions from compact state.
  Structural usage is not relabeled as variable importance; disabled importance
  and missing reinforcement/OOB histories are not reconstructed.

- Weighted-subspace forests expose original information-gain-ratio importance,
  total tree-node counts and individual-tree OOB errors. Original metrics are
  captured separately before prediction compaction; zero OOB placeholders are
  never treated as observed performance. Constant bridge features are identified.

- Regularized random forests now use the native forest diagnostic interface
  for original Gini importance, every tree's terminal-node count, and captured
  OOB error history. Zero-importance features and missing OOB values remain;
  fitting, regularization and prediction state are unchanged.

- Binary and survival PAMR models expose complete before/after centroid maps
  at the original fixed threshold, retaining raw units and native SD-scaled
  deviations. Native survival groups are not patient subgroup analyses or
  horizon probabilities; no threshold search or refitting is performed.

- Gehan elastic-net AFT models expose original coefficients and optionally
  captured full coefficient/CV paths. Input-unit paths and centering intercepts
  exactly preserve saved normalization; native standardized values remain
  available. Selection scores are not AUC or external-validation estimates.

- Buckley-James linear boosting exposes original log-time coefficients and,
  when already saved for horizon prediction, the residual Kaplan-Meier event
  steps. This training-derived residual distribution is not held-out calibration
  or a patient survival curve; no histories or uncertainty are reconstructed.

- Lin-Ying, LASSO and SSCAD additive-hazard models expose original input-unit
  coefficients. Opt-in diagnostic capture retains native penalized coefficient
  and BIC paths with the original selected lambda, without retaining the full
  native training object or changing compact prediction state. These effects
  are not log hazard ratios or the separate Cox-Breslow probability calibration.

- Aalen models expose unchanged cumulative feature effects and intercept
  processes; proportional-odds models expose native coefficients and baseline
  cumulative odds. Native steps and negative values are retained without
  smoothing, differentiation or invented uncertainty.

- MTLR exposes original feature/time-point weight maps and separate biases.
  Saved normalization, full effective time grid and exact original-input-unit
  equivalent weights are retained; coupled potentials are not mislabeled as
  independent logistic effects or hazard ratios.

- k-TSP exposes all original ordered pair rules and native pair-selection
  scores, with explicit half-votes for equality and unweighted vote averaging.
  Pair scores are not vote weights, validation accuracy or calibrated risks.

- Saved VSURF selectors expose original mean/SD importance and three-stage
  subset counts. Clinician reports accept optional named selector diagnostics
  in a separate Feature selection section, without refitting a selector or
  presenting selection evidence as validation performance.

- Six flexible parametric/spline survival options now expose original location
  coefficients on their actual parameter scales. Spline models additionally
  expose the original zero-predictor baseline within saved boundary knots,
  explicitly not a representative-patient or average survival curve.

- Original-model nomograms now support all four survreg AFT models, retaining
  their negative-log-time ranking score and native distribution/scale at the
  saved horizon. No Cox baseline or surrogate is used; without a saved horizon
  the output retains a labelled score axis rather than invented probabilities.

- The four survreg AFT options expose original log-time coefficients and
  conditional time ratios, explicitly distinct from hazard ratios. Compact
  parameters suffice; intercept and distribution scale remain separate.

- Conditional inference trees now expose every original branch and training-node
  membership count from compact binary and survival fits, retaining exact splits
  separately from rounded display labels. No pruning or new node tests are run.

- MARS/earth now exposes the original selected logistic basis coefficients and
  captured basis-pruning GCV. The separate least-squares coefficients are not
  mistaken for logistic slopes; degree-one MARS also supports exact nonlinear
  additive SHAP. Exact knots/directions remain available apart
  from rounded basis labels. GCV is not labelled as held-out logistic performance.

- nnet and multilayer neuralnet expose original topology and weight matrices
  from compact state, retaining bias terms, input scaling and nnet's native
  saturation semantics. No feature-importance, surrogate nomogram or invented
  training-history claim is made.

- randomGLM now exposes per-bag feature inclusion, original bag coefficients and
  bag-size distributions from compact state. Zero/omitted terms and intercept-only
  bags remain explicit. Probability averaging is not replaced by coefficient
  averaging or a surrogate nomogram; bag variation is not treated as a CI.

- SLOPE now retains its original complete CV summary only in opt-in captured
  diagnostics, including native early-stopped paths and native t-based bounds.
  Compact prediction state remains unchanged. Exact original SLOPE coefficients,
  nomograms and additive SHAP are available without capture or surrogate refits.

- All four LVQ variants now expose original learned prototype profiles and
  squared-distance maps from compact saved state. Training scaling and native
  prototype order/classes are preserved, with original-unit coordinates in the
  table; no patient prototypes, refitting or invented learning curves.

- Original LDA, QDA and regularized discriminant class means are available from
  compact saved models. QDA/RDA additionally expose prediction-covariance
  correlations, with both RDA regularizations retained. LDA now supports exact
  original-model coefficients, nomograms and additive SHAP. No covariance is
  invented from rank-reduced LDA scaling and no additive surrogate replaces
  quadratic QDA/RDA predictors.

- Exact coefficients, nomograms and additive SHAP now support stepwise GLM,
  GLMBoost, plsRglm, bias-reduced logistic and Cauchy-prior logistic models.
  Centered intercepts are converted to original input units, without replacing
  fitted slopes or penalties. Automatic diagnostics also reuse saved linear
  coefficients where no more specific registered native view exists.

- `sb_native_diagnostics()` detects implemented native adapters for a saved
  learner, returning the first available view by default or all registered
  candidates on request. Each extraction runs against the actual saved state;
  missing capture, unavailable views, unimplemented learners and errors stay
  distinct. Clinical report assembly accepts this bundle and retains reasons.

- Clinical reports also accept saved exact SHAP and nomogram objects. Global
  explanations and explicitly named dependence views share a topic; feature
  values remain ID-aligned, no explanation is recomputed, and no patient-card
  or default mass per-feature export is introduced. Figure heights adapt to
  the number of displayed explanation terms.

- `sb_clinician_views()` assembles saved ROC/PR, calibration, decision curves,
  probability distributions and explicitly requested threshold/confusion views,
  plus supplied native diagnostics. `sb_render_clinician_report()` exports them
  by topic with an availability table. Missing curves are not computed; partial
  contexts and extraction errors remain explicit. No SHAP surrogate or refit.

- The offline visual report can group saved plots into expandable topics,
  showing the first topic by default while keeping every figure available.
  Compact images now move original titles and explanatory labels into selectable
  figure legends; additional diagnostic captions can be supplied separately.
  No JavaScript, model refitting or automatic omission of detailed plots.

- GBM-Cox with a one-tree maximum now fits its sole candidate instead of
  hitting native GBM's one-row CV reduction error. Tuning explicitly records
  that no tree-count search ran; requesting a CV curve reports unavailable.
  Outer performance evaluation and multi-candidate searches are unchanged.

- `sb_naive_bayes_diagnostics()` adds original Gaussian and kernel class-density
  curves from compact fitted state, with no refitting or extra capture. Gaussian
  log-likelihood contributions reconstruct original posterior log odds when
  combined with the saved prior. Feature densities keep their own units; kernel
  likelihood ratios are explicitly unavailable rather than approximated.

- `sb_survival_svm_diagnostics()` supports all twelve existing survival-SVM
  formulation/kernel options: native pre-calibration score and kernel-row weight
  distributions, plus original-unit coefficients for the four linear variants.
  Censoring is not misrepresented as a binary outcome and kernel-row weights are
  not labeled feature or patient importance. Original artifacts require capture.

- Binary margin diagnostics now also cover three LIBLINEAR SVC variants, DWD
  and radial LS-SVM. Preserve original-unit linear weights and the exact saved
  score mappings; LS-SVM uses its two-class kernel-score contrast and does not
  invent sparse support-vector counts or linear input coefficients.

- `sb_svm_diagnostics()` adds native support-vector counts for five binary
  LIBSVM options, captured development decision distributions, and exact
  original-unit linear weights. Class orientation follows the native model;
  decision scores are not mislabeled as probabilities or geometric distances.

- PLS-Cox component-count CV is available from saved tuning without capture.
  Label native criterion 5 as integrated AUC (higher is better), not error.
  Preserve native spread without calling it conventional SE or a confidence
  interval; one-predictor fits explicitly have no CV curve.

- `sb_pls_diagnostics()` adds original loading/weight heatmaps and sample-score
  distributions for plsRglm, PLS-DA, sparse PLS-DA and fixed/CV plsRcox. Capture
  retains native binary component summaries without changing portable models.

- `sb_superpc_diagnostics()` displays captured native screening scores,
  singular-vector loadings and complete CV threshold/component results.
  Keep the original one-component selection and effective fold count, without
  mistaking native log-scale standard-error bands for confidence intervals.

- `sb_boosting_diagnostics()` adds captured hinge-loss boosting risk, component
  updates and original-margin coefficient paths. Default portable prediction
  state stays unchanged; historical compact fits do not get invented histories.

- `sb_boosting_diagnostics()` adds AdaBoost/LogitBoost split usage, AdaBoost
  voting weights and actual tree sizes, and LogitBoost stump sequences with
  exact thresholds and signed contributions. Native early convergence does
  not generate artificial extra iterations; no additional fitting is needed.

- `sb_boosting_diagnostics()` adds captured CoxBoost coefficient and nonzero
  paths, and the CV variant's saved partial-log-likelihood search without
  requiring capture. Preserve the original stopping step and prediction scale.

- `sb_boosting_diagnostics()` adds saved LightGBM/CatBoost importance, actual
  LightGBM tree-size distributions and native CatBoost interaction scores.
  Interaction display limits never discard the complete table, and configured
  rounds are not mistaken for the actual number of trees retained by LightGBM.
  Single-feature CatBoost models return no distinct pairs without triggering
  the native empty-result error; other native extraction failures retain their
  cause and never become invented zero interaction scores.
- `sb_boosting_diagnostics()` also reads GBM final-fit training errors, original
  GBM-Cox CV searches and native GBM/XGBoost feature scores. Keep linear XGBoost
  coefficients distinct from tree importance and retain zero-use features.
- Fix the final GBM-Cox refit when CV selects one tree: use the selected tree
  count without a redundant second CV that triggers a native gbm dimension error.
- `sb_boosting_diagnostics()` exposes native mboost in-bag loss sequences and
  component-update counts without refitting or changing stopping. Preserve
  objective signs, unused terms, intercepts and original feature labels, and
  distinguish component counts from predictive importance or selection stability.
- `sb_forest_diagnostics()` adds native importance and terminal-node-count
  distributions for ranger/randomForest, tree sizes for randomForestSRC, and
  captured native OOB results for randomForest/randomForestSRC. Disabled or
  unrecorded diagnostics are not reconstructed; a single OOB value stays a point.
- `sb_rule_ensemble()` exposes original PRE rule definitions, coefficients and
  coefficient-only validation paths directly from saved prediction state.
  Keep zero terms and distinguish native slopes from feature importance, and
  coefficient validation from validation of the complete rule-learning pipeline.
- Add explicit-threshold metrics and confusion plots, aggregate probability
  distributions, and ID-aligned SHAP dependence plots. Keep cohorts and folds
  separate, preserve unavailable groups, and never optimize a test-set threshold.
- Add original linear-model nomogram conversion with exact numerical prediction
  recovery. Zero coefficients stay in the model/table; their axes are optional.
  Short-axis labels are staggered and probability ticks are spaced for reading.
- Add opt-in native final-fit diagnostic capture without additional training.
  Captured artifacts can contain sensitive training data and engine objects.
- Add native regularization-path and ABESS GIC support-selection views. These
  development-fit diagnostics are not external validation performance estimates.
- Add saved glmnet inner-validation and nonzero-count paths, and retain native
  MCP/SCAD coefficient/BIC sequences for optional final-fit capture.
- Add Boruta importance-range, hit-rate and opt-in iteration-history plots,
  including shadow references and before/after rough-fix decisions. Preserve
  native post-rejection sentinels as raw evidence, not plotted numerical scores.
  Add varSelRF's saved OOB elimination history without relabelling it as CV.
- Add original mgcv GAM smooth-effect plots for binary and Cox models, with
  empirical training support and pointwise uncertainty conditional on fitted
  smoothing parameters. Preserve centered link-scale interpretation.
- Add complete native rpart tree diagrams for binary and survival learners,
  including node counts, incoming rules, terminal predictions and exact split
  tables. Preserve every branch and distinguish relative rates from probabilities.

# stormbreaker 0.3.10

- Rank all input genes from saved outer-fold selector counts, deduplicating
  shared selector configurations across downstream learners. Keep glmnet
  nonzero counts separate; report denominators, zero counts and unavailable
  records explicitly. Export rankings in the example development workflow.
- Simplify model reports with compact plots and optional gene-selection
  rankings. Retain units, intercepts, endpoint markers and sample controls.
- Remove repeated explanatory boilerplate from report pages and galleries.

# stormbreaker 0.3.9

- Unify SHAP views with a white, restrained blue/red theme. Show numerical
  importance labels, clearer within-feature colour legends and row guides.
- Use directional waterfall arrows, neutral background/final-score markers
  and a separate number column that keeps tiny contributions readable.
- Adapt offline SHAP plot heights to the displayed feature count and clarify
  the different colour meanings in beeswarm and waterfall views.
- No changes to fitting, attribution, predictions or survival-curve calculation.

# stormbreaker 0.3.8

- Inspect and plot saved linear coefficients in original units or per original-
  training standard deviation, retaining intercept and score-scale semantics.
- Recover full Breslow/Efron baselines from the original development task at
  fixed saved Cox coefficients. Check the saved horizon and predict individual
  step curves without reading new-patient outcomes or extrapolating follow-up.
- Render offline interactive model explanations with coefficient units, existing
  global SHAP views and a linked individual waterfall/survival sample selector.
  Existing fitted parameters, model selection and prediction behavior are unchanged.

# stormbreaker 0.3.7

- `sb_compare_cohorts()`, `sb_compare_cv()` and `sb_run_cv()` can show actual
  model, outer-fold and evaluation progress without changing predictions.
- `sb_read_cohort()` reads expression and clinical CSV files, preserves text
  sample IDs and joins by exact ID with clear mismatched-sample errors.
- `sb_select_model()` chooses a training-only fitted model by equally weighted
  validation ROC AUC or Harrell C, excludes test/external cohorts and returns
  an object ready to save and predict without extracting internal wrappers.

# stormbreaker 0.3.6

* Compare development, validation and test cohorts in one results matrix with
  `sb_compare_cohorts()` and `sb_view_results()`.
* Display standard algorithm names while retaining full configuration IDs
  and parameters in details and exports.
* Open results automatically after high-level comparisons in interactive R;
  save standalone HTML pages for scripted runs.
* Add multi-cohort metric views, validation-only summaries and CSV/SVG export.
  Native model implementations and saved metric calculations are unchanged.

# stormbreaker 0.3.5

## Initial public release

* Build binary and right-censored survival prediction workflows from a
  user-supplied gene set.
* Compare task-compatible models on shared folds with `sb_compare_models()`;
  use `sb_pipeline()` for feature-selection and learner combinations.
* Reuse identical fold-local Boruta and varSelRF selections across learners.
* Fit a final model, predict new patients and evaluate independent cohorts.
* Explore discrimination, calibration, decision curves, feature stability,
  cohort-level summaries and additive SHAP through saved-result interfaces.
* Browse leading cross-validated models in an offline, CSV-driven leaderboard.
* Compute FDA and penalized FDA posteriors on the log-odds scale to avoid
  numerical underflow while preserving their fitted models.
* Restore saved XGBoost models with both the 1.x and current R interfaces;
  write newly saved models in XGBoost's JSON format.
* Add English and Chinese quick starts, installation guidance, a complete
  license and a GitHub source/release layout.
- Added binary and survival PAMR selectors using native nonzero shrunken
  centroids at fixed threshold 1, with fold-local selection and no gene padding.
