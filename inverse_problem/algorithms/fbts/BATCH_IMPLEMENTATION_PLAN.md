# FBTS batch implementation plan

Scope agreed on 2026-09-24: acquire the full initial RX dataset once per target scene, select one random representative of each antenna-layout rotation/reflection class, and reconstruct each selected configuration with its own normal forward, adjoint, and sensitivity simulations. All active antennas transmit in ascending original-index order, starting with the first available antenna. The stages below are implemented; see VALIDATION.md for measured validation results and the remaining parallel-execution validation limit. The user subsequently confirmed that reflections are equivalent too; the default equivalence group is dihedral.

A **scene** means the complete fixed material distribution, potentially containing several target objects. Changing an object's shape, position, orientation, or material creates a different scene. Merely keeping the same DOI does not make two scenes share measurements.

1. **Scientific assumptions and limits**

   The current C solver adds source samples at antenna grid cells and reads Ez at receiver grid cells. It does not insert an antenna body or load into the material map merely because that antenna appears in the list. For a fixed scene, source waveform, and transmitter position, removing passive probes and zero-amplitude sources leaves the retained RX traces unchanged. Thus a full scan can supply measurements for reduced arrays. This is a source-code conclusion; numerical comparison is a required implementation gate below.

   The full scan consists of N separate FDTD runs, one transmitter at a time. Its tensor has dimensions `[original_tx, original_rx, time]`. For ascending active indices `A`, obtain `rx_full(A,A,:)`. Remove disabled transmitters and receivers from the reconstruction configuration as well. Do not pass zero-valued observations for disabled channels into the ordinary residual calculation.

   Every selected configuration starts from the same prescribed background estimate, with a fresh conjugate-gradient state. It still performs four FDTD calls per active transmitter per iteration: forward, adjoint, coarse sensitivity baseline, and coarse sensitivity perturbation. Do not warm-start one configuration from another configuration's reconstruction.

   Rotation/reflection classes describe the nominal circular antenna layout. They are not proof of equal images or equal reconstruction quality for a fixed asymmetric/off-center scene. Exact physical rotational equivalence would require rotating the entire experiment, including its material distribution and other non-invariant components. Rounded antenna coordinates, a Cartesian grid, and rectangular PML further limit numerical rotational symmetry. The dataset must identify this as sampling one orientation per layout class.

   The current implementation estimates only relative permittivity. Conductivity and permeability are held fixed in the inversion; copying nonzero true conductivity into the forward configuration would give it known information, not reconstruct it. The first batch implementation should validate the existing lossless, nonmagnetic setting. Broader known-property cases require an explicitly recorded model assumption; estimating additional properties is a separate algorithm change.

   Measurement reuse is not generally equivalent to physically removing modeled antenna bodies or changing their loads: those changes can alter the electromagnetic problem. If such antenna physics is added later, invalidate this assumption and revisit the acquisition model.

2. **Exact rotation/reflection-class rule**

   Let sorted disabled antenna indices be `d1 < ... < dk` on an N-position nominal uniformly spaced ring. For `k > 0`, form the clockwise gap tuple

   `g = [d2-d1, ..., dk-d(k-1), N+d1-dk]`.

   Gaps are counts of angular index steps, not rounded Euclidean distances. Their sum is N. The canonical key is `(N,k,lexicographically_smallest_cyclic_shift_of_g_or_reverse(g))`. For `k=1`, the gap tuple is `[N]`; `k=0` has a separate full-array key. Exclude `k=N`, which leaves no transmitter.

   Identify gap tuples under cyclic shifts and reversal, while retaining their adjacency. For example, on 12 positions, disabled sets `[1,2,4,7]` and `[1,2,5,7]` have gaps `[1,2,3,6]` and `[1,3,2,6]`. They have the same sorted gaps but are different rotation/reflection classes. Include the reversed tuple and all its cyclic shifts when choosing the canonical key. Arbitrary permutations of four or more gaps remain distinct unless a rotation or reflection relates them.

   For each class, construct its distinct rotations and reflections of the disabled-index set and select one uniformly. Symmetric layouts can have fewer than 2N distinct transforms. Store the distinct orbit size, canonical disabled set, selected disabled set, and applied index shift. These values specify the whole class without saving a long duplicate list.

   Use dedicated reproducible random streams for scene generation and representative selection. Choose representatives independently per scene using stable `(master seed, scene ID, class key)` inputs. Precompute and persist these choices before simulations, so resume and worker ordering cannot change them. Record RNG algorithm and state as well as seeds. For supplied existing configurations whose generation state is unknown, save their exact maps and mark the seed as unknown.

   Build one class catalog per compatible nominal array. Validate cyclic ordering, unique positions, nominal uniform angular spacing, and the correspondence between nominal and rounded positions. Check spacing in physical coordinates: when dx differs from dy, the current index-space circle becomes a physical ellipse. Reject the circular-symmetry assumption in that situation unless geometry is generated as a physical circle. Do not silently apply cyclic-index pruning to an arbitrary nonuniform array. Retain an exhaustive mode for validation and orientation studies.

   Run the full-array baseline first, followed by increasing disabled count and a stable canonical-key order. At every case, sort active original indices; local transmitter 1 is then the first available original antenna. Preserve all active transmitters and, initially, the current inclusion of self-receiver channels. Correct the conflicting comment in build_cfg.m that says only other antennas receive.

   Exhaustive independent enumeration gives:

   | N | Nonempty active subsets before pruning | Rotation/reflection classes, excluding the empty active array |
   | --- | --- | --- |
   | 4 | 15 | 5 |
   | 8 | 255 | 29 |
   | 12 | 4,095 | 223 |

   For N=12, class counts for disabled counts 0 through 11 are `[1,1,6,12,29,38,50,38,29,12,6,1]`. The gap-key classification was independently checked against explicit rotation orbits for all nonempty active subsets at N=1 through 12. This verifies the combinatorial design, not reconstruction equivalence under rotation.

   Orbit size also matters for later statistics. An unweighted mean over classes and a mean over all antenna subsets are different quantities. Save orbit sizes to support the appropriate analysis; one sampled orientation does not measure within-class orientation variability.

3. **Stage 1 — separate acquisition from inversion and prove equivalence**

   Refactor the source-pulse and temporal-weight setup out of runFbts.m into a shared helper, then extract its synthetic-measurement loop into `acquireFbtsMeasurements`. Keep the existing two-argument runFbts call working for the demo. Add a third optional validated measurement input, conceptually `runFbts(caseCfg, options, measurements)`, that skips acquisition when supplied.

   Define the input as an RX array plus metadata containing original TX/RX IDs, time samples, pulse identity, scene/system identity, and an acquisition fingerprint. Validate ordering and dimensions explicitly, including singleton arrays; avoid squeeze operations that lose axis meaning. Reject mismatched time grids, source waveforms, scenes, positions, incomplete scans, and nonfinite data.

   The acquisition fingerprint covers all inputs that affect measured signals: actual material arrays, grid and timestep, boundary/PML configuration, antenna coordinates and TX/RX order, source samples and excitation convention, any initial field state, and solver identity. It must not depend on cosmetic labels or reconstruction iteration count. A different scene, medium, waveform, grid, or solver requires a new acquisition. A different subset or iteration count does not.

   Acquire all N transmitters once per scene and save the acquisition before starting any reconstruction. A helper `selectFbtsMeasurements` selects TX and RX axes consistently with a reduced case configuration. The existing subset construction already updates positions, source dimensions, and local transmitter indices; preserve that behavior and add explicit original/local mappings.

   Acceptance gate: on a small real-MEX case, compare fresh reduced-array scans to slices of the full scan for the full array, singleton arrays, noncontiguous subsets, and subsets excluding original antenna 1. Compare several complete FBTS iterations using fresh versus supplied data, including image, costs, gradients, step sizes, restart flags, and sensitivity results at factors 1 and 4. Exclude runtime fields from equality checks. Use a documented numerical tolerance; report exact equality if observed. Also compare a reduced adjoint to a full-list adjoint whose disabled residual sources are zero. Failures block rollout of shared acquisition.

4. **Stage 2 — implement reproducible rotation pruning**

   Introduce pure helpers for canonical keys, class enumeration, and representative selection. Replace the batch's `2^N-1` case allocation and direct nchoosek traversal with the saved catalog. Keep exhaustive mode available.

   For N=8 or N=12, straightforward enumeration with one retained canonical class record is sufficient. Avoid retaining every equivalent subset as a large matrix. For larger arrays, use a generator/streaming approach or direct binary-necklace enumeration; rotation pruning still leaves an exponentially growing number of classes, so do not imply that it makes arbitrary N inexpensive.

   Acceptance gate: verify all counts above, class membership by explicit rotation, gap-order counterexamples, reflection merging and non-dihedral permutation separation, periodic layouts, k=0, k=1, k=N-1, N=1, and exclusion of the no-active-antenna case. Verify deterministic selection under resume and reordered processing, and that every selected representative belongs to its class. Do not write a test requiring rotated configurations to reconstruct identical fixed scenes.

5. **Stage 3 — connect the scene, acquisition, and case loops**

   Keep fbts_batch.m as a convenient entry script and move orchestration into a callable `runFbtsBatch` function suitable for headless execution. Separate batch controls from ordinary FBTS controls rather than making unrelated saving/scheduling options part of the numerical algorithm.

   Support the existing single cfg as a one-scene batch and an explicit collection of scene configurations for multiple scenes. Refactor build_cfg.m into a wrapper over a parameterized builder so generation can accept and preserve target options, a scene seed, and resolved defaults. Preserve an imported cfg's exact contents without regenerating its targets.

   The order is: resolve systems/scenes and options; build catalogs; persist selected case definitions; for each scene acquire/load its full RX dataset; for each selected class create the reduced configuration; start FBTS independently; save results and status. Default to one scene and serial execution, preserving the current scope until a larger collection is requested.

   Suggested batch controls: symmetry mode (`dihedral` by default, `rotation`, or `none`), master seed, scene inputs/generation specification, disabled-count range defaulting to `0:N-1`, headless/figure choice, resume policy, iteration checkpoint interval, and optional worker limit. Keep all iterations independent between configurations and keep the original base cfg unchanged.

   Acquisition durations belong to the scene, once. Each case records its cache reference, cache hit, and zero new acquisition calls. Update addSolverStatistics so an empty measurement-duration array gives count and total zero, with first/min/max/mean/median/variance unavailable; it currently accesses seconds(1). Preserve raw per-call timings for inversion phases. Batch wall time is measured directly, particularly when workers run concurrently.

   Acceptance gate: verify one N-transmitter acquisition per scene, the expected number of configuration reconstructions, unchanged source cfg, ascending original transmitter order, no repeated full-array baseline, shared source/time settings, and separate acquisition for two different scenes. Verify cache invalidation and incomplete-cache rejection.

6. **Stage 4 — save an analysis-ready dataset**

   Use MAT v7.3 files. Store small related metadata in scalar structs and struct arrays, and large numerical arrays as independent top-level variables. MATLAB matfile supports efficient partial access to v7.3 arrays but cannot index directly into fields of stored structure arrays. Keep portable relative references and content hashes; provide `loadFbtsBatchCase` to assemble one case with its shared system, scene, and acquisition.

   Proposed layout:

   ```text
   batch_XXXX/
     manifest.mat
     cases.csv
     provenance/
     systems/system_0001.mat
     scenes/scene_000001/
       truth.mat
       acquisition.mat
       cases/class_.../
         result.mat
         status.mat
         checkpoint.mat       # latest complete iteration, during work
         failure.mat          # if applicable
         figures/             # optional, off for server batches
   ```

   A scene references one system; multiple systems are allowed when grids, media, arrays, or other common parameters differ. System and scene data are saved once instead of copied into every result. The entire directory is the portable dataset; a case file alone requires its referenced files. The loader can later support an explicit standalone export.

   **manifest.mat — small metadata and catalog**

   | Variable/type | Saved contents |
   | --- | --- |
   | `batch` — scalar struct | Schema version, batch ID, timestamps, resolved batch and algorithm defaults, requested scene list, execution settings, master seeds/RNG algorithms, dataset conventions, source snapshot reference, and generation specifications. |
   | `classCatalog` — struct array | Per-system N, disabled count, canonical ordered-gap key, canonical disabled indices, orbit size, and rotation-and-reflection policy. |
   | `caseCatalog` — table or struct array | Stable case/scene/system/class IDs, randomly selected original active and disabled indices, selected rotation, RNG provenance, requested options, relative output paths, and fingerprints. |
   | `cases.csv` — separate index | One row per case: identifiers, disabled/active counts and IDs, class key/orbit size, completion/failure state, completed iterations, relative paths, and timing summaries. No image-quality metrics. Treat per-case status records as authoritative if interrupted before index refresh. |

   **systems/system_<hash>.mat — the common system model**

   | Variable/type | Saved contents |
   | --- | --- |
   | `model` — scalar struct | Nx, Ny, sizeZ, dx, dy, dt, Nt, deltaF, grid origin and units, coordinate/array-axis conventions, field staggering, snapshot settings, initial-state convention, boundary conditions, and full resolved PML type/settings/profile arrays. Include configured constants and effective solver constants. The present C solver sets eps0/mu0 internally; saving cfg constants alone would miss that distinction. |
   | `model.antennas` — struct | Original IDs and cyclic ordering, nominal angles/center/radius, exact rounded grid positions and corresponding physical coordinates, array construction metadata, and self-receiver inclusion policy. |
   | `model.sensitivity` — struct | Resampling method/factor, resolved coarse grid and PML, mapped antenna coordinates, coarse time grid, source samples, temporal weights, interpolation conventions, and collision checks. Case-specific overrides are saved in the case's resolved parameters. |
   | `x_m`, `y_m`, `time_s` — numeric vectors | Exact spatial coordinates and time samples. Document whether physical extents refer to cell centers or cell edges. |
   | `sourcePulse`, `timeWeight` — numeric vectors | Actual fine-grid pulse and objective weight K; pulse formula/identifier, tau, amplitudes and injection convention also live in metadata. Do not depend on a saved anonymous function handle. |
   | `epsr_background`, `cond_e_background`, `cond_m_background`, `murx_background`, `mury_background` — arrays | Complete background medium maps with their native grid dimensions; scalar descriptors and units in model. |
   | `doi_mask`, `valid_domain_mask`, `doi_linear_indices` — logical/integer arrays | DOI and valid physical-domain masks, excluding PML as appropriate, and the MATLAB linear-index mapping used by compact histories. |
   | `epsr_initial` — array | The actual initial estimate, defaulting to the background. Record which material properties are estimated and which are prescribed. Save any nonzero initial field arrays explicitly when supported. |

   **truth.mat — scene geometry and the true image**

   | Variable/type | Saved contents |
   | --- | --- |
   | `scene` — scalar struct | Scene/system IDs, fingerprints, target-generation settings and RNG state, overlap/placement policy, geometry rasterization method, 2x construction-grid description and sampling rule, and synthetic-data assumptions. |
   | `targets` — struct array, one entry per object | Shape name; center, dimensions, radius/bounds/vertices as applicable; orientation; geometry in index and physical coordinates; epsr, electrical/magnetic conductivity and relative permeability; object/class IDs and construction order. Record actual explicit geometry even when randomly generated. |
   | `epsr_true` — double `[Nx,Ny]` | The assembled true relative-permittivity map actually passed to the measurement FDTD solver, including background. This is the requested true image. It is authoritative over an idealized drawing of the targets. |
   | `cond_e_true`, `cond_m_true`, `murx_true`, `mury_true` — arrays | All remaining material maps actually passed to the solver, including correct staggered dimensions. |
   | `target_masks` — logical `[Nx,Ny,numObjects]` | Actual rasterized object membership masks aligned with the simulation grid, generated with the same geometry/sampling pipeline as the material maps. |
   | `target_union_mask`, `label_image` — logical/integer arrays | Union of target membership and object/material labels. Save overlap semantics and layering; the individual masks preserve overlapping memberships if that mode is later enabled. |

   **acquisition.mat — shared initial measurements**

   | Variable/type | Saved contents |
   | --- | --- |
   | `measurement` — scalar struct | System/scene IDs and acquisition fingerprint, original TX/RX axis IDs, tensor convention, source/time references, signal convention (current total Ez, not background-subtracted scattering), units/normalization convention, complete-transmitter flags, and acquisition status. |
   | `rx_full` — double `[N,N,Nt]` | Full acquired initial dataset, saved once per scene. Preserve raw amplitudes. Current data are noiseless; explicitly record that. |
   | `measurementTiming` — struct | Per-transmitter measurement-solver durations, acquisition wall time, timestamps, and actual solve count. |
   | Noise/preprocessing metadata, when applicable | If added later, save exact settings and random state, preserve clean RX separately, and save the actual noisy/processed full tensor used. Generate a shared noise realization before slicing to support paired comparisons. No noise is introduced by this plan. |

   **result.mat — one reconstructed configuration**

   | Variable/type | Saved contents |
   | --- | --- |
   | `caseInfo` — scalar struct | Stable IDs; relative references and hashes for system/truth/acquisition; schema version; original active/disabled indices; original-to-local mapping; TX order; first original transmitter; selected rotation, canonical gaps, orbit size, and selection RNG provenance. |
   | `parameters` — scalar struct | All resolved inversion settings: requested iterations, sensitivity coarsening, initialization, lower bound, finite-difference rule/scale, conjugate-gradient/restart/projection conventions, temporal-weight reference, channel policy, and any stopping rules. Record defaults currently hard-coded in runFbts as well as user options. |
   | `epsr_est` — double `[Nx,Ny]` | Final raw reconstructed material map, after the last completed update. No display normalization, cropping, transpose, or quantization. |
   | `epsr_history_doi` — double `[numDoiPixels,I+1]` | Initial estimate and every completed post-update DOI estimate, in the saved DOI index order. Outside-DOI values remain available from the prescribed background. This supports later image metrics versus iteration without storing those metrics. |
   | `history` — scalar struct of vectors/small arrays | Evaluated state indices, cost and per-TX cost, alpha, unconstrained alpha, Polak–Ribiere beta, restart flags/reasons, projection flags/cell counts, directional derivative, step_a, step_q, finite-difference h, sensitivity dt/Nt, raw/projected gradient norms, update norms, requested/completed iterations and stopping status. Preserve the raw objective, channel counts, and measurement energy needed for later normalization comparisons. `step_q` is an optimization coefficient, not the requested image-quality Q metric. |
   | `rx_model_history` — double `[M,M,Nt,I]` | Predicted RX at each evaluated pre-update state, with original TX/RX IDs from caseInfo. Together with shared measured RX and K this permits residual/convergence analysis without rerunning FDTD. Do not duplicate sliced measured RX in every case. |
   | `rx_sensitivity_coarse_history` — double `[M,M,NtCoarse,I]` | Actual coarse directional sensitivity traces used for the step calculation, plus coarse time/weight/interpolation references. Store committed slices if a run stops early. Fine-grid interpolated sensitivity can be regenerated from these. |
   | `gradient_raw_doi_last`, `gradient_projected_doi_last`, `direction_doi_last` — arrays | Last evaluated raw gradient, projected gradient, and search direction, explicitly labeled with their material-state index. |
   | `timing` — scalar struct | Raw per-transmitter/per-iteration durations for forward, adjoint, sensitivity baseline, and sensitivity perturbation; per-iteration runtimes; setup/cache-load/reconstruction elapsed times; shared acquisition reference and zero new acquisition count. Final write/plot/total durations are committed in status.mat after result writing completes. |

   Here `I` is the number of completed updates and `M` is the number of active antennas. If adaptive temporal grids are introduced later, sensitivity histories need per-grid blocks rather than an assumed fixed-size tensor; the current factor and grids are fixed per case.

   State alignment is mandatory. Define estimate states `e_0` (initial) through `e_I` (final). Loop iteration j evaluates cost, gradient, predicted RX and sensitivity at `e_(j-1)` and then produces `e_j`. Therefore `rx_model_history(:,:,:,j)` corresponds to `epsr_history_doi(:,j)`, while the final image is column `I+1`. The final saved pre-update RX/cost must not be described as RX/cost of the final image. An optional explicit final forward evaluation may be added later, with its own state label and timing; it is not necessary for computing image-quality metrics.

   **status, failure, checkpoints, and provenance**

   | File/type | Saved contents |
   | --- | --- |
   | `status.mat` — scalar status struct | Queued/running/succeeded/failed state, timestamps, completed iterations, committed result hash, final write/plot/total runtimes, cache status, and any plot-only failure. Writing this small record is excluded from the reported result-write duration. |
   | `failure.mat` — failure struct | Case identity, failed stage and iteration/TX if known, exception ID/message/stack, completed timing/history indices, last valid checkpoint reference, and recovery status. Never advertise a partially written result as successful. |
   | `checkpoint.mat` — scalar state struct plus arrays | At a completed iteration boundary: current estimate, previous projected gradient and direction, last projection flag, next iteration, relevant resolved options and fingerprints, RNG state, and committed history positions. The estimate alone is insufficient to resume the conjugate-gradient trajectory. |
   | `provenance/` — source snapshot plus metadata | Git revision and dirty status, hashes and copies of the relevant MATLAB/C/build/helper sources including relevant untracked files, MEX binary hash, compiler/build details when available, MATLAB/platform information, and relevant server CPU/RAM/worker settings. Mark unavailable build facts as unknown. Do not capture unrelated environment variables or credentials. |
   | Optional figures | Render from saved arrays after data commit. Store rendering settings if figures are exported. These are previews, not the analysis dataset. |

   Keep numeric output in double precision by default, masks logical, and labels/indices appropriate integers. Record units and dimension order. Do not save Q, DSC, SSIM, PCC, NRMSE, or the current relative_error_percent in the new batch dataset. The existing demo can retain its legacy result contract.

   Raw aligned truth/reconstruction arrays support intensity-based metrics. The saved object masks support DSC after choosing a reconstruction segmentation rule. DOI/domain masks, coordinate conventions, background maps and unnormalized values support explicit ROI, dynamic-range and normalization decisions for Q/SSIM/PCC/NRMSE. Q's exact definition, SSIM window/range, DSC threshold, and treatment of constant/empty regions belong to the later analysis protocol. Preserve 2-D image geometry for SSIM rather than flattening the DOI.

   Full space-time Ez/Hx/Hy movies are not part of the default dataset. They are much larger than the maps and receiver traces and are not needed for these image metrics. Add an explicit optional diagnostic mode for selected full field movies, all gradient maps, or intermediate coarse baseline/perturbation traces if that research need arises. The default saved items above are the complete agreed proposal, not a claim to retain every internal transient.

   Acceptance gate: save/load equality for maps, masks, IDs and traces; alignment of masks with actual rasterized truth; reconstruction of a runnable cfg via the loader; correct state association for all histories; no per-case full acquisition duplication; correct final iteration even after early failure; and portable dataset loading after moving its root folder.

7. **Stage 5 — make server execution durable and efficient**

   Save each result/checkpoint to a temporary file in the destination filesystem and commit it by rename after successful writing. Commit complete status only after the corresponding result is valid. Preserve the last valid checkpoint until its replacement is safely committed. Do not claim guarantees beyond what the server filesystem provides.

   Add explicit resume mode. The present output-directory helper rejects every nonempty directory, so resume needs a distinct path that validates the manifest, schema, source/config fingerprints and saved representative choices. Skip only verified successful matching cases. Resume matching checkpoints at iteration boundaries; redo an interrupted partial iteration. Never silently mix results made with different options or code identities.

   Update progress/status after every case, and checkpoint at a configurable completed-iteration interval, default one for long server runs. Keep one writer per case and acquisition. A single coordinator refreshes the manifest/CSV; alternatively rebuild summaries from per-case status files. A shared acquisition failure invalidates the dependent scene's cases, not unrelated scenes. A failed reconstruction or plot should not terminate later independent cases.

   Disable figures by default and separate numerical saving from plotting in saveFbtsResults.m. Currently plotting happens before MAT writing, so a plot failure can prevent saving an otherwise completed reconstruction. Make optional plot failures independently reportable.

   Add a dry-run summary of scene count, classes, active transmitter totals, expected solver calls, disk requirements and estimated peak memory. For fixed I, with one N-transmitter acquisition per scene, total solver calls are `N + 4*I*sum(M_case)`. Report fine and coarse calls separately when estimating time; raw call count is not a runtime estimate.

   First establish serial correctness and measure memory on the actual server. At the current 400x400 grid and approximately 600 time samples, one double-precision full Ez history alone is about 768 MB (0.715 GiB). Forward/adjoint fields, sliced copies and temporaries coexist, so this is not a process memory budget. Bound optional process-based workers by measured peak usage and available MATLAB resources; avoid nested parallelism and concurrent writes to shared MAT files. Independent cases are suitable worker jobs once acquisition is complete.

   Reuse immutable source/time/PML setup within cases where safe and stream histories to disk. Do not cache current-estimate forward/adjoint/sensitivity results across configurations. More invasive solver changes, shared first-iteration field caches, GPU work, or changes to the optimization method should remain separate follow-up work with their own validation.

   Acceptance gate: interrupt/resume a small run and compare its numerical trajectory to an uninterrupted run; reject stale/corrupt caches and checkpoints; preserve deterministic case choices under scheduling changes; verify no overlapping file writers; and confirm failure isolation and accurate timing accounting.

8. **Stage 6 — integrated validation and server pilot**

   Extend fbtsBatchTest.m and fbtsRefactorTest.m and add focused acquisition/catalog/storage/resume tests. Preserve standalone demo behavior while updating batch expectations from exhaustive output and figure-dependent saving to the new mode and schema. Existing tests that require CSV creation only at the end need updating for durable progress records.

   Run a small multi-scene end-to-end batch with both symmetry modes and compare each selected cached case against a fresh acquisition of that same representative. Then run a bounded full-resolution pilot to measure actual wall time, peak RAM, output size and the effect of the existing sensitivity factor. Check coarse-grid antenna collisions before expensive inversion work.

   Validate meaningful convergence separately from successful execution. Reduced arrays may be underdetermined, and FBTS may stall or reach a zero projected gradient; this does not invalidate measurement slicing. Record outcomes and preserve diagnostics. If a new convergence/stopping rule is introduced to treat stationary points cleanly, make it explicit and test it; do not silently change the algorithm while refactoring data flow.

   Document schema, commands, option defaults, dataset-loading examples, interpretation of rotation/reflection classes, self-channel behavior, and image/state ordering. Continue with the original current numerical method after the gates pass. Synthetic same-model data remain subject to the usual limitations of evaluating on the same forward model used for inversion; record that provenance when interpreting later metrics.

9. **Implementation touchpoints and evidence**

   | Existing file | Planned change |
   | --- | --- |
   | `fbts_batch.m` | Thin entry point over scene/catalog/shared-acquisition orchestration; independent cases; durable summaries/resume. |
   | `runFbts.m` | Optional supplied measurements, shared numerical setup, aligned history/checkpoint instrumentation; preserve iterative solves. |
   | `resolveFbtsOptions.m` | Resolve and validate algorithm controls and relevant explicit defaults; introduce a separate batch-options resolver. |
   | `build_cfg.m` | Parameterized builder wrapper with recorded generation settings and RNG provenance; correct receiver-policy comment. |
   | `saveFbtsResults.m`, `plot_fbts.m` | Data-first saving and optional independent plotting; new batch schema with demo compatibility. |
   | `createFbtsOutputDirectory.m` | Separate new-run allocation from validated resume behavior. |
   | `prepareCoarseSensitivityCfg.m` | Preserve its calculation; expose resolved metadata and preflight collision/time-grid validation. |
   | Tests and README | Equivalence, catalog counts, deterministic sampling, serialization, failure/resume, headless execution and revised usage. |

   Local implementation evidence: acquireFbtsMeasurements.m acquires the true-scene data; runFbts.m computes independent current-estimate forward/adjoint fields and both sensitivity solves. Results retain explicitly aligned pre-update signals and post-update images. solver.c and solver_tr.c inject samples only at listed coordinates; fdtd_mex.c reads receiver fields without receiver loading. buildCircularAntennaArrayIdx.m constructs nominal equal angles and rounds to grid cells. gridtmz.c sets the solver's electromagnetic constants internally. saveFbtsResults.m now saves numerical data before plotting.

   External format references: [MathWorks matfile documentation](https://www.mathworks.com/help/matlab/ref/matlab.io.matfile.html) describes partial access and structure-field limitations; [MAT-file versions](https://www.mathworks.com/help/matlab/import_export/mat-file-versions.html) describes v7.3 capabilities and storage tradeoffs. [Fhager et al., Antenna Modeling and Reconstruction Accuracy (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3625614/) discusses sensitivity of microwave reconstruction to physical antenna modeling.

   Verification status at planning time: source inspection and the independent combinatorial enumeration are complete. The isolated MATLAB feasibility probe passed outside the filesystem sandbox: all 15 nonempty four-antenna subsets matched fresh measurements exactly; eight three-iteration comparisons at sensitivity factors 1 and 4 matched every non-timing result; masked-full and reduced adjoints matched exactly. Production-interface regression checks and a bounded 400×400 pilot have subsequently passed; see VALIDATION.md. The MATLAB implementation has been updated; the C numerical solver is unchanged.
