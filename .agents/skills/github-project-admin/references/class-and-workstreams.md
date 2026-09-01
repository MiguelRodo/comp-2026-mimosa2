# Class and Workstream design

Use Class or Issue Type to describe **what kind of work item this is**. Use Workstream to describe **which stable functional lane of the Project it belongs to**. They should add different information.

These are starter conventions, not a universal schema. Preserve useful local vocabulary, and use only the values a Project actually needs.

## Class / Issue Type

A good default vocabulary is:

| Class | Preferred colour | Use for |
| --- | --- | --- |
| Epic | BLUE | A broad coordination outcome that stays useful while several distinct pieces of work are tracked separately. |
| Task | YELLOW | The normal concrete action or outcome when no more specific type adds value. |
| Deliverable | ORANGE | A bounded artefact or event such as a report, presentation, submission, assessment, release or handover. |
| Analysis | PURPLE | Work whose main outcome is an analytical result or inference from data, models or results. |
| Research | PINK | Investigation, design exploration or a spike whose purpose is to resolve uncertainty or choose a direction. |
| Enhancement | GREEN | A bounded improvement to an existing system, method, process or teaching material. |
| Bug | RED | A defect, regression or correctness failure. |
| Documentation | GRAY | Documentation-only maintenance or reference work rather than a substantive deliverable. |

Not every Project needs every value. Organisation-native Issue Types may use a smaller shared vocabulary; map the same semantics to the provider's exact names where useful.

### Epic is not the default parent type

Parenthood and Class are independent. A Task, Deliverable, Analysis, Research item or other type may have sub-issues and remain that type.

Use Epic when the grouping itself has durable coordination value, usually because the outcome spans several independently actionable pieces, phases or functional lanes. Do **not** create an Epic merely because:

- an item is top-level;
- an item has one or two small children;
- several related issues exist;
- the Project needs a single root node.

A Project may have several top-level issues and no root Epic. A presentation with preparation sub-issues can remain a Deliverable. A data analysis with preprocessing and figure sub-issues can remain Analysis. Promote to Epic only when the broader outcome is genuinely useful as a planning object in its own right.

Likewise, avoid turning subject matter into Class values. Values such as `Raw data`, `Processed data`, `Baseline`, `Robustness`, `D2 presentation` or a particular method usually describe a workstream, phase, scenario or deliverable rather than a reusable work-item type.

## Workstream

A Workstream should be a relatively stable functional lane that remains meaningful across many issues. Prefer a concise set, commonly around four to eight values for an ordinary Project. Larger umbrella Projects may reasonably need more.

Do not use Workstream to duplicate:

- Class or Issue Type;
- Priority or Status;
- a one-off milestone or deadline;
- a routing or `subproject:*` label;
- a scenario name that is better expressed in the issue title or hierarchy.

If an issue spans several lanes, choose the primary one. If the distinction matters enough to manage separately, split the work with parent/sub-issues rather than turning Workstream into a multi-purpose taxonomy.

### Research and scientific Projects

Start with the lanes that are actually needed:

| Workstream | Preferred colour |
| --- | --- |
| Study design | BLUE |
| Data | PURPLE |
| Methods | PINK |
| Implementation | YELLOW |
| Analysis | GREEN |
| Validation | RED |
| Reporting | ORANGE |
| Administration | GRAY |

Useful optional additions include `Simulation` (ORANGE) for simulation-heavy studies and `Collaboration` (PINK) where coordination itself creates substantive work.

Common consolidations are:

- `Diagnostics` into `Validation`;
- method-specific names such as `Covariates` into `Methods`;
- `Real data` into `Analysis` or `Data` depending on the work;
- `Baseline` and `Robustness` into `Simulation`, `Validation` or `Analysis`;
- milestone-specific streams such as `D1 Proposal`, `D2 Presentation`, `D3 Write-up` and `D4 Poster` into `Reporting`, with the issues themselves typed as Deliverables where appropriate.

### Teaching Projects

| Workstream | Preferred colour |
| --- | --- |
| Delivery | BLUE |
| Assessment | RED |
| Materials | GREEN |
| Student support | PURPLE |
| Administration | GRAY |

Add `Supervision` (PINK) when supervision is genuinely part of the Project rather than a separate sub-project. An umbrella teaching Project that also contains research supervision can use the relevant research lanes alongside these rather than forcing all work into one profile.

### Personal and household Projects

| Workstream | Preferred colour |
| --- | --- |
| Finance | GREEN |
| Administration | GRAY |
| Communication | BLUE |
| Digital | PURPLE |
| Maintenance | ORANGE |

Use routing or sub-project labels for areas such as household or family when those labels already define the grouping; do not repeat them as Workstreams unless they genuinely describe a different functional lane.

### Software and product Projects

| Workstream | Preferred colour |
| --- | --- |
| Design | BLUE |
| Implementation | YELLOW |
| Integrations | PURPLE |
| Testing | RED |
| Documentation | GREEN |
| Release / operations | ORANGE |
| Governance | GRAY |
| User experience | PINK |

Component-specific streams such as a particular CLI, API, UI or provider are still reasonable when the component is a long-lived lane with enough work to justify it. Prefer the standard names when they carry the same meaning; use local names when they convey real structure rather than mere novelty.

## Colours are presentational

Preferred colours make repeated names easier to recognise across Projects, but colour is not semantic state. Preserve useful existing colours unless the requested outcome includes changing them.

GitHub's standard single-select palette is `BLUE`, `GRAY`, `GREEN`, `ORANGE`, `PINK`, `PURPLE`, `RED` and `YELLOW`. If there are more categories than distinct colours, reuse colours. If another provider exposes additional colours, using them is fine. Lack of a unique colour must never block classification, validation or ordinary Project administration unless a repository has explicitly declared an exact palette as a local constraint.

## Migration and refinement

When organising an existing Project, inspect the actual issues before proposing changes. Prefer consolidation over churn:

1. preserve values that already carry a useful stable distinction;
2. identify duplicates between Class and Workstream;
3. fold milestone-, scenario- or method-specific Workstreams into stable lanes where that improves navigation;
4. avoid mass reclassification merely to match these examples;
5. show the proposed vocabulary and issue moves before changing live fields when the request is a broad reorganisation.
