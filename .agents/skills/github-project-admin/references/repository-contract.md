# Repository contract

Repository-specific GitHub Project configuration lives under `.projects/`. It contains only facts that differ between repositories. Common operating behaviour belongs in the skill.

## Source precedence

Use a local replacement skill only when this exact file exists:

```text
.projects/skills/github-project-admin/SKILL.md
```

The `.projects/` directory by itself never overrides the canonical skill.

Always read `.projects/project.md` after selecting the skill.

## Single-Project form

Use this form when one repository resolves to one Project:

```markdown
# GitHub Project configuration

| Key | Value |
| --- | --- |
| Contract version | 1 |
| Mode | single |
| Issue repository | octo-org/example |
| Project owner | octo-org |
| Project number | 12 |
| Project title | Example planning |
| Routing | linked repository |
| Privacy | repository |
| Issue write-up style | tidy |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- | --- |
| Class | organization issue type | Issue Type |
| Priority | organization issue field | Priority |
| Status | project field | Status |
| Workstream | project field | Workstream |
| Due date | project field | Target date |

## Priority mapping

| Common value | Provider value |
| --- | --- |
| P0 | Urgent |
| P1 | High |
| P2 | Medium |
| P3 | Low |

## Governance

- Exact user-requested changes may be applied without retrieving scope-design sources.
- Keep private material out of this repository.
```

Setup discovers whether `Project owner` is a user or organisation from GitHub. An optional `Owner type` row may assert `user` or GitHub's provider spelling `organization`; setup fails if that assertion disagrees with the live owner. `Routing` may name a linked repository, one exact routing label, or another deterministic repository-specific rule.

## Multi-Project form

Use `.projects/project.md` as a dispatcher:

```markdown
# GitHub Project dispatcher

| Key | Value |
| --- | --- |
| Contract version | 1 |
| Mode | dispatcher |
| Issue repository | octo-user/issues |
| Privacy | private repository |
| Governance | collaborative |

## Routes

| Project key | Routing label | Project number | Contract |
| --- | --- | --- | --- |
| alpha | project:alpha | 4 | .projects/projects/alpha.md |
| beta | project:beta | 5 | .projects/projects/beta.md |
```

Each referenced file uses the single-Project form with `Mode` set to `project` and adds a `Project key` metadata row. Its key, `label:` routing value, Project number and issue repository must match the dispatcher row exactly. Route keys, routing labels and Project numbers must each be unique. A supplied label, key and number must resolve to the same row.

The guided initializer may create this dispatcher with only the route-table header. That zero-route form is a valid saved onboarding state, but it cannot resolve ordinary administration. Rerun the initializer to add one Project at a time. Each addition discovers the live Project, writes one child contract, updates the dispatcher and validates the combined result before preserving it. Onboarding records routing labels in the contracts but does not create or apply them on GitHub.

## Issue write-up style

A resolved Project contract may contain an `Issue write-up style` metadata row. This controls how much an agent expands issue titles and bodies when creating an issue or substantially rewriting one. It does not override a more recent explicit user instruction.

Supported values are:

- `unrestricted`: the agent may add useful grounded structure, context, implementation detail, acceptance criteria or decomposition when helpful;
- `tidy`: light expansion and tidying only, and the default when the row is absent;
- `minimal`: minimal and direct wording, with no inferred context, rationale, implementation detail, acceptance criteria or extra subtasks.

For a multi-Project repository, put the setting in the resolved `.projects/projects/*.md` child contract so different Projects can use different defaults. Users may edit this row directly when they want a different style.

## Class and Workstream vocabularies

Class or Issue Type describes the kind of work item. Workstream describes the stable functional lane. Follow [the Class and Workstream design reference](class-and-workstreams.md) when proposing or refining these values.

The shared profiles are defaults and examples, not implicit contract values. A repository may keep a smaller or more domain-specific vocabulary. Existing useful local values should not be renamed merely to match the examples.

Important defaults are:

- `Task` is the ordinary fallback Class.
- `Epic` is reserved for a broad coordination outcome. Top-level placement or having children does not by itself make an item an Epic.
- A Task, Deliverable, Analysis or other type may have sub-issues and remain that type.
- Workstreams should be stable lanes and should not duplicate Class, Priority, Status, one-off milestones, routing labels or sub-project labels.
- A contract may list exact Class or Workstream values when the Project genuinely requires them; otherwise the agent should inspect the issue set and propose a concise vocabulary before live changes.

## Field locations and mappings

For each dimension, record the semantic name, provider location and exact provider field name. Typical locations are:

- repository issue metadata;
- organisation Issue Type;
- organisation issue field;
- Project field;
- repository label;
- native parent/sub-issue relationship.

Do not store transient GraphQL node IDs, REST option IDs or credentials. Discover IDs and live options at operation time.

The completed Priority table must contain P0, P1, P2 and P3 exactly once, with four distinct, non-empty provider values. Omit no value. When the provider uses `Urgent`, `High`, `Medium` and `Low`, use the default table. A repository may use an exact one-to-one override such as P0, P1, P2 and P3.

The guided initializer does not change live Priority options or ask a non-technical operator to interpret them. Until an agent has inspected the live field, the initial contract may use this exact section instead:

```markdown
## Priority mapping

Priority mapping status: pending
```

This is a safe incomplete state, not a default mapping. Its Field locations row may use `pending live inspection` until the provider location is confirmed. The repository may use other configured dimensions, but an agent must not rank, read semantically or change Priority until it records that location and replaces the marker with a complete one-to-one table. Adding, removing or renaming a provider option remains a separate live mutation that requires explicit authority.

## Option colours

A repository may make a single-select palette exact by using an `Option` and `Colour` table in the field's values section:

```markdown
## Class values

| Option | Colour |
| --- | --- |
| Epic | BLUE |
| Task | YELLOW |
| Bug | RED |
| Documentation | GREEN |
```

Supported GitHub colours are `BLUE`, `GRAY`, `GREEN`, `ORANGE`, `PINK`, `PURPLE`, `RED` and `YELLOW`. Reusing a colour is allowed.

If a contract lists values without colours, colour is not a contract constraint. When creating or organising a Project, an agent may choose stable colours without asking if the choice is purely presentational. The preferred colours in [the Class and Workstream design reference](class-and-workstreams.md) are reusable defaults, not semantic state.

If there are more categories than distinct provider colours, reuse colours. If another provider exposes additional colours, those may be used. A lack of unique colours must not block ordinary administration or classification. Only an explicitly declared exact palette is a local contract constraint. Preserve useful existing colours unless the requested outcome includes changing them.

## Governance and source rules

Record only local constraints, for example:

- whether issues may contain private material;
- whether the repository is personal, shared or public;
- whether assignment defaults exist;
- whether a source must be consulted before inventing or restructuring scope;
- whether routing labels or sub-project labels are required;
- which external mirror is read-only.

Do not repeat fresh inspection, narrow writes, stale refusal, preservation or readback rules. The skill already owns them.

## Exceptional setup

Use `.projects/setup.sh` only for prerequisites unique to this repository. The shared `scripts/setup.sh` discovers it automatically from the repository root.

By default the local script extends the shared setup and runs after the common GitHub checks. It must not call or copy the shared setup.

To replace the common setup entirely, put this exact marker within the first 20 lines:

```bash
# github-project-admin: override
```

In override mode the shared entry point disables shell tracing, finds the repository and immediately runs `.projects/setup.sh`; it does not install `gh`, check authentication or validate the contract. The local script receives `PROJECTS_REPOSITORY_ROOT` and `PROJECTS_SETUP_MODE` in its environment.

Keep local setup idempotent so rerunning it after a partial failure is safe. Skill installation and updates must never edit or delete `.projects/setup.sh`. Repository language runtimes, such as R, remain separate from GitHub Project administration unless a real Project operation depends on them.
