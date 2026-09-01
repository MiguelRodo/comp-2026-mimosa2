---
description: Administer GitHub issues and Projects from short outcome requests. Use for Project-aware inspection, prioritisation, creation, updates, assignment, routing, hierarchy, field changes, or when a surface must return minimal executable gh commands instead of writing directly.
metadata:
    github-path: skills/github-project-admin
    github-ref: refs/heads/main
    github-repo: https://github.com/MiguelRodo/projects
    github-tree-sha: fab254d77687b7e117e70d5373e42c41ef77c40b
name: github-project-admin
---
# GitHub Project administration

Treat the user's request as the desired outcome. Do not require the user to repeat the operating procedure in this skill.

## Select the skill and local contract

1. Find the target repository root.
2. If `.projects/skills/github-project-admin/SKILL.md` exists there and is not this skill, load that deliberate local replacement and stop applying this copy. Do not merge two skill definitions implicitly.
3. Otherwise use this canonical skill.
4. Read `.projects/project.md`. If it is missing, stop and identify the missing repository contract. Do not infer topology from a repository name, issue title, Project title, or previous run.
5. In a shell-capable environment, run `scripts/validate-contract.sh REPOSITORY_ROOT` before relying on the contract.

Read [the repository contract reference](references/repository-contract.md) when creating, migrating, validating, or interpreting `.projects/` files.

## Resolve exactly one context

- A `single` contract resolves directly.
- A `dispatcher` contract must resolve exactly one row by an explicit Project number, routing label, Project key, or other exact identifier supplied by the request or live issue state. Then read the referenced `.projects/projects/*.md` contract.
- All identifiers supplied by the user, contract, issue and Project must agree. Stop on zero matches, multiple matches, or disagreement.
- Treat live GitHub as authority for current issue, membership, field, option and hierarchy state. Treat `.projects/` as authority for topology, mappings, governance and source requirements.
- Retrieve an external source only when the local contract requires it for this operation. An exact requested change is not a scope-design task.

## Interpret short requests

A request such as `Set example#313 to P2.` supplies an exact target and desired value. It authorises only that delta. A request such as `What are the highest-priority open items in example?` is read-only.

Ask a question only when a missing fact would change the target or outcome. Do not ask the user to restate fresh inspection, narrow mutation, preservation, stale checks or readback requirements.

Use this default common Priority vocabulary unless the resolved contract declares a complete override:

| Common value | Provider value |
| --- | --- |
| P0 | Urgent |
| P1 | High |
| P2 | Medium |
| P3 | Low |

Mappings must be one-to-one. Never collapse two provider values into one common value or use different read and write mappings. Discover the exact live provider option before a write and stop if it does not match the contract.

For single-select option colours, follow any exact palette in the repository contract. If none is declared and the choice is only presentational, choose stable colours that fit the option meanings without asking the user to design a palette. Preserve useful existing colours. Read [the repository contract reference](references/repository-contract.md) when creating or changing option definitions.

## Choose an execution surface

Use the first capable surface:

1. a future `projects` CLI, or `projectscli` if that is the installed name;
2. proven native `gh` commands, versioned GitHub REST, and GraphQL;
3. an equivalent authenticated provider connector that can perform the same inspection and independent readback.

The future CLI is optional and does not exist merely because this skill mentions it. Until it is available, use the direct operations in [the GitHub operations reference](references/github-operations.md).

If the current surface cannot perform an authorised mutation, inspect as far as safely possible and return the smallest executable command block that completes the operation. Use placeholders only for facts that cannot be discovered. Do not claim that returned commands ran.

Run `scripts/setup.sh` when preparing an environment or when `gh` prerequisites are missing. The host must provide credentials and network access. Never print, persist, transform or request a token in a prompt.

When adopting the skill in a repository that does not yet have `.projects/project.md`, run `scripts/init-project.sh` from that repository. It discovers live GitHub facts, asks only for local choices, writes the first single-Project contract and adds a bounded `AGENTS.md` routing section. It never mutates live issues or Projects. For a multi-Project answer, it stops with an agent handoff instead of inventing routing.

## Inspect and plan

For every target issue, read the current:

- stable issue identity, URL and state;
- Project identity and membership;
- affected field definition, location, option and current value;
- labels, assignees, milestone and relevant native hierarchy;
- any other state the proposed operation could replace or remove.

State the exact requested delta. Separate it from necessary implied changes. Adding Project membership is an implied change only when the requested Project field cannot exist without it and the contract authorises membership.

For a read-only ranking, translate provider values through the resolved mapping and order P0, P1, P2, then P3. Apply any requested tie-breaker; otherwise report the tie rather than inventing precision.

## Apply a mutation safely

Before writing:

1. discover current IDs and options rather than reusing IDs from documentation or a previous run;
2. re-read the stale-sensitive target immediately before the write;
3. compare it with the inspected state and stop if the target, membership, affected value or preservation set changed;
4. use the narrowest endpoint or command that changes only the authorised field;
5. preserve every unrequested field, label, assignee, milestone, relationship, body, comment and Project membership.

Do not replace a whole collection when an additive or field-specific operation exists. Do not remove Project membership without explicit authority because doing so also removes Project-local values. Do not retry an uncertain create or mutation blindly.

After writing, perform a separate targeted read. A successful mutation response is not verification. Report success only when the independently observed value equals the requested provider-native value and the preservation set remains intact.

Stop on ambiguity, missing permission, stale state, unavailable required sources, an unexpected collateral change, or failed readback.

## Report

For a read, report the resolved context, relevant values, ranking rule and any unavailable field.

For a write, report separately:

- requested delta;
- applied delta;
- independent readback;
- preserved state relevant to the mutation;
- any remaining manual action.

Keep private issue content within its authorised repository and surface. Never copy credentials or private task details into a public issue, PR, log or evidence report.

When authoring ChatGPT Project or another provider's standing instructions, follow [the provider instruction reference](references/provider-project-instructions.md). Keep that layer as a pointer, not a second copy of this procedure.
