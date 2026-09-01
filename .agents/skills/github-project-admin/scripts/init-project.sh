#!/usr/bin/env bash

# Guided, repository-local onboarding for github-project-admin. This script
# writes repository configuration only. It never changes live GitHub issues or
# Project fields. Keep this file compatible with Bash 3.2.
set +x
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
skill_dir="$(cd "$script_dir/.." && pwd)"
generated_contract=""
section_number=0
colour_reset=""
colour_bold=""
colour_blue=""
colour_green=""
colour_yellow=""
colour_red=""

if [[ -t 1 && -z "${NO_COLOR:-}" && "${TERM:-dumb}" != "dumb" ]]; then
  colour_reset="$(printf '\033[0m')"
  colour_bold="$(printf '\033[1m')"
  colour_blue="$(printf '\033[34m')"
  colour_green="$(printf '\033[32m')"
  colour_yellow="$(printf '\033[33m')"
  colour_red="$(printf '\033[31m')"
fi

section() {
  section_number=$((section_number + 1))
  printf '\n%s%s%d. %s%s\n' \
    "$colour_bold" "$colour_blue" "$section_number" "$1" "$colour_reset"
  printf '%s\n' '----------------------------------------'
}

success() {
  printf '%s%s[OK]%s %s\n' "$colour_bold" "$colour_green" "$colour_reset" "$1"
}

note() {
  printf '%s[INFO]%s %s\n' "$colour_blue" "$colour_reset" "$1"
}

warning() {
  printf '%s%s[WARNING]%s %s\n' \
    "$colour_bold" "$colour_yellow" "$colour_reset" "$1" >&2
}

die() {
  printf '%s%s[ERROR]%s %s\n' \
    "$colour_bold" "$colour_red" "$colour_reset" "$1" >&2
  exit 1
}

cleanup() {
  if [[ -n "$generated_contract" && -f "$generated_contract" ]]; then
    rm -f -- "$generated_contract"
  fi
}
trap cleanup EXIT

ask_default() {
  local prompt="$1" default="$2" answer
  read -r -p "$prompt [$default]: " answer ||
    die "input ended before setup was complete"
  printf '%s' "${answer:-$default}"
}

ask_yes_no() {
  local prompt="$1" default="$2" answer hint
  if [[ "$default" == "yes" ]]; then
    hint="Y/n"
  else
    hint="y/N"
  fi
  while true; do
    read -r -p "$prompt [$hint]: " answer ||
      die "input ended before setup was complete"
    answer="${answer:-$default}"
    answer="$(printf '%s' "$answer" | tr '[:upper:]' '[:lower:]')"
    case "$answer" in
      y|yes) return 0 ;;
      n|no) return 1 ;;
      *) echo "Please answer yes or no." >&2 ;;
    esac
  done
}

ask_choice_default() {
  local prompt="$1" minimum="$2" maximum="$3" default="$4" answer
  while true; do
    answer="$(ask_default "$prompt" "$default")"
    if [[ "$answer" =~ ^[0-9]+$ ]] &&
       ((answer >= minimum && answer <= maximum)); then
      printf '%s' "$answer"
      return 0
    fi
    echo "Choose a number from $minimum to $maximum." >&2
  done
}

require_text() {
  local prompt="$1" answer
  while true; do
    read -r -p "$prompt: " answer ||
      die "input ended before setup was complete"
    if [[ -n "$answer" && "$answer" != *"|"* && "$answer" != *$'\n'* ]]; then
      printf '%s' "$answer"
      return 0
    fi
    echo "Enter a value without a table separator (|)." >&2
  done
}

append_agents_pointer() {
  local agents_file="$repository_root/AGENTS.md"
  local marker='<!-- github-project-admin:start -->'
  if [[ -f "$agents_file" ]] && grep -Fq "$marker" "$agents_file"; then
    note "AGENTS.md already contains the GitHub Project starting point."
    return 0
  fi

  if [[ -s "$agents_file" ]]; then
    printf '\n' >>"$agents_file"
  fi
  cat >>"$agents_file" <<'EOF'
<!-- github-project-admin:start -->
## GitHub issues and Projects

For GitHub issue or Project administration, use
`.agents/skills/github-project-admin/SKILL.md` and read
`.projects/project.md` before acting.
<!-- github-project-admin:end -->
EOF
  success "Added the GitHub Project starting point to AGENTS.md."
}

print_auth_help() {
  cat <<'EOF'

Run these commands, then start the initializer again:

  gh auth login --web --scopes "project,read:org"
  gh auth status

If you were already signed in but Project access is missing, run:

  gh auth refresh --scopes "project,read:org"
EOF
}

relative_skill_path() {
  case "$skill_dir" in
    "$repository_root"/*)
      printf '%s' "${skill_dir#"$repository_root"/}"
      ;;
    *)
      return 1
      ;;
  esac
}

save_onboarding_files() {
  local skill_path="" branch=""
  local -a save_paths
  save_paths=()

  skill_path="$(relative_skill_path 2>/dev/null || true)"
  if [[ -n "$skill_path" && -d "$skill_path" ]]; then
    save_paths+=("$skill_path")
  fi
  if [[ -e "$repository_root/.projects/project.md" ]]; then
    save_paths+=(".projects/project.md")
  fi
  if [[ -e "$repository_root/AGENTS.md" ]]; then
    save_paths+=("AGENTS.md")
  fi

  section "Save the repository setup"
  if ((${#save_paths[@]} == 0)) ||
     [[ -z "$(git status --porcelain -- "${save_paths[@]}")" ]]; then
    success "The onboarding files are already committed."
    return 0
  fi

  echo "The setup changed only these onboarding paths:"
  printf '  %s\n' "${save_paths[@]}"
  echo
  if ! ask_yes_no "May I stage, commit and push these onboarding files?" no; then
    note "The files were left uncommitted."
    echo "Review them with:"
    echo
    printf '  git status --short --'
    printf ' %q' "${save_paths[@]}"
    printf '\n'
    return 0
  fi

  if ! git add -- "${save_paths[@]}"; then
    warning "Git could not stage the onboarding files. Nothing was committed or pushed."
    echo "After fixing the reported problem, run:"
    echo
    printf '  git add --'
    printf ' %q' "${save_paths[@]}"
    printf '\n'
    return 0
  fi

  if git diff --cached --quiet -- "${save_paths[@]}"; then
    success "There was nothing new to commit."
    return 0
  fi

  if ! git commit -m "Configure GitHub Project administration" -- \
       "${save_paths[@]}"; then
    warning "Git could not create the commit. The onboarding files remain staged."
    echo "Fix the error above, then run:"
    echo
    echo "  git commit -m \"Configure GitHub Project administration\""
    echo "  git push"
    return 0
  fi
  success "Created the onboarding commit."

  branch="$(git symbolic-ref --quiet --short HEAD 2>/dev/null || true)"
  if [[ -z "$branch" ]]; then
    warning "The repository is in detached-HEAD state, so the commit was not pushed."
    echo "Create or switch to a branch, then run git push."
    return 0
  fi

  if git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}' \
       >/dev/null 2>&1; then
    if git push; then
      success "Pushed the onboarding commit."
    else
      warning "The commit is safe locally, but Git could not push it."
      echo "After fixing access, the remote or branch protection, run:"
      echo
      echo "  git push"
    fi
  else
    if git push -u origin "$branch"; then
      success "Pushed the onboarding commit."
    else
      warning "The commit is safe locally, but Git could not push it."
      echo "After fixing access, the remote or branch protection, run:"
      echo
      printf '  git push -u origin %s\n' "$branch"
    fi
  fi
}

print_chatgpt_setup() {
  section "Use the repository with ChatGPT"
  cat <<EOF
1. Open https://chatgpt.com/projects and create or open a ChatGPT Project.
2. Make $repository available to that Project through its GitHub connection.
3. Paste this into the ChatGPT Project instructions:

  For work concerning a GitHub repository, especially reading or updating
  GitHub issues or Projects, first retrieve and follow the target repository's
  AGENTS.md. Follow the skill and configuration files it references. If the
  repository or AGENTS.md is unavailable, say so rather than guessing.

  Treat my prompt as the desired outcome. If this chat cannot make a required
  GitHub change, return the smallest safe command block for me to paste into a
  terminal, including a check of the result.

After that, ask for the outcome you want in ordinary language. ChatGPT will do
what its connection allows and return terminal commands for anything it cannot
do directly.
EOF
}

print_codex_setup() {
  section "Use the repository with Codex cloud"
  cat <<EOF
1. Open https://chatgpt.com/codex/settings/environments and create an environment.
2. Choose the $repository repository.
3. Use this setup command:

  bash .agents/skills/github-project-admin/scripts/setup.sh

4. Create a classic GitHub personal access token at:

  https://github.com/settings/tokens/new

   Give it an expiry and the repo, read:org and project scopes. If your
   organisation uses SSO, authorise the token for that organisation.
5. Add the token to the environment as an environment variable named GH_TOKEN.
   Do not add it as a setup-only secret because the agent needs it after setup.
6. Enable internet access for the agent phase and allow:

  github.com
  api.github.com

Codex can then read AGENTS.md, run the GitHub commands and verify the result.
EOF
}

print_single_project_first_request() {
  local class_name="Class" ending choice
  if [[ "$project_owner_type" == "organization" ]]; then
    class_name="Issue Type"
  fi

  section "Choose a useful first request"
  if ! ask_yes_no \
    "Would you like an agent to organise existing issues with $class_name and Workstream next?" \
    yes; then
    note "Setup is complete. You can now make ordinary requests when you need them."
    return 0
  fi

  echo "How should the agent proceed?"
  echo "  1. Show you its plan before changing GitHub"
  echo "  2. Carry out the work and verify it"
  choice="$(ask_choice_default "Choose 1 or 2" 1 2 1)"
  if [[ "$choice" == "1" ]]; then
    ending="Give me an overview of what you plan to do based on this request. Do not make changes until I approve the plan."
  else
    ending="Do this now, then independently verify and summarise the changes."
  fi

  cat <<EOF

Use this as your first request in ChatGPT or Codex:

  Start from AGENTS.md. Inspect the current issues and GitHub Project. Confirm
  the pending Priority location and mapping from the Project's existing field
  without adding, removing or renaming options. Set up or refine $class_name and
  Workstream with sensible values and colours based on the existing issues,
  then organise those issues using the Project fields and useful native
  parent/sub-issue relationships. $ending
EOF
}

print_multi_project_first_request() {
  section "Finish the multi-Project routing"
  cat <<EOF
This repository uses several Projects, so an agent must help define the routing.
Use this as the first request in ChatGPT or Codex:

  Start from AGENTS.md. Finish setting up $repository as a $governance
  multi-Project repository. Inspect what GitHub and the repository already
  provide. Ask me only for Project owners, Project numbers and routing choices
  that cannot be discovered safely. Create and validate the dispatcher and
  Project contracts. Give me an overview of what you plan to do before changing
  live issues or Projects.

After that configuration is committed, you can ask the agent to organise the
existing issues and set up Issue Type or Class, Workstream and useful native
parent/sub-issue relationships across the resolved Projects.

To add another Project later, ask:

  Start from AGENTS.md. Add Project <number> owned by <user or organisation> to
  this repository's existing multi-Project configuration. Preserve current
  routes, ask me only for the new routing decision, and validate the result.
EOF
}

command -v git >/dev/null 2>&1 || die "git is required"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 ||
  die "run this command inside the repository you want to configure"
repository_root="$(git rev-parse --show-toplevel)"
cd "$repository_root"

section "Check GitHub"
if ! command -v gh >/dev/null 2>&1; then
  echo "Install GitHub CLI from https://cli.github.com/ and run:" >&2
  print_auth_help >&2
  die "GitHub CLI (gh) is not installed"
fi
if ! gh auth status >/dev/null 2>&1; then
  print_auth_help >&2
  die "GitHub CLI is not authenticated"
fi
authenticated_login="$(gh api user --jq .login 2>/dev/null)" ||
  die "could not read the active GitHub login"
if ! gh project list --owner "$authenticated_login" --limit 1 \
     >/dev/null 2>&1; then
  print_auth_help >&2
  die "the active GitHub login does not have Project access"
fi
success "GitHub CLI is authenticated as $authenticated_login with Project access."

repository_record="$(gh repo view --json nameWithOwner,visibility \
  --jq '[.nameWithOwner,.visibility] | @tsv')" ||
  die "could not discover the current GitHub repository"
IFS=$'\t' read -r repository visibility <<<"$repository_record"
[[ "$repository" =~ ^[^/[:space:]]+/[^/[:space:]]+$ ]] ||
  die "GitHub returned an invalid repository identity"
repository_owner="${repository%%/*}"

contract_file="$repository_root/.projects/project.md"
if [[ -e "$contract_file" ]]; then
  section "Check the existing repository setup"
  note "A repository contract already exists at .projects/project.md."
  bash "$script_dir/validate-contract.sh" "$repository_root"
  append_agents_pointer
  success "The existing repository setup was not replaced."
  save_onboarding_files
  echo
  note "For usage and update instructions, read .agents/skills/github-project-admin/README.md."
  exit 0
fi

section "Configure the repository"
cat <<EOF
I will ask a few questions to configure $repository so that ChatGPT and coding
agents can understand its GitHub Project and work with it safely. GitHub will
supply facts it already knows; the questions are only about local choices.

No live issue or Project value will be changed during this setup.
EOF

if ask_yes_no \
  "Will anyone else work with you on Projects managed from this repository?" \
  no; then
  governance="collaborative"
else
  governance="personal"
fi

if ! ask_yes_no "Does this repository use one GitHub Project?" yes; then
  append_agents_pointer
  note "Several Projects need an explicit routing rule, so no contract was guessed."
  save_onboarding_files
  print_chatgpt_setup
  print_codex_setup
  print_multi_project_first_request
  echo
  success "The repository is ready for the multi-Project configuration handoff."
  echo "The initializer did not change any live GitHub issue or Project value."
  exit 0
fi

project_owner="$(ask_default \
  "GitHub user or organisation that owns the Project" "$repository_owner")"
[[ "$project_owner" =~ ^[A-Za-z0-9_.-]+$ ]] ||
  die "invalid Project owner login"

echo
echo "Find the Project number after /projects/ in its web address:"
echo "  https://github.com/users/example/projects/40       means 40"
echo "  https://github.com/orgs/example-org/projects/12    means 12"
project_number="$(require_text "Project number")"
[[ "$project_number" =~ ^[1-9][0-9]*$ ]] ||
  die "Project number must be a positive integer"

observed_owner_type="$(gh api "users/$project_owner" --jq .type 2>/dev/null)" ||
  die "could not discover whether the Project owner is a person or organisation"
case "$observed_owner_type" in
  User)
    project_owner_type="user"
    project_selector='.data.user.projectV2'
    project_query='query($login: String!, $number: Int!) { user(login: $login) { projectV2(number: $number) { number title public } } }'
    class_location="project field"
    class_field="Class"
    ;;
  Organization)
    project_owner_type="organization"
    project_selector='.data.organization.projectV2'
    project_query='query($login: String!, $number: Int!) { organization(login: $login) { projectV2(number: $number) { number title public } } }'
    class_location="organization issue type"
    class_field="Issue Type"
    ;;
  *) die "unsupported GitHub owner type: $observed_owner_type" ;;
esac

project_record="$(gh api graphql -f query="$project_query" \
  -F login="$project_owner" -F number="$project_number" \
  --jq "$project_selector | [.number,.title,(.public|tostring)] | @tsv")" ||
  die "could not read Project $project_owner/$project_number"
IFS=$'\t' read -r observed_number project_title project_public <<<"$project_record"
[[ "$observed_number" == "$project_number" && -n "$project_title" ]] ||
  die "GitHub did not return the expected Project"
[[ "$project_title" != *"|"* && "$project_title" != *$'\n'* ]] ||
  die "this Project title needs agent-assisted contract generation"

visibility_lower="$(printf '%s' "$visibility" | tr '[:upper:]' '[:lower:]')"
if [[ "$visibility_lower" == "private" ]]; then
  privacy="private repository"
elif [[ "$project_public" == "false" ]]; then
  privacy="$visibility_lower repository with a private Project"
else
  privacy="$visibility_lower repository"
fi

success "Found $project_title, owned by the GitHub $project_owner_type $project_owner."

generated_contract="$(mktemp)"
cat >"$generated_contract" <<EOF
# GitHub Project configuration

| Key | Value |
| --- | --- |
| Contract version | 1 |
| Mode | single |
| Issue repository | $repository |
| Project owner | $project_owner |
| Project number | $project_number |
| Project title | $project_title |
| Routing | Project membership; no routing label |
| Privacy | $privacy |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- | --- |
| Class | $class_location | $class_field |
| Priority | pending live inspection | Priority |
| Status | project field | Status |
| Workstream | project field | Workstream |
| Due date | project field | Target date |
| Parent | native issue relationship | Parent issue |

## Priority mapping

Priority mapping status: pending

The initializer left the Project's existing Priority field and options unchanged.
Before using Priority, an agent must confirm its provider location, inspect the
live options and replace the pending status with a complete one-to-one P0, P1,
P2 and P3 mapping.

## Class and Workstream

Class and Workstream option sets are intentionally not fixed by onboarding.
An agent may inspect the existing issues and suggest a concise, useful vocabulary
before the live Project is changed.

## Governance

- This is a $governance Project.
- Project membership determines Project scope; no routing label is required.
- Labels must not duplicate Class, Priority, Status or Workstream.
- Assignment is explicit only unless a later repository decision says otherwise.
- Exact requested administration requires no separate scope-design source.
EOF

mkdir -p "$repository_root/.projects"
mv "$generated_contract" "$contract_file"
generated_contract=""
bash "$script_dir/validate-contract.sh" "$repository_root"
append_agents_pointer
success "Added .projects/project.md and the AGENTS.md starting point."

save_onboarding_files
print_chatgpt_setup
print_codex_setup
print_single_project_first_request

echo
success "Repository onboarding is complete."
echo "The initializer did not change any live GitHub issue or Project value."
