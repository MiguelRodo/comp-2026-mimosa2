#!/usr/bin/env bash

# Guided, repository-local onboarding for github-project-admin. This script
# writes repository configuration only. It never changes live GitHub issues or
# Project fields.
set +x
set -Eeuo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
generated_contract=""

die() {
  echo "ERROR: $*" >&2
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
  read -r -p "$prompt [$default]: " answer || die "input ended before setup was complete"
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
    read -r -p "$prompt [$hint]: " answer || die "input ended before setup was complete"
    answer="${answer:-$default}"
    case "${answer,,}" in
      y|yes) return 0 ;;
      n|no) return 1 ;;
      *) echo "Please answer yes or no." >&2 ;;
    esac
  done
}

ask_choice() {
  local prompt="$1" minimum="$2" maximum="$3" answer
  while true; do
    read -r -p "$prompt: " answer || die "input ended before setup was complete"
    if [[ "$answer" =~ ^[0-9]+$ ]] && ((answer >= minimum && answer <= maximum)); then
      printf '%s' "$answer"
      return 0
    fi
    echo "Choose a number from $minimum to $maximum." >&2
  done
}

require_text() {
  local prompt="$1" answer
  while true; do
    read -r -p "$prompt: " answer || die "input ended before setup was complete"
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
    echo "AGENTS.md already contains the GitHub Project routing section."
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
  echo "Added the GitHub Project routing section to AGENTS.md."
}

command -v git >/dev/null 2>&1 || die "git is required"
command -v gh >/dev/null 2>&1 || die "gh is required; install the skill with GitHub CLI first"
git rev-parse --is-inside-work-tree >/dev/null 2>&1 ||
  die "run this command inside the repository you want to configure"
repository_root="$(git rev-parse --show-toplevel)"
cd "$repository_root"

gh auth status >/dev/null 2>&1 ||
  die "GitHub is not authenticated; authenticate gh or configure GH_TOKEN"

contract_file="$repository_root/.projects/project.md"
if [[ -e "$contract_file" ]]; then
  echo "A repository contract already exists at .projects/project.md."
  echo "It was not changed. Ask an agent to review or migrate it if necessary."
  if [[ -x "$script_dir/validate-contract.sh" || -f "$script_dir/validate-contract.sh" ]]; then
    bash "$script_dir/validate-contract.sh" "$repository_root"
  fi
  append_agents_pointer
  exit 0
fi

repository_record="$(gh repo view --json nameWithOwner,visibility \
  --jq '[.nameWithOwner,.visibility] | @tsv')" ||
  die "could not discover the current GitHub repository"
IFS=$'\t' read -r repository visibility <<<"$repository_record"
[[ "$repository" =~ ^[^/[:space:]]+/[^/[:space:]]+$ ]] ||
  die "GitHub returned an invalid repository identity"
repository_owner="${repository%%/*}"

echo "Repository: $repository"
echo
if ask_yes_no "Will people other than you collaborate on Projects managed from this repository?" no; then
  governance="collaborative"
else
  governance="personal"
fi

if ! ask_yes_no "Does this repository use exactly one GitHub Project?" yes; then
  echo
  echo "This is a $governance, multi-Project repository."
  echo "No repository configuration was generated because each issue needs one"
  echo "deterministic route to the correct Project. An agent can set that up without"
  echo "asking you to hand-write dispatcher files."
  echo
  echo "Ask ChatGPT or Codex:"
  echo
  echo "  Set up $repository as a $governance multi-Project repository. Use the"
  echo "  installed github-project-admin skill, preserve the existing AGENTS.md,"
  echo "  discover what GitHub can, and ask me only for Project numbers and routing"
  echo "  choices that cannot be determined safely. Do not change live issues or"
  echo "  Projects until I approve the generated repository configuration."
  exit 2
fi

project_owner="$(ask_default "GitHub login that owns the Project" "$repository_owner")"
[[ "$project_owner" =~ ^[A-Za-z0-9_.-]+$ ]] || die "invalid Project owner login"

project_number="$(require_text "Project number")"
[[ "$project_number" =~ ^[1-9][0-9]*$ ]] || die "Project number must be a positive integer"

observed_owner_type="$(gh api "users/$project_owner" --jq .type 2>/dev/null)" ||
  die "could not discover whether the Project owner is a person or organisation"
case "$observed_owner_type" in
  User)
    project_owner_type="user"
    project_selector='.data.user.projectV2'
    project_query='query($login: String!, $number: Int!) { user(login: $login) { projectV2(number: $number) { number title public } } }'
    ;;
  Organization)
    project_owner_type="organization"
    project_selector='.data.organization.projectV2'
    project_query='query($login: String!, $number: Int!) { organization(login: $login) { projectV2(number: $number) { number title public } } }'
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

echo "Project: $project_title ($project_owner/$project_number)"
echo "Owner type: $project_owner_type (discovered automatically)"
field_names="$(gh project field-list "$project_number" --owner "$project_owner" \
  --format json --jq '.fields[].name' 2>/dev/null |
  awk 'BEGIN { separator = "" } { printf "%s%s", separator, $0; separator = ", " } END { print "" }')"
if [[ -n "$field_names" ]]; then
  echo "Current Project fields: $field_names"
fi
echo

if [[ "$project_owner_type" == "organization" ]] &&
   ask_yes_no "Does this Project use organisation Issue Type and Priority fields?" yes; then
  class_location="organization issue type"
  class_field="Issue Type"
  priority_location="organization issue field"
  priority_field="Priority"
else
  class_location="project field"
  class_field="Class"
  priority_location="project field"
  priority_field="Priority"
fi

if ask_yes_no "Does Project membership alone determine which issues belong here?" yes; then
  routing="Project membership; no routing label"
else
  routing_label="$(require_text "Exact routing label, including any prefix")"
  routing="label:$routing_label"
fi

echo
echo "Priority mapping:"
echo "  1. Urgent, High, Medium, Low (recommended default)"
echo "  2. P0, P1, P2, P3"
echo "  3. Enter four provider values"
priority_choice="$(ask_choice "Choose the provider's Priority values" 1 3)"
case "$priority_choice" in
  1) p0="Urgent"; p1="High"; p2="Medium"; p3="Low" ;;
  2) p0="P0"; p1="P1"; p2="P2"; p3="P3" ;;
  3)
    p0="$(require_text "Provider value for P0")"
    p1="$(require_text "Provider value for P1")"
    p2="$(require_text "Provider value for P2")"
    p3="$(require_text "Provider value for P3")"
    ;;
esac
[[ "$p0" != "$p1" && "$p0" != "$p2" && "$p0" != "$p3" &&
   "$p1" != "$p2" && "$p1" != "$p3" && "$p2" != "$p3" ]] ||
  die "Priority provider values must be distinct"

echo
echo "How should an agent handle the existing issues?"
echo "  1. Do not classify them now"
echo "  2. Suggest Class and Workstream values only"
echo "  3. Propose values, wait for approval, then apply and verify them"
organisation_choice="$(ask_choice "Choose the next step" 1 3)"

if [[ "${visibility^^}" == "PRIVATE" ]]; then
  privacy="private repository"
elif [[ "$project_public" == "false" ]]; then
  privacy="${visibility,,} repository with a private Project"
else
  privacy="${visibility,,} repository"
fi

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
| Routing | $routing |
| Privacy | $privacy |

## Field locations

| Common dimension | Provider location | Provider field |
| --- | --- | --- |
| Class | $class_location | $class_field |
| Priority | $priority_location | $priority_field |
| Status | project field | Status |
| Workstream | project field | Workstream |
| Due date | project field | Target date |
| Parent | native issue relationship | Parent issue |

## Priority mapping

| Common value | Provider value |
| --- | --- |
| P0 | $p0 |
| P1 | $p1 |
| P2 | $p2 |
| P3 | $p3 |

## Class and Workstream

Class and Workstream option sets are intentionally not fixed by onboarding.
An agent may inspect the existing issues and suggest a concise, useful vocabulary
before the live Project is changed.

## Governance

- This is a $governance Project.
- Project membership or the exact routing rule above determines Project scope.
- Labels must not duplicate Class, Priority, Status or Workstream.
- Assignment is explicit only unless a later repository decision says otherwise.
- Exact requested administration requires no separate scope-design source.
EOF

echo
echo "The following repository contract will be created:"
echo
sed 's/^/  /' "$generated_contract"
echo
if ! ask_yes_no "Create .projects/project.md and add the AGENTS.md routing section?" yes; then
  echo "No files were changed."
  exit 0
fi

mkdir -p "$repository_root/.projects"
mv "$generated_contract" "$contract_file"
generated_contract=""
bash "$script_dir/validate-contract.sh" "$repository_root"
append_agents_pointer

echo
echo "Repository onboarding files are ready. Review and commit:"
echo
echo "  git diff -- .projects/project.md AGENTS.md .agents/skills/github-project-admin"
echo
echo "For a ChatGPT Project, use this stable instruction:"
echo
echo "  For work concerning a GitHub repository, especially reading or updating"
echo "  GitHub issues or Projects, first retrieve and follow the target"
echo "  repository's AGENTS.md. Follow the skill and configuration it references."
echo
case "$organisation_choice" in
  1)
    echo "No issue-classification prompt was selected."
    ;;
  2)
    echo "Give ChatGPT or Codex this request:"
    echo
    echo "  Start from AGENTS.md. Review the current issues and suggest sensible Class"
    echo "  and Workstream options and values. Do not change GitHub yet."
    ;;
  3)
    echo "Give ChatGPT or Codex this request:"
    echo
    echo "  Start from AGENTS.md. Review the current issues and propose sensible Class"
    echo "  and Workstream options and values. Wait for my approval, then apply and"
    echo "  independently verify the approved changes."
    ;;
esac

echo
echo "The onboarding script did not change any live GitHub issue or Project value."
