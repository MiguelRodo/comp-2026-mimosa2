#!/usr/bin/env bash

set -Eeuo pipefail

test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
skill_dir="$(cd "$test_dir/.." && pwd)"
validator="$skill_dir/scripts/validate-contract.sh"
setup="$skill_dir/scripts/setup.sh"
initializer="$skill_dir/scripts/init-project.sh"

bash -n "$validator"
bash -n "$setup"
bash -n "$initializer"

bash "$validator" "$test_dir/fixtures/single"
bash "$validator" "$test_dir/fixtures/single-user"
bash "$validator" "$test_dir/fixtures/dispatcher"
if bash "$validator" "$test_dir/fixtures/invalid" >/dev/null 2>&1; then
  echo "ERROR: lossy Priority mapping unexpectedly validated" >&2
  exit 1
fi
if bash "$validator" "$test_dir/fixtures/invalid-dispatcher" >/dev/null 2>&1; then
  echo "ERROR: mismatched dispatcher leaf unexpectedly validated" >&2
  exit 1
fi
if bash "$validator" "$test_dir/fixtures/invalid-colour" >/dev/null 2>&1; then
  echo "ERROR: unsupported option colour unexpectedly validated" >&2
  exit 1
fi
if bash "$validator" "$test_dir/fixtures/invalid-pending" >/dev/null 2>&1; then
  echo "ERROR: mixed pending and partial Priority mapping unexpectedly validated" >&2
  exit 1
fi

grep -Fqx 'name: github-project-admin' "$skill_dir/SKILL.md"
test -f "$skill_dir/README.md"
grep -Fq 'Set example#313 to P2.' "$test_dir/short-requests.md"
grep -Fq 'P3 | Low' "$skill_dir/SKILL.md"
grep -Fq 'Priority mapping status: pending' "$skill_dir/SKILL.md"
grep -Fq 'gh auth login --web --scopes "project,read:org"' "$skill_dir/README.md"
grep -Fq 'https://chatgpt.com/codex/settings/environments' "$skill_dir/README.md"
if grep -Eq '\$\{[^}]+(,,|\^\^)\}|(^|[^[:alnum:]_])(mapfile|readarray)([^[:alnum:]_]|$)|declare[[:space:]]+-A' \
  "$initializer"; then
  echo "ERROR: initializer uses a Bash feature newer than Bash 3.2" >&2
  exit 1
fi

test_tmp_dir="$(mktemp -d)"
cleanup() {
  if [[ -n "$test_tmp_dir" && "$test_tmp_dir" != "/" && -d "$test_tmp_dir" ]]; then
    rm -rf -- "$test_tmp_dir"
  fi
}
trap cleanup EXIT

mkdir -p "$test_tmp_dir/bin"
cat >"$test_tmp_dir/bin/gh" <<'EOF'
#!/usr/bin/env bash
set -Eeuo pipefail
case "${1:-}" in
  version)
    echo "gh version 2.98.0 (test)"
    ;;
  project)
    if [[ "${2:-}" == "--help" ]]; then
      exit 0
    elif [[ "${2:-}" == "list" ]]; then
      echo '{"projects":[]}'
    elif [[ "${2:-}" == "field-list" ]]; then
      printf 'Status\nPriority\nWorkstream\n'
    else
      exit 2
    fi
    ;;
  auth)
    [[ "${2:-}" == "status" ]]
    ;;
  api)
    if [[ "${2:-}" == "user" ]]; then
      echo "octocat"
    elif [[ "${2:-}" == users/* ]]; then
      case "${2#users/}" in
        octo-user) echo "User" ;;
        octo-org) echo "Organization" ;;
        *) exit 2 ;;
      esac
    elif [[ "${2:-}" == "graphql" ]]; then
      if [[ -n "${GH_GRAPHQL_LOG:-}" ]]; then
        printf '%s\n' "$*" >>"$GH_GRAPHQL_LOG"
      fi
      if [[ "$*" == *"user(login:"* && "$*" == *"organization(login:"* ]]; then
        echo "combined owner query is forbidden" >&2
        exit 90
      elif [[ "$*" == *"@tsv"* ]]; then
        printf '12\tExample planning\ttrue\n'
      elif [[ "$*" == *".number"* ]]; then
        if [[ "$*" == *"organization(login:"* ]]; then echo "12"; else echo "4"; fi
      elif [[ "$*" == *"organization(login:"* ]]; then
        echo "Example planning"
      else
        echo "User planning"
      fi
    else
      exit 2
    fi
    ;;
  repo)
    [[ "${2:-}" == "view" ]]
    if [[ "$*" == *"nameWithOwner,visibility"* ]]; then
      printf 'octo-org/example\tPUBLIC\n'
    elif [[ "$*" == *"octo-user/example"* ]]; then
      echo "octo-user/example"
    else
      echo "octo-org/example"
    fi
    ;;
  issue)
    [[ "${2:-}" == "list" ]]
    echo "[]"
    ;;
  skill)
    [[ "${2:-}" == "install" ]]
    if [[ -n "${GH_SKILL_LOG:-}" ]]; then
      printf '%s\n' "$*" >>"$GH_SKILL_LOG"
    fi
    ;;
  *)
    exit 2
    ;;
esac
EOF
chmod 0755 "$test_tmp_dir/bin/gh"

secret_value="test-token-must-not-appear"
PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  bash -x "$setup" --skip-install --no-contract \
  --repository octo-org/example \
  --project-owner octo-org --project-owner-type organization --project-number 12 \
  --project-title "Example planning" \
  --install-skill-from octo-org/project-skills --agent codex \
  >"$test_tmp_dir/setup.log" 2>&1

if grep -Fq "$secret_value" "$test_tmp_dir/setup.log"; then
  echo "ERROR: setup output exposed GH_TOKEN" >&2
  exit 1
fi
grep -Fq 'preflight passed' "$test_tmp_dir/setup.log"
grep -Fq 'Verified repository: octo-org/example.' "$test_tmp_dir/setup.log"
grep -Fq 'Verified Project: octo-org/12.' "$test_tmp_dir/setup.log"

PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  bash "$setup" --skip-install --no-contract --no-repository \
  --project-owner octo-user --project-number 4 \
  --project-title "User planning" \
  >"$test_tmp_dir/discovered-user.log" 2>&1
grep -Fq 'Verified Project: octo-user/4.' "$test_tmp_dir/discovered-user.log"

if PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  GH_GRAPHQL_LOG="$test_tmp_dir/mismatch-graphql.log" \
  bash "$setup" --skip-install --no-contract --no-repository \
  --project-owner octo-user --project-owner-type organization --project-number 4 \
  >"$test_tmp_dir/mismatch.log" 2>&1; then
  echo "ERROR: mismatched declared owner type unexpectedly passed" >&2
  exit 1
fi
grep -Fq 'declared Project owner type disagrees with GitHub' "$test_tmp_dir/mismatch.log"
if [[ -s "$test_tmp_dir/mismatch-graphql.log" ]]; then
  echo "ERROR: mismatched owner type reached the Project GraphQL query" >&2
  exit 1
fi

PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  GH_SKILL_LOG="$test_tmp_dir/local-skill.log" \
  bash "$setup" --skip-install --no-contract --no-repository \
  --install-skill-from "$test_dir/fixtures/single" \
  >"$test_tmp_dir/local-skill-output.log" 2>&1
grep -Fq -- "skill install $test_dir/fixtures/single github-project-admin --agent universal --scope user --force --from-local" \
  "$test_tmp_dir/local-skill.log"

PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  bash "$setup" --skip-install \
  --contract-root "$test_dir/fixtures/single" \
  >"$test_tmp_dir/contract-setup.log" 2>&1
grep -Fq 'Verified repository: octo-org/example.' "$test_tmp_dir/contract-setup.log"
grep -Fq 'Verified Project: octo-org/12.' "$test_tmp_dir/contract-setup.log"

PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  bash "$setup" --skip-install \
  --contract-root "$test_dir/fixtures/single-user" \
  >"$test_tmp_dir/user-contract-setup.log" 2>&1
grep -Fq 'Verified repository: octo-user/example.' "$test_tmp_dir/user-contract-setup.log"
grep -Fq 'Verified Project: octo-user/4.' "$test_tmp_dir/user-contract-setup.log"

mkdir -p "$test_tmp_dir/extend/.projects"
cat >"$test_tmp_dir/extend/.projects/setup.sh" <<'EOF'
#!/usr/bin/env bash
printf '%s:%s\n' "$PROJECTS_SETUP_MODE" "$PROJECTS_REPOSITORY_ROOT" >"$LOCAL_SETUP_LOG"
EOF
extend_setup_sha="$(sha256sum "$test_tmp_dir/extend/.projects/setup.sh" | awk '{print $1}')"
(
  cd "$test_tmp_dir/extend"
  PATH="$test_tmp_dir/bin:$PATH" LOCAL_SETUP_LOG="$test_tmp_dir/extend.log" \
    bash "$setup" --skip-install --no-contract --no-repository \
    >"$test_tmp_dir/extend-output.log" 2>&1
)
grep -Fq "extend:$test_tmp_dir/extend" "$test_tmp_dir/extend.log"
grep -Fq 'Running repository setup (extend)' "$test_tmp_dir/extend-output.log"
[[ "$(sha256sum "$test_tmp_dir/extend/.projects/setup.sh" | awk '{print $1}')" == "$extend_setup_sha" ]]

mkdir -p "$test_tmp_dir/override/.projects"
cat >"$test_tmp_dir/override/.projects/setup.sh" <<'EOF'
#!/usr/bin/env bash
# github-project-admin: override
printf '%s:%s\n' "$PROJECTS_SETUP_MODE" "$PROJECTS_REPOSITORY_ROOT" >"$LOCAL_SETUP_LOG"
EOF
(
  cd "$test_tmp_dir/override"
  LOCAL_SETUP_LOG="$test_tmp_dir/override.log" \
    bash "$setup" --skip-install --no-contract --no-repository \
    >"$test_tmp_dir/override-output.log" 2>&1
)
grep -Fq "override:$test_tmp_dir/override" "$test_tmp_dir/override.log"
grep -Fq 'Repository override setup completed.' "$test_tmp_dir/override-output.log"
if grep -Fq 'preflight passed' "$test_tmp_dir/override-output.log"; then
  echo "ERROR: override setup unexpectedly ran the shared preflight" >&2
  exit 1
fi

mkdir -p "$test_tmp_dir/init-single"
git -C "$test_tmp_dir/init-single" init -q
printf '# Existing repository guidance\n\nKeep this text.\n' >"$test_tmp_dir/init-single/AGENTS.md"
(
  cd "$test_tmp_dir/init-single"
  PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-output.log" 2>&1 <<'EOF'



12



EOF
)
bash "$validator" "$test_tmp_dir/init-single"
grep -Fq 'Keep this text.' "$test_tmp_dir/init-single/AGENTS.md"
grep -Fq '<!-- github-project-admin:start -->' "$test_tmp_dir/init-single/AGENTS.md"
grep -Fq '| Issue repository | octo-org/example |' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq '| Project title | Example planning |' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq '| Class | organization issue type | Issue Type |' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq '| Priority | pending live inspection | Priority |' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq '| Routing | Project membership; no routing label |' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fxq 'Priority mapping status: pending' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq 'This is a personal Project.' "$test_tmp_dir/init-single/.projects/project.md"
grep -Fq 'I will ask a few questions to configure' "$test_tmp_dir/init-output.log"
grep -Fq 'after /projects/ in its web address' "$test_tmp_dir/init-output.log"
grep -Fq 'GitHub user or organisation that owns the Project' "$initializer"
grep -Fq '1. Check GitHub' "$test_tmp_dir/init-output.log"
grep -Fq 'Use the repository with ChatGPT' "$test_tmp_dir/init-output.log"
grep -Fq 'Use the repository with Codex cloud' "$test_tmp_dir/init-output.log"
grep -Fq 'https://chatgpt.com/codex/settings/environments' "$test_tmp_dir/init-output.log"
grep -Fq 'Issue Type and' "$test_tmp_dir/init-output.log"
grep -Fq 'Give me an overview of what you plan to do' "$test_tmp_dir/init-output.log"
grep -Fq 'The files were left uncommitted.' "$test_tmp_dir/init-output.log"
if grep -Eq 'Current Project fields|following repository contract|Choose the provider.s Priority values|Does Project membership alone' \
  "$test_tmp_dir/init-output.log"; then
  echo "ERROR: initializer printed a removed diagnostic or question" >&2
  exit 1
fi
PATH="$test_tmp_dir/bin:$PATH" GH_TOKEN="$secret_value" \
  bash "$setup" --skip-install --contract-root "$test_tmp_dir/init-single" \
  >"$test_tmp_dir/init-contract-setup.log" 2>&1
grep -Fq 'Verified Project: octo-org/12.' "$test_tmp_dir/init-contract-setup.log"

init_contract_sha="$(sha256sum "$test_tmp_dir/init-single/.projects/project.md" | awk '{print $1}')"
init_agents_sha="$(sha256sum "$test_tmp_dir/init-single/AGENTS.md" | awk '{print $1}')"
(
  cd "$test_tmp_dir/init-single"
  printf '\n' | PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-rerun.log" 2>&1
)
grep -Fq 'existing repository setup was not replaced' "$test_tmp_dir/init-rerun.log"
[[ "$(sha256sum "$test_tmp_dir/init-single/.projects/project.md" | awk '{print $1}')" == "$init_contract_sha" ]]
[[ "$(sha256sum "$test_tmp_dir/init-single/AGENTS.md" | awk '{print $1}')" == "$init_agents_sha" ]]

mkdir -p "$test_tmp_dir/init-multiple"
git -C "$test_tmp_dir/init-multiple" init -q
(
  cd "$test_tmp_dir/init-multiple"
  printf '\nn\n\n' | PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-multiple.log" 2>&1
)
test ! -e "$test_tmp_dir/init-multiple/.projects/project.md"
grep -Fq '<!-- github-project-admin:start -->' "$test_tmp_dir/init-multiple/AGENTS.md"
grep -Fq 'as a personal' "$test_tmp_dir/init-multiple.log"
grep -Fq 'Finish the multi-Project routing' "$test_tmp_dir/init-multiple.log"
grep -Fq 'Add Project <number>' "$test_tmp_dir/init-multiple.log"

mkdir -p "$test_tmp_dir/init-collaborative-multiple"
git -C "$test_tmp_dir/init-collaborative-multiple" init -q
(
  cd "$test_tmp_dir/init-collaborative-multiple"
  printf 'y\nn\n\n' | PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-collaborative-multiple.log" 2>&1
)
test ! -e "$test_tmp_dir/init-collaborative-multiple/.projects/project.md"
grep -Fq 'as a collaborative' "$test_tmp_dir/init-collaborative-multiple.log"

mkdir -p "$test_tmp_dir/init-collaborative-single"
git -C "$test_tmp_dir/init-collaborative-single" init -q
(
  cd "$test_tmp_dir/init-collaborative-single"
  PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-collaborative-single.log" 2>&1 <<'EOF'
y


12

n
EOF
)
bash "$validator" "$test_tmp_dir/init-collaborative-single"
grep -Fq 'This is a collaborative Project.' \
  "$test_tmp_dir/init-collaborative-single/.projects/project.md"

mkdir -p "$test_tmp_dir/init-commit" "$test_tmp_dir/init-remote.git"
git -C "$test_tmp_dir/init-remote.git" init -q --bare
git -C "$test_tmp_dir/init-commit" init -q -b main
git -C "$test_tmp_dir/init-commit" config user.name "Setup Test"
git -C "$test_tmp_dir/init-commit" config user.email "setup@example.invalid"
printf 'base\n' >"$test_tmp_dir/init-commit/base.txt"
git -C "$test_tmp_dir/init-commit" add base.txt
git -C "$test_tmp_dir/init-commit" commit -qm "Base"
git -C "$test_tmp_dir/init-commit" remote add origin "$test_tmp_dir/init-remote.git"
git -C "$test_tmp_dir/init-commit" push -qu origin main
printf 'keep staged\n' >"$test_tmp_dir/init-commit/unrelated.txt"
git -C "$test_tmp_dir/init-commit" add unrelated.txt
(
  cd "$test_tmp_dir/init-commit"
  PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-commit.log" 2>&1 <<'EOF'



12
y
n
EOF
)
[[ "$(git -C "$test_tmp_dir/init-commit" log -1 --format=%s)" == \
   "Configure GitHub Project administration" ]]
grep -Fqx 'unrelated.txt' <(git -C "$test_tmp_dir/init-commit" diff --cached --name-only)
if git -C "$test_tmp_dir/init-commit" show --format= --name-only HEAD | grep -Fqx 'unrelated.txt'; then
  echo "ERROR: onboarding commit included an unrelated staged file" >&2
  exit 1
fi
[[ "$(git -C "$test_tmp_dir/init-commit" rev-parse HEAD)" == \
   "$(git -C "$test_tmp_dir/init-remote.git" rev-parse refs/heads/main)" ]]
grep -Fq 'Pushed the onboarding commit.' "$test_tmp_dir/init-commit.log"

mkdir -p "$test_tmp_dir/init-push-failure"
git -C "$test_tmp_dir/init-push-failure" init -q -b main
git -C "$test_tmp_dir/init-push-failure" config user.name "Setup Test"
git -C "$test_tmp_dir/init-push-failure" config user.email "setup@example.invalid"
(
  cd "$test_tmp_dir/init-push-failure"
  PATH="$test_tmp_dir/bin:$PATH" bash "$initializer" \
    >"$test_tmp_dir/init-push-failure.log" 2>&1 <<'EOF'



12
y
n
EOF
)
grep -Fq 'The commit is safe locally, but Git could not push it.' \
  "$test_tmp_dir/init-push-failure.log"
grep -Fq 'git push -u origin main' "$test_tmp_dir/init-push-failure.log"
[[ "$(git -C "$test_tmp_dir/init-push-failure" log -1 --format=%s)" == \
   "Configure GitHub Project administration" ]]

echo "github-project-admin tests passed"
