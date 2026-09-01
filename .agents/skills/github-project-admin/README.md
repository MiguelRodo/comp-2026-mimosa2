# GitHub Project admin

This skill lets an agent manage GitHub issues and Projects from ordinary requests. It checks the current state before changing anything and checks the result afterwards.

## Set up a repository

### 1. Create or find the GitHub Project

Open the **Projects** tab on your GitHub profile or organisation. Create the Project if it does not exist yet.

The Project number is the number after `/projects/` in its web address. For example, these are Project 40 and Project 12:

```text
https://github.com/users/example/projects/40
https://github.com/orgs/example-org/projects/12
```

Have the number for each relevant Project ready before starting.

### 2. Install the skill and answer the questions

From the repository, run:

```bash
gh skill install MiguelRodo/projects github-project-admin \
  --agent universal \
  --scope project

bash .agents/skills/github-project-admin/scripts/init-project.sh
```

The initializer asks whether the setup is personal or collaborative and whether the repository uses one or several Projects. It discovers privacy, repository identity, Project title and owner type from GitHub.

- For one Project, it generates and validates `.projects/project.md` and adds a bounded `AGENTS.md` section.
- For several Projects, it returns a tailored request for an agent to ask only for the Project numbers and routing choices it cannot discover safely.

It also asks whether ChatGPT or Codex should leave existing issues alone, suggest Class and Workstream values, or propose them and wait for approval before applying them. The initializer itself never changes live issues or Project fields.

Review and commit:

```text
.agents/skills/github-project-admin/
.projects/
AGENTS.md
```

### 3. Choose how you want to use it

#### ChatGPT or another conversational service

Create a provider Project or workspace with access to the repository. Use this standing instruction:

> For work concerning a GitHub repository, especially reading or updating GitHub issues or Projects, first retrieve and follow the target repository's `AGENTS.md`. Follow the skill and configuration files it references. If the repository or `AGENTS.md` is unavailable, say so rather than guessing.

Ask for the result you want. If the service cannot make the change, it should return a concise command block for you to paste into a terminal.

#### Codex Cloud or another shell-capable agent

Create an environment for the repository and use:

```bash
bash .agents/skills/github-project-admin/scripts/setup.sh
```

Provide `GH_TOKEN` as an environment variable with repository and Project access. Allow `github.com` and `api.github.com`. The agent can then run and verify the commands directly.

## Repository shapes

The onboarding flow handles all four combinations:

| Governance | Topology | Onboarding result |
| --- | --- | --- |
| Personal | One Project | Generate the contract directly |
| Collaborative | One Project | Generate the contract directly |
| Personal | Several Projects | Give a personal multi-Project routing handoff |
| Collaborative | Several Projects | Give a collaborative multi-Project routing handoff |

Privacy is discovered separately from GitHub. A private repository can still be collaborative, and a public repository can still use a Project managed by one person.

## Repository-specific setup

If the repository needs extra tools, add `.projects/setup.sh`. It runs automatically after the shared setup and is not replaced when the skill is updated.

To replace common setup completely, place this within the first 20 lines:

```bash
# github-project-admin: override
```

## Update the skill

Run this inside the repository, then commit the changed skill files:

```bash
gh skill update github-project-admin
```

## If something fails

Paste the terminal output into the chat, or say which command failed. The agent should inspect what already succeeded and give you only the corrected or remaining commands.

If the failure is reusable, the agent may offer to improve `MiguelRodo/projects`. It should open an issue or pull request only after you agree.
