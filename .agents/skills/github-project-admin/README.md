# GitHub Project admin

This skill helps ChatGPT, Codex and other agents understand and manage a repository's GitHub issues and Projects. It checks the current state before changing anything and checks the result afterwards.

## Before you start

Install [GitHub CLI](https://cli.github.com/) 2.96 or newer, then check it:

```text
gh --version
```

If you have not signed in, run:

```text
gh auth login --web --scopes "project,read:org"
```

If you are already signed in but need Project access, run:

```text
gh auth refresh --scopes "project,read:org"
```

Check the result:

```text
gh auth status
gh api user --jq .login
```

The commands in this guide work in Bash, Git Bash and WSL. They also work from PowerShell when `bash` is installed and available as a command.

## 1. Create or find the GitHub Project

Open the **Projects** tab on your GitHub profile or organisation. Create the Project if needed.

Its number is after `/projects/` in the web address:

```text
https://github.com/users/example/projects/40       Project number 40
https://github.com/orgs/example-org/projects/12    Project number 12
```

Have the number for each relevant Project ready.

## 2. Install and configure the repository

From the repository, run these one-line commands:

```text
gh skill install MiguelRodo/projects github-project-admin --agent universal --scope project
bash .agents/skills/github-project-admin/scripts/init-project.sh
```

The initializer first explains that it will configure the repository so chats and agents can understand the Project. It discovers GitHub facts and asks only about collaboration, whether the repository uses one or several Projects, and the owner and number of a single Project.

This covers personal or collaborative repositories with one Project or several. Repository and Project privacy are discovered separately from GitHub.

For one Project, it creates `.projects/project.md` and adds a small section to `AGENTS.md`. It does not ask you to edit them or print the full contract.

For several Projects, it prepares the repository and gives you an agent request to finish the routing. The correct route for each issue cannot be guessed safely.

The initializer does not change live issues or Project fields. It leaves the existing Priority field and options alone, even if they do not include P3. Priority is marked as pending until an agent confirms its location, inspects the options and records a complete mapping.

It then asks whether it may stage, commit and push only the onboarding files. A failed commit or push leaves the work in a recoverable local state and prints the next command.

## 3. Use it from ChatGPT

Create or open a [ChatGPT Project](https://chatgpt.com/projects), make the repository available to it, and paste this into the Project instructions:

> For work concerning a GitHub repository, especially reading or updating GitHub issues or Projects, first retrieve and follow the target repository's `AGENTS.md`. Follow the skill and configuration files it references. If the repository or `AGENTS.md` is unavailable, say so rather than guessing.
>
> Treat my prompt as the desired outcome. If this chat cannot make a required GitHub change, return the smallest safe command block for me to paste into a terminal, including a check of the result.

Ask for the outcome you want. If ChatGPT cannot make a change itself, it should return a short terminal-ready command block.

## 4. Use it from Codex cloud

Open [Codex environments](https://chatgpt.com/codex/settings/environments), create an environment and choose the repository. Use:

```text
bash .agents/skills/github-project-admin/scripts/setup.sh
```

Create a [classic GitHub personal access token](https://github.com/settings/tokens/new) with an expiry and the `repo`, `read:org` and `project` scopes. Authorise it for organisation SSO if required.

Add it to the Codex environment as an environment variable named `GH_TOKEN`, not a setup-only secret. Enable agent internet access and allow:

```text
github.com
api.github.com
```

See the [official Codex environment guide](https://developers.openai.com/codex/environments/cloud-environment) for how environment variables, setup and agent internet access work.

## 5. Start with the current issues

The initializer offers a tailored first request after the ChatGPT and Codex instructions. For a single Project, the request can:

- confirm the local Priority location and mapping from the existing live field without changing it;
- set up or refine Issue Type or Class and Workstream, with sensible colours;
- organise existing issues and useful native parent/sub-issue relationships.

Choose whether the agent should show you a plan first or carry out and verify the work immediately.

For several Projects, first use the generated request to define the dispatcher and routing. To add another Project later, ask:

> Start from `AGENTS.md`. Add Project 12 owned by `example-org` to this repository's existing multi-Project configuration. Preserve current routes, ask me only for the new routing decision, and validate the result.

## Repository-specific setup

If the repository needs extra tools, add `.projects/setup.sh`. It runs automatically after the shared setup and is not replaced when the skill is updated.

To replace common setup completely, place this within the first 20 lines:

```bash
# github-project-admin: override
```

## Update the skill

Run this inside the repository, then commit the changed skill files:

```text
gh skill update github-project-admin
```

The update does not replace `.projects/project.md` or `.projects/setup.sh`.

## If something fails

Paste the terminal output into the chat, or say which section failed. The agent should inspect what already succeeded and give you only the corrected or remaining commands.

If the failure is reusable, the agent may offer to improve `MiguelRodo/projects`. It should open an issue or pull request only after you agree.
