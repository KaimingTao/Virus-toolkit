# AGENTS.md

## General behavior

- Prefer simple, maintainable, and readable code over clever code.
- Make the smallest safe change that solves the problem.
- Preserve existing project structure and conventions unless explicitly asked to refactor.
- Explain important design decisions briefly.
- Ask before introducing major new dependencies.
- Do not silently ignore errors or exceptions.
- Do not hardcode machine-specific paths or credentials.
- Prefer deterministic and reproducible workflows.

---

## Repository inspection

Before making changes:

- Inspect the repository structure first.
- Identify the existing language, package manager, formatter, linter, and test framework.
- Follow the repository's existing style when possible.
- Reuse existing utilities and helper functions before creating new ones.

---

## Python environment requirements

- Prefer `uv` for Python project and dependency management.
- If the repository does not contain a `pyproject.toml`, initialize the project with:

```bash
uv init
