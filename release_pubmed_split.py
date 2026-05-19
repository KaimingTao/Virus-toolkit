#!/usr/bin/env python3
"""Release and purge a virus-specific PubMed split dataset from git history."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Zip a virus PubMed split folder, create a GitHub release with gh, "
            "stop tracking the split files, rewrite git history to purge them, "
            "and force-push the rewritten branch and release tag."
        )
    )
    parser.add_argument(
        "virus_dir",
        type=Path,
        help="Virus directory such as HIV-1 or HCV.",
    )
    parser.add_argument(
        "--split-dir",
        type=Path,
        default=None,
        help="Explicit PubMed split directory. Defaults to the only *_pubmed_split directory under virus_dir.",
    )
    parser.add_argument(
        "--remote",
        default="origin",
        help="Git remote to push to (default: origin).",
    )
    parser.add_argument(
        "--branch",
        default=None,
        help="Git branch to force-push. Defaults to the current branch.",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Release tag name. Defaults to <virus>-<split-dir>-YYYY-MM-DD.",
    )
    parser.add_argument(
        "--release-title",
        default=None,
        help="GitHub release title. Defaults to '<virus> PubMed split CSVs'.",
    )
    parser.add_argument(
        "--release-notes",
        default=None,
        help="GitHub release notes. Defaults to a short generated note.",
    )
    parser.add_argument(
        "--commit-message",
        default=None,
        help="Commit message for removing the split folder from tracking.",
    )
    parser.add_argument(
        "--yes",
        action="store_true",
        help="Skip the destructive-action confirmation prompt.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned actions without changing the repository.",
    )
    return parser.parse_args(argv)


def run(
    args: list[str],
    *,
    cwd: Path,
    dry_run: bool,
    env: dict[str, str] | None = None,
) -> None:
    print(f"+ {' '.join(args)}")
    if dry_run:
        return
    subprocess.run(args, cwd=cwd, check=True, env=env)


def git_output(repo_root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def resolve_repo_root(cwd: Path) -> Path:
    try:
        return Path(git_output(cwd, "rev-parse", "--show-toplevel")).resolve()
    except subprocess.CalledProcessError as exc:
        raise SystemExit("This script must run inside a git repository.") from exc


def ensure_clean_tracked_state(repo_root: Path) -> None:
    status = git_output(repo_root, "status", "--porcelain")
    tracked_changes = [
        line for line in status.splitlines() if line and not line.startswith("?? ")
    ]
    if tracked_changes:
        joined = "\n".join(tracked_changes)
        raise SystemExit(
            "Refusing to proceed with tracked working tree changes:\n"
            f"{joined}\n"
            "Commit, stash, or discard them first."
        )


def resolve_split_dir(virus_dir: Path, split_dir: Path | None) -> Path:
    if split_dir is not None:
        resolved = split_dir.resolve()
        if not resolved.is_dir():
            raise SystemExit(f"Split directory not found: {resolved}")
        return resolved

    matches = sorted(path for path in virus_dir.iterdir() if path.is_dir() and path.name.endswith("_pubmed_split"))
    if not matches:
        raise SystemExit(f"No *_pubmed_split directory found under {virus_dir}")
    if len(matches) > 1:
        names = ", ".join(path.name for path in matches)
        raise SystemExit(
            f"Multiple *_pubmed_split directories found under {virus_dir}: {names}. "
            "Pass --split-dir explicitly."
        )
    return matches[0].resolve()


def ensure_gh_available() -> None:
    if shutil.which("gh") is None:
        raise SystemExit("GitHub CLI `gh` is required but not installed or not on PATH.")


def create_zip(split_dir: Path, zip_path: Path, dry_run: bool) -> None:
    print(f"Creating archive: {zip_path}")
    if dry_run:
        return

    if zip_path.exists():
        zip_path.unlink()

    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for file_path in sorted(path for path in split_dir.rglob("*") if path.is_file()):
            archive.write(file_path, file_path.relative_to(split_dir.parent))


def ensure_ignore_rules(repo_root: Path, split_dir: Path, zip_path: Path, dry_run: bool) -> None:
    gitignore_path = repo_root / ".gitignore"
    lines = gitignore_path.read_text(encoding="utf-8").splitlines()

    split_rule = f"{split_dir.relative_to(repo_root).as_posix()}/"
    zip_rule = zip_path.relative_to(repo_root).as_posix()
    additions = [rule for rule in (split_rule, zip_rule) if rule not in lines]
    if not additions:
        return

    print(f"Updating ignore rules in {gitignore_path}")
    if dry_run:
        return

    new_text = gitignore_path.read_text(encoding="utf-8")
    if new_text and not new_text.endswith("\n"):
        new_text += "\n"
    new_text += "\n".join(additions) + "\n"
    gitignore_path.write_text(new_text, encoding="utf-8")


def confirm_or_exit(message: str, assume_yes: bool) -> None:
    if assume_yes:
        return
    response = input(f"{message} [y/N]: ").strip().lower()
    if response not in {"y", "yes"}:
        raise SystemExit("Aborted.")


def build_defaults(virus_dir: Path, split_dir: Path, branch: str | None) -> tuple[str, str, str, str]:
    branch_name = branch or git_output(virus_dir.parent, "branch", "--show-current")
    tag = f"{virus_dir.name.lower()}-{split_dir.name.replace('_', '-')}-{date.today().isoformat()}"
    title = f"{virus_dir.name} PubMed split CSVs"
    notes = f"Archived release asset for {split_dir.relative_to(virus_dir.parent).as_posix()}."
    return branch_name, tag, title, notes


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    repo_root = resolve_repo_root(Path.cwd().resolve())
    ensure_clean_tracked_state(repo_root)
    ensure_gh_available()

    virus_dir = (repo_root / args.virus_dir).resolve()
    if not virus_dir.is_dir():
        raise SystemExit(f"Virus directory not found: {virus_dir}")
    if repo_root not in virus_dir.parents:
        raise SystemExit(f"Virus directory must be inside the repository: {virus_dir}")

    split_dir = resolve_split_dir(virus_dir, args.split_dir)
    if repo_root not in split_dir.parents:
        raise SystemExit(f"Split directory must be inside the repository: {split_dir}")

    zip_path = split_dir.with_suffix(".zip")
    branch_name, default_tag, default_title, default_notes = build_defaults(
        virus_dir, split_dir, args.branch
    )
    tag_name = args.tag or default_tag
    release_title = args.release_title or default_title
    release_notes = args.release_notes or default_notes
    commit_message = args.commit_message or f"Stop tracking {virus_dir.name} PubMed split CSVs"

    print(f"repo_root={repo_root}")
    print(f"virus_dir={virus_dir.relative_to(repo_root)}")
    print(f"split_dir={split_dir.relative_to(repo_root)}")
    print(f"zip_path={zip_path.relative_to(repo_root)}")
    print(f"branch={branch_name}")
    print(f"tag={tag_name}")

    confirm_or_exit(
        "This will create a release, rewrite git history for all refs, and force-push the result.",
        args.yes,
    )

    create_zip(split_dir, zip_path, args.dry_run)
    ensure_ignore_rules(repo_root, split_dir, zip_path, args.dry_run)

    run(["git", "tag", "-f", tag_name, "HEAD"], cwd=repo_root, dry_run=args.dry_run)
    run(
        [
            "gh",
            "release",
            "create",
            tag_name,
            str(zip_path.relative_to(repo_root)),
            "--title",
            release_title,
            "--notes",
            release_notes,
        ],
        cwd=repo_root,
        dry_run=args.dry_run,
    )
    run(
        ["git", "rm", "--cached", "-r", "--ignore-unmatch", str(split_dir.relative_to(repo_root))],
        cwd=repo_root,
        dry_run=args.dry_run,
    )
    run(["git", "add", ".gitignore"], cwd=repo_root, dry_run=args.dry_run)
    run(["git", "commit", "-m", commit_message], cwd=repo_root, dry_run=args.dry_run)
    run(
        [
            "git",
            "filter-branch",
            "--force",
            "--index-filter",
            f"git rm -r --cached --ignore-unmatch {split_dir.relative_to(repo_root).as_posix()}",
            "--prune-empty",
            "--tag-name-filter",
            "cat",
            "--",
            "--all",
        ],
        cwd=repo_root,
        dry_run=args.dry_run,
        env={
            **os.environ,
            "FILTER_BRANCH_SQUELCH_WARNING": "1",
        },
    )

    refs_original = repo_root / ".git" / "refs" / "original"
    if refs_original.exists():
        print(f"Removing backup refs: {refs_original}")
        if not args.dry_run:
            shutil.rmtree(refs_original)

    run(
        ["git", "reflog", "expire", "--expire=now", "--all"],
        cwd=repo_root,
        dry_run=args.dry_run,
    )
    run(
        ["git", "gc", "--prune=now", "--aggressive"],
        cwd=repo_root,
        dry_run=args.dry_run,
    )
    run(["git", "tag", "-f", tag_name, "HEAD"], cwd=repo_root, dry_run=args.dry_run)
    run(
        ["git", "push", "--force-with-lease", args.remote, branch_name],
        cwd=repo_root,
        dry_run=args.dry_run,
    )
    run(
        ["git", "push", "--force", args.remote, f"refs/tags/{tag_name}"],
        cwd=repo_root,
        dry_run=args.dry_run,
    )

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
