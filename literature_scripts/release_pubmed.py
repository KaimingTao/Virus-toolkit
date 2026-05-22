#!/usr/bin/env python3
"""Create a GitHub release for a virus-specific PubMed dataset."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Zip a virus PubMed folder and create a GitHub release with gh."
        )
    )
    parser.add_argument(
        "virus_dir",
        type=Path,
        help="Virus directory such as HIV-1 or HCV.",
    )
    parser.add_argument(
        "--pubmed",
        type=Path,
        default=None,
        help="Explicit PubMed directory. Defaults to virus_dir/pubmed_search.",
    )
    parser.add_argument(
        "--tag",
        default=None,
        help="Release tag name. Defaults to <virus>-<pubmed-dir>-YYYY-MM-DD.",
    )
    parser.add_argument(
        "--release-title",
        default=None,
        help="GitHub release title. Defaults to '<virus> PubMed CSVs'.",
    )
    parser.add_argument(
        "--release-notes",
        default=None,
        help="GitHub release notes. Defaults to a short generated note.",
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


def resolve_repo_root(cwd: Path) -> Path:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--show-toplevel"],
            cwd=cwd,
            check=True,
            capture_output=True,
            text=True,
        )
        return Path(result.stdout.strip()).resolve()
    except subprocess.CalledProcessError as exc:
        raise SystemExit("This script must run inside a git repository.") from exc


def git_output(repo_root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def resolve_pubmed_dir(virus_dir: Path, pubmed_dir: Path | None) -> Path:
    if pubmed_dir is not None:
        resolved = pubmed_dir.resolve()
        if not resolved.is_dir():
            raise SystemExit(f"PubMed directory not found: {resolved}")
        return resolved

    resolved = (virus_dir / "pubmed_search").resolve()
    if not resolved.is_dir():
        raise SystemExit(f"Default PubMed directory not found: {resolved}")
    return resolved


def ensure_gh_available() -> None:
    if shutil.which("gh") is None:
        raise SystemExit("GitHub CLI `gh` is required but not installed or not on PATH.")


def create_zip(pubmed_dir: Path, zip_path: Path, dry_run: bool) -> None:
    print(f"Creating archive: {zip_path}")
    if dry_run:
        return

    if zip_path.exists():
        zip_path.unlink()

    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
        for file_path in sorted(path for path in pubmed_dir.rglob("*") if path.is_file()):
            archive.write(file_path, file_path.relative_to(pubmed_dir.parent))


def confirm_or_exit(message: str, assume_yes: bool) -> None:
    if assume_yes:
        return
    response = input(f"{message} [y/N]: ").strip().lower()
    if response not in {"y", "yes"}:
        raise SystemExit("Aborted.")


def build_defaults(virus_dir: Path, pubmed_dir: Path) -> tuple[str, str, str]:
    tag = f"{virus_dir.name.lower()}-{pubmed_dir.name.replace('_', '-')}-{date.today().isoformat()}"
    title = f"{virus_dir.name} PubMed CSVs"
    notes = f"Archived release asset for {pubmed_dir.relative_to(virus_dir.parent).as_posix()}."
    return tag, title, notes


def resolve_release_target(repo_root: Path) -> str:
    branch = git_output(repo_root, "rev-parse", "--abbrev-ref", "HEAD")
    if branch and branch != "HEAD":
        return branch
    return git_output(repo_root, "rev-parse", "HEAD")


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    repo_root = resolve_repo_root(Path.cwd().resolve())
    ensure_gh_available()

    virus_dir = (repo_root / args.virus_dir).resolve()
    if not virus_dir.is_dir():
        raise SystemExit(f"Virus directory not found: {virus_dir}")
    if repo_root not in virus_dir.parents:
        raise SystemExit(f"Virus directory must be inside the repository: {virus_dir}")

    pubmed_dir = resolve_pubmed_dir(virus_dir, args.pubmed)
    if repo_root not in pubmed_dir.parents:
        raise SystemExit(f"PubMed directory must be inside the repository: {pubmed_dir}")

    zip_path = pubmed_dir.with_name(f"{virus_dir.name}_{pubmed_dir.name}.zip")
    default_tag, default_title, default_notes = build_defaults(virus_dir, pubmed_dir)
    tag_name = args.tag or default_tag
    release_title = args.release_title or default_title
    release_notes = args.release_notes or default_notes
    release_target = resolve_release_target(repo_root)

    print(f"repo_root={repo_root}")
    print(f"virus_dir={virus_dir.relative_to(repo_root)}")
    print(f"pubmed_dir={pubmed_dir.relative_to(repo_root)}")
    print(f"zip_path={zip_path.relative_to(repo_root)}")
    print(f"tag={tag_name}")
    print(f"target={release_target}")

    confirm_or_exit(
        "This will create a zip archive and publish a GitHub release.",
        args.yes,
    )

    create_zip(pubmed_dir, zip_path, args.dry_run)

    run(
        [
            "gh",
            "release",
            "create",
            tag_name,
            "--target",
            release_target,
            str(zip_path.relative_to(repo_root)),
            "--title",
            release_title,
            "--notes",
            release_notes,
        ],
        cwd=repo_root,
        dry_run=args.dry_run,
    )

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
