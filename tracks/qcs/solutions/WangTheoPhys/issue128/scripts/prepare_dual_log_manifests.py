#!/usr/bin/env python3
"""Prepare or verify canonical exact free-log word manifests."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from trottercert.cubic_field import fourth_order_suzuki_cubic_stages
from trottercert.dual_word_manifest import (
    load_manifest_index,
    verify_manifest_index,
    write_manifest_set,
)

ISSUE_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATHS = (
    ISSUE_ROOT / "src/trottercert/cubic_field.py",
    ISSUE_ROOT / "src/trottercert/cubic_local.py",
    ISSUE_ROOT / "src/trottercert/dual_word_manifest.py",
)


def _source_hashes() -> dict[str, str]:
    return {
        str(path.relative_to(ISSUE_ROOT)): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in SOURCE_PATHS
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()

    if args.verify is not None:
        if any(
            value is not None
            for value in (args.degree, args.shard_count, args.output_dir)
        ):
            raise SystemExit("--verify cannot be combined with generation options")
        index = load_manifest_index(args.verify)
        verify_manifest_index(index, args.verify.parent)
        print(
            "dual word manifest set valid: "
            f"degree={index.degree} groups={index.total_groups} "
            f"words={index.total_words} shards={index.shard_count}",
            flush=True,
        )
        return

    if args.degree is None or args.shard_count is None or args.output_dir is None:
        raise SystemExit(
            "generation requires --degree, --shard-count, and --output-dir"
        )
    index_path = write_manifest_set(
        args.output_dir,
        fourth_order_suzuki_cubic_stages(4),
        args.degree,
        args.shard_count,
        implementation_sources=_source_hashes(),
    )
    index = load_manifest_index(index_path)
    print(
        "wrote dual word manifest set: "
        f"degree={index.degree} groups={index.total_groups} "
        f"words={index.total_words} shards={index.shard_count} "
        f"index={index_path}",
        flush=True,
    )


if __name__ == "__main__":
    main()
