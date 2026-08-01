#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from fractions import Fraction

from trottercert.crosscheck import (
    small_exact_crosscheck,
    small_open_exact_crosscheck,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=2)
    parser.add_argument("--width", type=int)
    parser.add_argument("--height", type=int)
    parser.add_argument("--tolerance", type=int, default=10**6)
    args = parser.parse_args()
    if (args.width is None) != (args.height is None):
        parser.error("--width and --height must be supplied together")
    result = (
        small_exact_crosscheck(args.length, Fraction(1, args.tolerance))
        if args.width is None
        else small_open_exact_crosscheck(
            args.width,
            args.height,
            Fraction(1, args.tolerance),
        )
    )
    print(
        json.dumps(
            result,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
