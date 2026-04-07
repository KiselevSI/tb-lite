#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate text tables listed in a manifest file."
    )
    parser.add_argument(
        "--input-list",
        required=True,
        help="Text file with one input path per line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file path.",
    )
    parser.add_argument(
        "--keep-header",
        action="store_true",
        help="Keep the first line only from the first non-empty file.",
    )
    parser.add_argument(
        "--prepend-line",
        help="Write this line to the output before concatenating inputs.",
    )
    return parser.parse_args()


def read_input_list(path: Path) -> list[Path]:
    items: list[Path] = []
    seen: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            key = str(Path(line))
            if key in seen:
                continue

            seen.add(key)
            items.append(Path(line))

    return items


def main():
    args = parse_args()

    inputs = read_input_list(Path(args.input_list))
    kept_header = False

    with Path(args.output).open("w", encoding="utf-8") as out_handle:
        if args.prepend_line is not None:
            out_handle.write(args.prepend_line)
            out_handle.write("\n")

        for path in inputs:
            with path.open("r", encoding="utf-8", errors="replace") as in_handle:
                if args.keep_header:
                    first_line = in_handle.readline()
                    if not first_line:
                        continue
                    if not kept_header:
                        out_handle.write(first_line)
                        kept_header = True
                    for line in in_handle:
                        out_handle.write(line)
                else:
                    for line in in_handle:
                        out_handle.write(line)


if __name__ == "__main__":
    main()
