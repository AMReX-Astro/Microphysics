#!/usr/bin/env python3

"""Print a release's Markdown notes from CHANGES.md."""

import argparse
from pathlib import Path
import re


def release_notes(changelog, version):
    headings = list(re.finditer(r"^##\s+(\d{2}\.\d{2})\s*$", changelog, re.MULTILINE))
    for index, heading in enumerate(headings):
        if heading.group(1) == version:
            end = headings[index + 1].start() if index + 1 < len(headings) else len(changelog)
            return changelog[heading.end():end].strip()
    raise ValueError(f"No changelog entry for {version}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("version", help="release version (YY.MM)")
    args = parser.parse_args()
    try:
        print(release_notes(Path("CHANGES.md").read_text(), args.version))
    except ValueError as exc:
        parser.error(str(exc))
