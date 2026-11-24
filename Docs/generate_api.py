"""Generate Breathe pages with unambiguous paths for alternative implementations."""

import argparse
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET


def generate_api(xml_dir, output_dir):
    # index.xml contains only the basename, even with FULL_PATH_NAMES enabled.
    # Use the source location in each compound XML to disambiguate file pages.
    repo_root = Path(__file__).resolve().parent.parent
    index = ET.parse(xml_dir / "index.xml")
    for compound in index.getroot():
        if compound.get("kind") != "file":
            continue
        refid = compound.get("refid")
        definition = ET.parse(xml_dir / f"{refid}.xml").getroot().find("compounddef")
        location = Path(definition.find("location").get("file"))
        # STRIP_FROM_PATH makes relative locations repository-relative.
        # Older Doxygen versions can still emit absolute source locations.
        if location.is_absolute():
            location = location.relative_to(repo_root)
        filename = location.as_posix()
        compound.find("name").text = filename

    # Breathe 5 also uses the index names to find files, so normalize the XML
    # index as well as the generated directives. Keep Doxygen's unique refids.
    index.write(xml_dir / "index.xml", encoding="utf-8", xml_declaration=True)
    subprocess.run(
        ["breathe-apidoc", "-f", "-o", str(output_dir), str(xml_dir), "-g", "class,file"],
        check=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xml_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    generate_api(args.xml_dir, args.output_dir)
