#!/usr/bin/env python3
"""Update all the pynucastro networks in the current directory."""

import subprocess
import sys
from pathlib import Path

cwd = Path.cwd()
for net_file in sorted(cwd.glob("**/pynucastro.net")):
    network_dir = net_file.parent

    # find a python file that calls write_network
    for path in network_dir.glob("*.py"):
        with open(path, "r") as f:
            if "write_network" in f.read():
                update_script = path.name
                break
    else:
        # no matching python file found
        continue

    print(f"updating {network_dir.relative_to(cwd)} with {update_script}...")
    result = subprocess.run(
        [sys.executable, update_script],
        cwd=network_dir,
        capture_output=False,
        check=False,
    )
    if result.returncode != 0:
        print(f"error: python exited with status {result.returncode}")
    else:
        print("updated successfully")
    print()
