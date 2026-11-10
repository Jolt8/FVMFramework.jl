"""Export the ROM geometry JSON from the FreeCAD assembly without opening FreeCAD.

    python export_rom_geometry.py
    python export_rom_geometry.py --cells "*=1" --cells fuel_grain=20 --cells fluid_domain=20
    python export_rom_geometry.py --cad path/to/hybrid_assembly_for_cfd.FCStd --output path/to/out.json

Opens hybrid_assembly_for_cfd.FCStd (read only: it is never saved), recomputes
the ROMModel inside it, validates, and writes rom_geometry/hybrid_rom_geometry.json
next to this script. Validation errors refuse the export (exit code 1) and are
printed; warnings and notes are printed too.

--cells overrides a region's NAxialCells for this export only ("*" = every
region, applied first); the .FCStd file is not changed.

The work is done by the ROMGeometry workbench's scripts/export_headless.py,
which re-launches itself under FreeCAD 1.1's bundled Python. The workbench
repository is looked up in this order: --workbench, the ROM_GEOMETRY_REPO
environment variable, then the usual folders (see WORKBENCH_CANDIDATES).
"""

import argparse
import os
import subprocess
import sys


SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
USER_HOME = os.path.expanduser("~")
APPDATA_DIRECTORY = os.environ.get("APPDATA", "")

DEFAULT_OUTPUT_PATH = os.path.join(SCRIPT_DIRECTORY, "hybrid_rom_geometry.json")

# The working copy of the workbench comes first: an installed copy in the
# FreeCAD Mod folder may be an older snapshot.
WORKBENCH_CANDIDATES = [
    os.path.join(USER_HOME, "Desktop", "ROMwb"),
    os.path.join(USER_HOME, "Desktop", "freecad-rom-geometry"),
    os.path.join(APPDATA_DIRECTORY, "FreeCAD", "v1-1", "Mod", "freecad-rom-geometry"),
    os.path.join(APPDATA_DIRECTORY, "FreeCAD", "v1-1", "Mod", "ROMwb"),
]

CAD_CANDIDATES = [
    os.path.join(USER_HOME, "Desktop", "hybrid-rocket-CAD", "src", "hybrid_assembly_for_cfd.FCStd"),
    os.path.join(USER_HOME, "OneDrive", "Desktop", "hybrid-rocket-CAD", "src", "hybrid_assembly_for_cfd.FCStd"),
]


def first_existing_path(candidate_paths, description):
    for candidate_path in candidate_paths:
        if candidate_path and os.path.exists(candidate_path):
            return candidate_path
    searched_paths = "\n  ".join(path for path in candidate_paths if path)
    raise FileNotFoundError(f"could not find the {description}; searched:\n  {searched_paths}")


def main():
    argument_parser = argparse.ArgumentParser(description="Export the hybrid engine ROM geometry JSON from FreeCAD, headless.")
    argument_parser.add_argument("--cad", help="FreeCAD assembly holding the ROMModel (default: hybrid-rocket-CAD/src/hybrid_assembly_for_cfd.FCStd)")
    argument_parser.add_argument("--output", default=DEFAULT_OUTPUT_PATH, help="JSON file to write (default: hybrid_rom_geometry.json)")
    argument_parser.add_argument("--workbench", help="path of the ROMGeometry workbench repository")
    argument_parser.add_argument("--model", help="label of the ROMModel, if the file has more than one")
    argument_parser.add_argument("--cells", action="append", default=[], metavar="REGION=N", help="override NAxialCells for this export (repeatable, '*' = all regions)")
    arguments = argument_parser.parse_args()

    workbench_candidates = [arguments.workbench, os.environ.get("ROM_GEOMETRY_REPO")] + WORKBENCH_CANDIDATES
    workbench_directory = first_existing_path([path for path in workbench_candidates if path], "ROMGeometry workbench (pass --workbench or set ROM_GEOMETRY_REPO)")
    export_script_path = os.path.join(workbench_directory, "scripts", "export_headless.py")
    if not os.path.exists(export_script_path):
        raise FileNotFoundError(f"{workbench_directory} has no scripts/export_headless.py; is it the ROMGeometry workbench?")

    cad_path = arguments.cad or first_existing_path(CAD_CANDIDATES, "hybrid_assembly_for_cfd.FCStd (pass --cad)")
    output_path = os.path.abspath(arguments.output)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    export_command = [sys.executable, export_script_path, os.path.abspath(cad_path), output_path]
    if arguments.model:
        export_command += ["--model", arguments.model]
    for cell_override in arguments.cells:
        export_command += ["--cells", cell_override]

    print(f"workbench: {workbench_directory}")
    print(f"assembly:  {cad_path}")
    print(f"output:    {output_path}", flush=True)
    export_result = subprocess.run(export_command)
    sys.exit(export_result.returncode)


if __name__ == "__main__":
    main()
