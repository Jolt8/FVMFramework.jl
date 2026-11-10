"""Export every aliased cell of the FreeCAD parameter table to JSON.

    python export_freecad_table_data.py
    python export_freecad_table_data.py --cad path/to/table_hybrid.FCStd --output path/to/out.json

Opens table_hybrid.FCStd (read only: it is never saved), recomputes it so
formula cells are current, and writes freecad_table_data.json next to this
script. Cells without an alias are ignored. Parameters are listed down each
column, columns left to right, and "easily_viewed_parameters" repeats them
as one "name": "value unit" line each, for searching with Ctrl+F.

Units are kept in FreeCAD's internal units, which are the ones the table is
written in: mm, mm^2, mm^3, kg, s, K and degrees. A quantity's
Unit.Signature gives the exponent of each base unit, so

    33.32 mm        -> {"value": 33.32, "unit": "mm"}
    4 mm^2          -> {"value": 4.0, "unit": "mm^2"}
    30 deg          -> {"value": 30.0, "unit": "°"}
    8 (no unit)     -> {"value": 8, "unit": ""}

The unit strings are written so that Unitful.uparse can read them; see
read_freecad_table_data.jl in this folder.
"""

import argparse
import json
import os
import subprocess
import sys

# FreeCAD 1.1 is built against Python 3.11. If this script is started with
# another interpreter, re-launch it with FreeCAD's bundled python.exe.
freecad_bin = r"C:\Program Files\FreeCAD 1.1\bin"
freecad_lib = r"C:\Program Files\FreeCAD 1.1\lib"
freecad_python = os.path.join(freecad_bin, "python.exe")

if os.path.normcase(sys.executable) != os.path.normcase(freecad_python):
    if not os.path.exists(freecad_python):
        raise RuntimeError(f"FreeCAD's Python executable was not found at '{freecad_python}'. Please ensure FreeCAD 1.1 is installed.")
    relaunch_result = subprocess.run([freecad_python, os.path.abspath(__file__)] + sys.argv[1:])
    sys.exit(relaunch_result.returncode)

if hasattr(os, "add_dll_directory") and os.path.isdir(freecad_bin):
    os.add_dll_directory(freecad_bin)
sys.path.append(freecad_bin)
sys.path.append(freecad_lib)
import FreeCAD as App

# Imported after FreeCAD on purpose: FreeCAD's start-up rebinds the global
# name `datetime` in this script to the datetime class.
from datetime import datetime as DateTime


SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
USER_HOME = os.path.expanduser("~")
DEFAULT_OUTPUT_PATH = os.path.join(SCRIPT_DIRECTORY, "freecad_table_data.json")

TABLE_CANDIDATES = [
    os.path.join(USER_HOME, "Desktop", "hybrid-rocket-CAD", "src", "table_hybrid.FCStd"),
    os.path.join(USER_HOME, "OneDrive", "Desktop", "hybrid-rocket-CAD", "src", "table_hybrid.FCStd"),
]

# FreeCAD's internal unit for each entry of Unit.Signature, written the way
# Unitful.uparse reads it.
FREECAD_INTERNAL_UNIT_NAMES = ("mm", "kg", "s", "A", "K", "mol", "cd", "°")


def internal_unit_string(unit_signature):
    """e.g. (2,0,0,0,0,0,0,0) -> "mm^2", (-1,1,-2,0,0,0,0,0) -> "mm^-1*kg*s^-2"."""
    unit_parts = []
    for base_unit_name, exponent in zip(FREECAD_INTERNAL_UNIT_NAMES, unit_signature):
        if exponent == 0:
            continue
        if exponent == 1:
            unit_parts.append(base_unit_name)
        else:
            unit_parts.append(f"{base_unit_name}^{exponent}")
    return "*".join(unit_parts)


def cell_sort_key(cell_address):
    """"B12" -> (1, "B", 12): down each column, columns left to right (A..Z, then AA..)."""
    column_letters = "".join(character for character in cell_address if character.isalpha())
    row_number = int("".join(character for character in cell_address if character.isdigit()))
    return (len(column_letters), column_letters, row_number)


def main():
    argument_parser = argparse.ArgumentParser(description="Export the aliased cells of the FreeCAD parameter table to JSON.")
    argument_parser.add_argument("--cad", help="FreeCAD file holding the spreadsheet (default: hybrid-rocket-CAD/src/table_hybrid.FCStd)")
    argument_parser.add_argument("--output", default=DEFAULT_OUTPUT_PATH, help="JSON file to write (default: freecad_table_data.json next to this script)")
    arguments = argument_parser.parse_args()

    table_path = arguments.cad
    if table_path is None:
        existing_candidates = [candidate_path for candidate_path in TABLE_CANDIDATES if os.path.exists(candidate_path)]
        if not existing_candidates:
            raise FileNotFoundError("could not find table_hybrid.FCStd (pass --cad); searched:\n  " + "\n  ".join(TABLE_CANDIDATES))
        table_path = existing_candidates[0]
    table_path = os.path.abspath(table_path)

    table_document = App.openDocument(table_path)
    table_document.recompute()  # formula cells (e.g. nozzle_throat_area) are evaluated, not just read

    spreadsheets = [document_object for document_object in table_document.Objects if document_object.TypeId == "Spreadsheet::Sheet"]
    if not spreadsheets:
        raise RuntimeError(f"{table_path} contains no spreadsheet")

    table_entries = {}
    skipped_text_cells = []

    for spreadsheet in spreadsheets:
        aliased_cells = [cell_address for cell_address in spreadsheet.getUsedCells() if spreadsheet.getAlias(cell_address)]

        for cell_address in sorted(aliased_cells, key=cell_sort_key):
            alias = spreadsheet.getAlias(cell_address)
            if alias in table_entries:
                raise RuntimeError(f"alias '{alias}' appears in both '{table_entries[alias]['sheet']}' and '{spreadsheet.Label}'")

            cell_value = spreadsheet.get(cell_address)

            if hasattr(cell_value, "Unit"):
                # Quantity: value and unit both in FreeCAD's internal units.
                entry_value = cell_value.Value
                entry_unit = internal_unit_string(tuple(cell_value.Unit.Signature))
            elif isinstance(cell_value, bool):
                entry_value = int(cell_value)
                entry_unit = ""
            elif isinstance(cell_value, (int, float)):
                entry_value = cell_value
                entry_unit = ""
            else:
                skipped_text_cells.append(f"{alias} ({cell_address}) = {cell_value!r}")
                continue

            table_entries[alias] = {
                "value": entry_value,
                "unit": entry_unit,
                "cell": cell_address,
                "sheet": spreadsheet.Label,
                "contents": spreadsheet.getContents(cell_address),
            }

    # One line per parameter ("name": "33.32 mm"), for finding a name with Ctrl+F.
    easily_viewed_parameters = {}
    for alias, table_entry in table_entries.items():
        easily_viewed_parameters[alias] = f"{table_entry['value']} {table_entry['unit']}".strip()

    table_export = {
        "source_document": table_path.replace("\\", "/"),
        "exported_at": DateTime.now().astimezone().isoformat(timespec="seconds"),
        "unit_note": "values are in FreeCAD's internal units: mm, kg, s, A, K, mol, cd, degrees",
        "easily_viewed_parameters": easily_viewed_parameters,
        "parameters": table_entries,
    }

    output_path = os.path.abspath(arguments.output)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as output_file:
        json.dump(table_export, output_file, indent=1, ensure_ascii=False)
        output_file.write("\n")

    print(f"table:  {table_path}")
    print(f"wrote:  {output_path} ({len(table_entries)} aliased parameters)")
    for skipped_text_cell in skipped_text_cells:
        print(f"skipped text cell: {skipped_text_cell}")


if __name__ == "__main__":
    main()
