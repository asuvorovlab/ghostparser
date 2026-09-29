"""Recursively consolidate tabular files with pandas.

The script scans an input directory (including all nested subdirectories),
matches files by start/end name pattern + extension, appends a source-folder
column per row, and writes one consolidated output file.
"""

import argparse
from pathlib import Path

import pandas as pd


def normalize_extension(ext: str) -> str:
    """Normalize an extension to a dotted lowercase form (e.g. '.tsv')."""
    ext = ext.strip().lower()
    if not ext:
        raise ValueError("Extension cannot be empty.")
    if not ext.startswith("."):
        ext = f".{ext}"
    return ext


def infer_delimiter(extension: str) -> str:
    """Infer delimiter from extension with sensible defaults."""
    return "\t" if extension == ".tsv" else ","


def find_matching_files(
    root_dir: Path,
    extension: str,
    file_starts_with: str,
    file_ends_with: str,
    match_subfolder: str,
) -> list[Path]:
    """Find files recursively matching stem start/end and extension."""
    start = file_starts_with.lower()
    end = file_ends_with.lower()
    subfolder_fragment = match_subfolder.strip().strip("/").lower()

    matches: list[Path] = []
    for path in root_dir.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix.lower() != extension:
            continue

        if subfolder_fragment:
            # Match against relative parent path so filtering is scoped under input-folder.
            relative_parent = path.parent.resolve().relative_to(root_dir.resolve()).as_posix().lower()
            wrapped_parent = f"/{relative_parent}/"
            wrapped_fragment = f"/{subfolder_fragment}/"
            if wrapped_fragment not in wrapped_parent:
                continue

        stem = path.stem.lower()
        if start and not stem.startswith(start):
            continue
        if end and not stem.endswith(end):
            continue
        matches.append(path)

    return sorted(matches)


def consolidate_files(
    files: list[Path],
    output_file: Path,
    delimiter: str,
    input_root: Path,
    source_column_name: str,
) -> tuple[int, int]:
    """Concatenate the files into one output file with a source-folder column.

    Args:
        files: Files to read, in the order they are stacked.
        output_file: Path of the consolidated file to write.
        delimiter: Field delimiter of the input and output files.
        input_root: Root the source-folder column is written relative to.
        source_column_name: Name of the source-folder column.

    Returns:
        A tuple ``(files_used, rows_written)``.
    """
    dataframes: list[pd.DataFrame] = []
    expected_columns: list[str] | None = None

    for file_path in files:
        df = pd.read_csv(file_path, sep=delimiter)

        # Skip files with header only and no rows.
        if df.empty:
            continue

        current_columns = df.columns.tolist()
        if expected_columns is None:
            expected_columns = current_columns
        elif current_columns != expected_columns:
            raise ValueError(
                f"Header mismatch in '{file_path}'. "
                f"Expected {expected_columns}, found {current_columns}."
            )

        source_folder = str(file_path.parent.resolve())
        try:
            source_folder = str(file_path.parent.resolve().relative_to(input_root.resolve()))
        except ValueError:
            pass

        df[source_column_name] = source_folder
        dataframes.append(df)

    if not dataframes:
        raise ValueError("No non-empty matching files were found to consolidate.")

    merged = pd.concat(dataframes, ignore_index=True)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_file, sep=delimiter, index=False)

    return len(dataframes), int(len(merged))


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Recursively find files by extension, remove duplicate headers, and "
            "consolidate data rows into one output file."
        )
    )
    parser.add_argument(
        "--input-folder",
        required=True,
        help="Root folder to scan recursively.",
    )
    parser.add_argument(
        "--extension",
        default="tsv",
        help="File extension to match (e.g. tsv, .tsv, csv, .csv).",
    )
    parser.add_argument(
        "--file-starts-with",
        default="",
        help="Only include files whose stem starts with this value.",
    )
    parser.add_argument(
        "--file-ends-with",
        default="",
        help="Only include files whose stem ends with this value.",
    )
    parser.add_argument(
        "--match-subfolder",
        default="",
        help=(
            "Optional subfolder path fragment to restrict search scope "
            "(e.g. 'fdr_corrected'). Only files under matching subfolders are considered."
        ),
    )
    parser.add_argument(
        "--output-folder",
        required=True,
        help="Folder where the consolidated output file will be written.",
    )
    parser.add_argument(
        "--output-file-name",
        default="",
        help=(
            "Optional output file name. Defaults to 'consolidated<extension>' "
            "(e.g. consolidated.tsv)."
        ),
    )
    parser.add_argument(
        "--source-column-name",
        default="source_folder",
        help="Name of the column that records the source folder per row.",
    )
    return parser


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    input_folder = Path(args.input_folder).expanduser().resolve()
    output_folder = Path(args.output_folder).expanduser().resolve()
    extension = normalize_extension(args.extension)

    if not input_folder.exists() or not input_folder.is_dir():
        raise FileNotFoundError(f"Input folder does not exist or is not a directory: {input_folder}")

    matching_files = find_matching_files(
        root_dir=input_folder,
        extension=extension,
        file_starts_with=args.file_starts_with,
        file_ends_with=args.file_ends_with,
        match_subfolder=args.match_subfolder,
    )
    if not matching_files:
        raise FileNotFoundError(
            "No files matched the requested pattern under "
            f"{input_folder}. extension={extension}, "
            f"starts_with='{args.file_starts_with}', ends_with='{args.file_ends_with}', "
            f"match_subfolder='{args.match_subfolder}'"
        )

    output_file_name = args.output_file_name or f"consolidated{extension}"
    output_file = output_folder / output_file_name
    delimiter = infer_delimiter(extension)

    files_used, rows_written = consolidate_files(
        files=matching_files,
        output_file=output_file,
        delimiter=delimiter,
        input_root=input_folder,
        source_column_name=args.source_column_name,
    )

    print(f"Scanned folder: {input_folder}")
    print(
        "Pattern: "
        f"starts_with='{args.file_starts_with}', "
        f"ends_with='{args.file_ends_with}', "
        f"match_subfolder='{args.match_subfolder}', "
        f"extension='{extension}'"
    )
    print(f"Matched files: {len(matching_files)}")
    print(f"Files consolidated: {files_used}")
    print(f"Data rows written: {rows_written}")
    print(f"Source column: {args.source_column_name}")
    print(f"Output file: {output_file}")


if __name__ == "__main__":
    main()
