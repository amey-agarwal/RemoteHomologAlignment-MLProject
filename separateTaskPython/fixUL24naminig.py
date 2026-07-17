from pathlib import Path
import shutil

def recover_zip_files(root_dir, separator="[forward_slash_character]"):
    """
    separator:
        "／" : Full-width slash (recommended)
        "∕" : Unicode division slash
        "_" : Underscore
    """
    root = Path(root_dir)

    for folder in root.iterdir():
        if not folder.is_dir():
            continue

        # Find zip files directly inside this folder
        zip_files = list(folder.glob("*.zip"))

        if len(zip_files) == 0:
            print(f"No zip found in {folder.name}")
            continue

        if len(zip_files) > 1:
            print(f"Multiple zip files found in {folder.name}, skipping.")
            continue

        zip_path = zip_files[0]

        new_name = f"{folder.name}{separator}{zip_path.name}"
        destination = root / new_name

        print(f"Moving:\n  {zip_path}\n-> {destination}")
        shutil.move(str(zip_path), str(destination))

# Example usage
recover_zip_files("reservePDBfilegroups/rawpredictionzipFiles-af2/localSeqAlignment/uL24_AF2_files_bundle")