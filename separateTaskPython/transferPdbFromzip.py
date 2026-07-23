import os
import zipfile
import shutil
import tempfile

def extract_pdbs_from_nested_zip(source_zip, destination_dir):
    """
    Extract all .pdb files from ZIP files contained inside a source ZIP.

    Parameters
    ----------
    source_zip : str
        Path to the outer ZIP archive.
    destination_dir : str
        Directory where all .pdb files will be copied.
    """

    os.makedirs(destination_dir, exist_ok=True)

    extracted = 0

    with tempfile.TemporaryDirectory() as temp_dir:

        # Extract the outer ZIP
        with zipfile.ZipFile(source_zip, 'r') as outer_zip:
            outer_zip.extractall(temp_dir)

        # Find every inner ZIP
        for root, _, files in os.walk(temp_dir):
            for file in files:
                if not file.lower().endswith(".zip"):
                    continue

                inner_zip_path = os.path.join(root, file)

                with tempfile.TemporaryDirectory() as inner_temp:
                    try:
                        with zipfile.ZipFile(inner_zip_path, 'r') as inner_zip:
                            inner_zip.extractall(inner_temp)

                        # Find every .pdb file inside
                        for r, _, inner_files in os.walk(inner_temp):
                            for f in inner_files:
                                if f.lower().endswith(".pdb"):
                                    src = os.path.join(r, f)

                                    # Prevent filename collisions
                                    base = os.path.splitext(file)[0]
                                    dst = os.path.join(destination_dir,
                                                       f"{base}_{f}")

                                    shutil.copy2(src, dst)
                                    extracted += 1
                                    print(f"Copied: {dst}")

                    except zipfile.BadZipFile:
                        print(f"Skipping invalid ZIP: {inner_zip_path}")

    print(f"\nFinished! Extracted {extracted} PDB files.")

import os
import zipfile

def renamePDBwithProteinName(pdbfolder,protein_name):
    pdb_files = os.listdir(pdbfolder)
    for filename in pdb_files:
        os.rename(f"{pdbfolder}/{filename}",f"{pdbfolder}/{protein_name}_{filename}")
    return 

def extract_pdb_files(zip_file, destination_dir,protein_name=None):
    """
    Extract all .pdb files directly from a ZIP archive.

    Parameters
    ----------
    zip_file : str
        Path to the ZIP file.

    destination_dir : str
        Directory where PDB files will be extracted.

    Returns
    -------
    int
        Number of PDB files extracted.
    """

    os.makedirs(destination_dir, exist_ok=True)

    extracted = 0

    with zipfile.ZipFile(zip_file, "r") as z:

        for member in z.namelist():

            # Ignore directories
            if member.endswith("/"):
                continue

            # Only extract PDB files
            if not member.lower().endswith(".pdb"):
                continue

            # Remove any path information inside the zip
            filename = os.path.basename(member)
            if os.path.exists(os.path.join(destination_dir, filename)):
                continue

            if protein_name:
                filename = f"{protein_name}_{filename}"

            destination = os.path.join(destination_dir, filename)

            
            with z.open(member) as source, open(destination, "wb") as target:
                target.write(source.read())

            extracted += 1
            print(f"Extracted: {filename}")

    if protein_name:
        print("renamed PDB files with protein name")
    print(f"\nFinished! Extracted {extracted} PDB files.")

    return extracted
# source_zip = "/Users/ameyagarwal/Downloads/uL02_AF2_files_bundle.zip"
# destination_dir = "/Users/ameyagarwal/Desktop/Lab/remoteHomologProject/reservePDBfilegroups/uL02_allSeqs_AF2Pred"

# extract_pdbs_from_nested_zip(source_zip, destination_dir)