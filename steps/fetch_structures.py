from typing import Dict, List, Union, Tuple
import requests
import gzip
import os

from shared import load_config, load_structure_information

def fetch_structures():
    """
    Fetches protein structures from the PDBe database and unzips them.

    The function downloads the assembly files for the specified PDB codes,
    unzips them, and saves them in the specified directory.

    Returns:
        None
    """

    # Load configuration
    config, success, errors = load_config()
    
    #TODO: get list of PDB codes from config

    structures = load_structure_information()

    pdb_codes = list(set([structure['pdb_code'].split('_')[0].lower() for structure in structures]))

    # create directory for raw structures
    tmp_path = config['paths']['tmp']
    raw_structures_path = f"{tmp_path}/structures/raw"
    
    if not os.path.exists(raw_structures_path):
        os.makedirs(raw_structures_path)
        

    # Fetch and unzip structures
    for pdb_code in pdb_codes:
        success = True

        # Construct the URL for the assembly file from PDBe
        pdbe_assembly_url = f"https://www.ebi.ac.uk/pdbe/static/entry/download/{pdb_code}-assembly1.cif.gz"

        # Construct the local path for the downloaded assembly file
        downloaded_assembly_path = f"{raw_structures_path}/{pdb_code}-assembly1.cif.gz"

        # Check if the file already exists
        if not os.path.exists(downloaded_assembly_path):
            # If the file does not exist, download it from PDBe
            print(f"Fetching {pdbe_assembly_url}...")
            response = requests.get(pdbe_assembly_url, stream=True)

            # Check if the request was successful
            if response.status_code == 200:
                with open(downloaded_assembly_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"Downloaded {pdbe_assembly_url} to {downloaded_assembly_path}")
                
            else:
                # If the request failed, set success to False and print an error message
                success = False
                print(f"Failed to fetch {pdbe_assembly_url}. Status code: {response.status_code}")
        
        # Unzip the downloaded assembly file
        if success:
            unzipped_assembly_path = f"{raw_structures_path}/{pdb_code}-assembly1.cif"

            # Check if the unzipped file already exists
            if os.path.exists(unzipped_assembly_path):
                print(f"Unzipped file already exists: {unzipped_assembly_path}")
            else:
                print(f"Unzipping {downloaded_assembly_path}...")
                with gzip.open(downloaded_assembly_path, 'rb') as f_in:
                    with open(unzipped_assembly_path, 'wb') as f_out:
                        f_out.write(f_in.read())
                print(f"Unzipped {downloaded_assembly_path} to {unzipped_assembly_path}")



if __name__ == "__main__":
    fetch_structures()
    print("Structures fetched.")