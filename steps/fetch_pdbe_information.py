from typing import Dict, List, Tuple, Union
import requests
import os
import csv
import json

from shared import load_config, load_structures_data


def fetch_pdbe_molecules_data(config:Dict, pdb_code:str, tmp_dir:str) -> Dict:
    """
    Fetch the PDBE molecule information for a given PDB code.

    Args:
        config (Dict): The configuration dictionary.
        pdb_code (str): The PDB code to fetch the summary for.
        tmp_dir (str): The directory to save the summary data.

    Returns:
        Dict: The PDBE molecule data.
    """

    tmp_file = f"{tmp_dir}/{pdb_code}.json"
    output_file = f"{config['paths']['pdbe_information']}/{pdb_code}_molecules.json"

    content = None
    if os.path.exists(tmp_file):
        with open(tmp_file, 'r') as f:
            content = json.load(f)
            
    else:
        print (f"Fetching PDBE molecules data for PDB code {pdb_code}")
                
        url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_code}"
        response = requests.get(url)
        if response.status_code == 200:
            content = response.json()
            if len(content) > 0:
                # Parse the content to extract the summary data
                print (content)
                with open(tmp_file, 'w') as f:
                    json.dump(content, f, indent=4)
                
            else:
                print (f"No molecules data found for PDB code {pdb_code}.")

    return content


def parse_protein_data(entry:Dict) -> Dict:
    """
    Parse the PDBE protein data to extract relevant information.

    Args:
        content (Dict): The data on that polypeptide chain.
    

    Returns:
        Dict: The parsed PDBE protein data.
    """
    parsed_data = {}
    
    parsed_data['chains'] = entry['in_chains']
    parsed_data['name'] = entry['molecule_name'][0]
    if entry['gene_name'] is not None:
        parsed_data['gene_name'] = entry['gene_name'][0]
    else:
        parsed_data['gene_name'] = None
    parsed_data['species'] = entry['source'][0]['organism_scientific_name']
    parsed_data['start_residue'] = entry['source'][0]['mappings'][0]['start']['residue_number']
    parsed_data['end_residue'] = entry['source'][0]['mappings'][0]['end']['residue_number']
    parsed_data['sequence'] = entry['sequence']
    parsed_data['pdb_sequence'] = entry['pdb_sequence']

    
    return parsed_data


def parse_bound_molecule_data(entry:Dict) -> Dict:
    """
    Parse the PDBE bound molecules data to extract relevant information.

    Args:
        entry (Dict): The PDBE bound molecule data.

    Returns:
        Dict: The parsed PDBE bound molecule data.
    """

    parsed_data = {}
    
    parsed_data['name'] = entry['molecule_name'][0]
    parsed_data['with_chains'] = entry['in_chains']
    parsed_data['pdb_molecule_code'] = entry['chem_comp_ids'][0]

    return parsed_data


def parse_molecules_data(content:Dict, pdb_code:str) -> Dict:
    """
    Parse the PDBE molecule data to extract relevant information.

    Args:
        content (Dict): The PDBE molecule data.
        pdb_code (str): The PDB code.

    Returns:
        Dict: The parsed PDBE molecule data.
    """

    parsed_data = {
        'protein_chains': [],
        'ligands': [],
        'species': []
    }
    for entry in content[pdb_code]:
        if 'polypeptide' in entry['molecule_type']:
            polypeptide = parse_protein_data(entry)
            parsed_data['protein_chains'].append(polypeptide)
            parsed_data['species'].append(polypeptide['species'])
        elif 'bound' in entry['molecule_type']:
            bound_molecule = parse_bound_molecule_data(entry)
            parsed_data['ligands'].append(bound_molecule)
        elif 'water' in entry['molecule_type']:
            parsed_data['water'] = {'with_chains': entry['in_chains']}

    parsed_data['species'] = list(set(parsed_data['species']))



    if len(parsed_data['species']) > 1:
        # Some structures have VHH nanobodies from Lama glama (llama)
        if not 'Lama glama' in parsed_data['species']:
            print (f"Multiple species found for PDB code {pdb_code}: {parsed_data['species']}")

            print (json.dumps(parsed_data, indent=4))
            print ('')
    elif 'Homo sapiens' not in parsed_data['species']:
            print (f"Non-human structure found for  {pdb_code}: {parsed_data['species']}")

            print (json.dumps(parsed_data, indent=4))
            print ('')
    return parsed_data





def fetch_pdbe_molecules() -> List[Dict]:
    """
    Fetch PDBE molecules data for all PDB codes in the CSV file.

    Returns:
        List[Dict]: A list of dictionaries containing PDBE molecules data.
    """
    config, success, errors = load_config()
    if not success:
        print (f"Failed to load configuration: {errors}")
        return []

    tmp_dir = f"{config['paths']['tmp']}/pdbe"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    

    structures = load_structures_data(parsed=False)

    fetched_pdb_codes = []

    for structure in structures:
        pdb_code = structure['pdb_code'].lower()
        if '_' not in pdb_code:
            fetched_pdb_codes.append(pdb_code)
            content = fetch_pdbe_molecules_data(config, pdb_code, tmp_dir)
            if content:
                parsed_data = parse_molecules_data(content, pdb_code)
                with open(f"{config['paths']['pdbe_information']}/{pdb_code}_molecules.json", 'w') as f:
                    json.dump(parsed_data, f, indent=4)
            else:
                print (f"No molecules data found for PDB code {pdb_code}.")

    return []



if __name__ == "__main__":
    # Fetch STCRDAB summaries
    fetch_pdbe_molecules()