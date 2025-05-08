from typing import Dict, List, Tuple, Union
import requests
import os
import csv
import json

from shared import load_config, load_structures_data


def fetch_stcrdab_summary(config:Dict, pdb_code:str, tmp_dir:str) -> Dict:
    """
    Fetch the STCRDAB summary for a given PDB code.

    Args:
        config (Dict): The configuration dictionary.
        pdb_code (str): The PDB code to fetch the summary for.
        tmp_dir (str): The directory to save the summary data.

    Returns:
        Dict: The STCRDAB summary data.
    """

    tmp_file = f"{tmp_dir}/{pdb_code}.txt"
    output_file = f"{config['paths']['stcrdab_information']}/{pdb_code}.json"


    if os.path.exists(tmp_file):
        print (f"STCRDAB summary for {pdb_code} already exists.")
        with open(tmp_file, 'r') as f:
            content = f.read()
            
    else:
        print (f"Fetching STCRDAB summary for PDB code {pdb_code}")
                
        url = f"https://opig.stats.ox.ac.uk/webapps/stcrdab-stcrpred/summary/{pdb_code}"
        response = requests.get(url)
        if response.status_code == 200:
            content = response.text
            if len(content) > 0:
                # Parse the content to extract the summary data
                print (content)
                with open(f"{tmp_dir}/{pdb_code}.txt", 'w') as f:
                    f.write(content)
            else:
                print (f"No summary data found for PDB code {pdb_code}.")

    content_rows = []

    mappings = {
        'Achain': 'alpha',
        'Bchain': 'beta',
        'Gchain': 'gamma',
        'Dchain': 'delta',
    }
    with open(tmp_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            processed_row = {}
            for key in mappings:
                mapped_key = mappings[key]
                if key in row and row[key] != 'NA':
                    processed_row[f"tcr_{mapped_key}"] = row[key]
                    subgroup_label = f"{mapped_key}_subgroup"
                    processed_row[subgroup_label] = row[subgroup_label] 
            content_rows.append(processed_row)

    with open(output_file, 'w') as f:
        json.dump(content_rows, f, indent=4)
    print (content_rows)        
    return content_rows



def fetch_stcrdab_summaries() -> List[Dict]:
    """
    Fetch STCRDAB summaries for all PDB codes in the CSV file.

    Returns:
        List[Dict]: A list of dictionaries containing STCRDAB summary data.
    """
    config, success, errors = load_config()
    if not success:
        print (f"Failed to load configuration: {errors}")
        return []

    tmp_dir = f"{config['paths']['tmp']}/stcrdab"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    

    structures = load_structures_data(parsed=False)

    fetched_pdb_codes = []

    for structure in structures:
        if structure['receptors'] != '' and structure['receptors'] != 'vhh_nanobody':
            if '_' in structure['pdb_code']:
                pdb_code = structure['pdb_code'].split('_')[0].lower()
            else:
                pdb_code = structure['pdb_code'].lower()
            if pdb_code not in fetched_pdb_codes:
                fetched_pdb_codes.append(pdb_code)

                fetch_stcrdab_summary(config, pdb_code, tmp_dir)
 

    return []



if __name__ == "__main__":
    # Fetch STCRDAB summaries
    fetch_stcrdab_summaries()