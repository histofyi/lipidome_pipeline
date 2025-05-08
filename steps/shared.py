from typing import Dict, List, Union, Tuple
import csv
import os
import json

from Bio.PDB.Structure import Structure
from Bio.PDB.PDBParser import PDBParser


def load_json(file_path:str) -> Tuple[Union[Dict, List], bool, List[str]]:
    """
    Load a JSON file and return its content.
    
    Args:
        file_path (str): Path to the JSON file.

    Returns:
        Tuple[Union[Dict, List], bool, List[str]]: A tuple containing:
            - The content of the JSON file (as a dictionary or list).
            - A boolean indicating success (True) or failure (False).
            - A list of error messages (if any).
    """
    errors = []
    success = True
    if os.path.exists(file_path):
        with open(file_path, 'r') as f:
            data = json.load(f)
    else:
        data = None
        success = False
        errors.append(f"File {file_path} does not exist.")
    return data, success, errors


def load_config() -> Tuple[Dict, bool, List[str]]:
    """
    Load the configuration file and return its content.
    
    Returns:
        Tuple[Dict, bool, List[str]]: A tuple containing:
            - The content of the configuration file (as a dictionary).
            - A boolean indicating success (True) or failure (False).
            - A list of error messages (if any).
    """
    config_path = "config.json"
    config, success, errors = load_json(config_path)
    if not success:
        errors.append(f"Failed to load configuration from {config_path}.")
    return config, success, errors


def generate_facet_folder_path(config: Dict[str, Union[str, List[str]]], facet: str, isoform:str) -> str:
    """
    Generate the folder path for a given facet and isoform

    Args:
        config (Dict[str, Union[str, List[str]]]): Configuration dictionary
        facet (str): The facet for which to generate the folder path
        isoform (str): The isoform for which to generate the folder path
    
    Returns:
        str: The generated folder path
    """
    if facet not in config['structure_facets']:
        pathkey = facet
    else:
        pathkey = f"{facet}_structures"
    base_path = config['paths'][pathkey]
    folder_path = f"{base_path}/{isoform}"
    return folder_path


def load_structures_data(parsed:bool=True) -> List[Dict]:
    """
    Load structure information from a TSV file.

    Returns:
        List[Dict]: A list of dictionaries containing structure information.
    """
    # Load structure information from a TSV file
    information_filepath = 'input/structures.tsv'
    if os.path.exists(information_filepath):
        with open(information_filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            structures = [row for row in reader]
        if parsed:
            structures = [parse_structure_information(structure) for structure in structures if not '_' in structure['pdb_code']]
    else:
        structures = []
        print(f"File {information_filepath} does not exist.")
    return structures


def parse_structure_information(structure: Dict) -> Dict:
    """
    Parse structure information and return relevant details.

    Args:
        structure (Dict): A dictionary containing structure information.

    Returns:
        Dict: A dictionary containing parsed structure information.
    """
    parsed_structure = {}
    parsed_structure['pdb_code'] = structure['pdb_code'].lower()
    parsed_structure['isoform'] = structure['isoform']
    parsed_structure['lipid_names'] = []
    if structure['antigenic_pdb_ligand_code'] != '' or structure['antigenic_pdb_lipid_code'] != 'ENDOG':
        parsed_structure['lipid_names'].append(structure['antigenic_pdb_ligand_code'])
    if structure['spacer_pdb_ligand_code1'] != '':
        parsed_structure['lipid_names'].append(structure['spacer_pdb_ligand_code1'])
    if structure['spacer_pdb_ligand_code2'] != '':
        parsed_structure['lipid_names'].append(structure['spacer_pdb_ligand_code2'])
    if structure['receptors'] != '':
        parsed_structure['receptors'] = structure['receptors']
    else:
        parsed_structure['receptors'] = None
    return parsed_structure


def load_structure(facet:str, isoform:str, pdb_code:str) -> Structure:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, f"output/structures/{facet}/{isoform}/{pdb_code}.pdb")
    if structure:
        return structure
    else:
        return None