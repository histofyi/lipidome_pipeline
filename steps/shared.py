from typing import Dict, List, Union, Tuple
import csv
import os
import json


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


def load_structure_information() -> List[Dict]:
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
    else:
        structures = []
        print(f"File {information_filepath} does not exist.")
    return structures