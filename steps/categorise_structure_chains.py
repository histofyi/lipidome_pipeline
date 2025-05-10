from typing import Dict, List, Union, Tuple

import os
import json

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from Bio.PDB.Polypeptide import three_to_one

from Levenshtein import distance

from shared import load_config, load_structures_data, load_structure, generate_facet_folder_path


def assign_receptor_chains(config:Dict, chain_list:List, pdb_code:str) -> Dict:
    stcrdab_info_file = f"{generate_facet_folder_path(config, 'stcrdab_information', None)}/{pdb_code}.json"
    stcrdab_info = []
    if os.path.exists(stcrdab_info_file):
        stcrdab_info = json.load(open(stcrdab_info_file))
    
    receptor_chains = {}
    for structure_complex in stcrdab_info:
        chain_keys = [key for key in structure_complex.keys() if key.startswith('tcr')]
        for chain_key in chain_keys:
            chain_id = structure_complex[chain_key]
            if structure_complex[chain_key] in chain_list:
                chain_type = chain_key.split('_')[1]
                receptor_chains[chain_id] = [{'roles': [chain_key], 'molecule': structure_complex[f"{chain_type}_subgroup"]}]

    return receptor_chains


def assign_protein_chain(config:Dict, chain_info:Dict, isoform:str, chain_id:str, pdb_code:str) -> Dict:

    errors = []

    hard_cutoff = 20
    soft_cutoff = 40

    chain_role = []
    
    # First check if the chain is a b2m chain (they're shorter than the CD1 chains)

    chain_sequence = chain_info['sequence']
    chain_length = chain_info['length']

    if chain_length <= 200:
        sequence_distance = distance(chain_sequence, config['sequences']['b2m']['sequence'])
        if sequence_distance <= hard_cutoff:
            chain_role.append('b2m')
        elif sequence_distance <= soft_cutoff:
            print (f"Under soft cutoff for {chain_id} in {pdb_code}. Score: {sequence_distance} for beta2m")
    
    # Next we check if the chain is a CD1 chain
    else:
        sequence_distance = distance(chain_sequence, config['sequences'][isoform]['sequence'])
        if sequence_distance <= hard_cutoff:
            chain_role.append(isoform)
        elif sequence_distance <= soft_cutoff:
            print (f"Under soft cutoff  cutoff for {chain_id} in {pdb_code}. Score: {sequence_distance} for {isoform}")
            if config['sequences'][isoform]['sequence'][10:200] in chain_info['sequence']:
                chain_role.append(isoform)
        else:
            if config['sequences'][isoform]['sequence'][10:200] in chain_info['sequence']:
                chain_role.append(isoform)
            if config['sequences']['b2m']['sequence'][5:90] in chain_sequence:
                chain_role.append('b2m')
                chain_role.append('single_chain_construct')

    return chain_role


def assign_lipid_chain(config:Dict, chain_info:Dict, lipid_names:List[str], lipid_roles:List[str], chain_id:str, pdb_code:str) -> Dict:
    errors = []
    chain_lipids = []
    lipids_found = False
    for ligand in chain_info['ligands']:
        if ligand in lipid_names:
            lipids_found = True
            lipid_role = lipid_roles[ligand]
            chain_lipids.append({
                'roles': [lipid_role],
                'molecule': ligand
            })
    if lipids_found: 
        return chain_lipids
    else:
        return None

def generate_structure_information(config: Dict[str, Union[str, List[str]]], structure:Structure, pdb_code:str, isoform:str, lipid_names:List[str], lipid_roles:List[str]) -> Tuple[Dict[str, Union[str, List[str]]], List[str], str, str]:
    chains_info = {}
    errors = []

    cd1_chain_id = None
    b2m_chain_id = None

    # First we'll generate a list of all the chains in the structure
    chain_list = [chain.id for chain in structure.get_chains()]

    # Next we'll see if there are any receptor chains in the structure, we do this next so we don't have to check every protein chain
    receptor_chains = assign_receptor_chains(config, chain_list, pdb_code)
    
    lipid_chains = {}
    cd1_chains = {}
    for chain in structure.get_chains():
        lipid_types = None
        chain_info = {}
        chain_info['residues'] = [{'aa':three_to_one(residue.get_resname()), 'position':residue.id[1]} for residue in chain.get_residues() if residue.id[0] == ' ' and residue.get_resname() != 'UNK']
        chain_info['waters'] = [{'number': residue.id[1]} for residue in chain.get_residues() if residue.id[0] == 'W']
        chain_info['ligands'] = [residue.get_resname() for residue in chain.get_residues() if residue.id[0].startswith('H_')]
        if len(chain_info['residues']) > 0: 
            chain_info['sequence'] = ''.join([residue['aa'] for residue in chain_info['residues']])
            chain_info['length'] = len(chain_info['residues']) 
        else:
            chain_info['sequence'] = ''
            chain_info['length'] = 0

        if chain.id not in receptor_chains:
            if chain_info['length'] > 0:
                chain_type = assign_protein_chain(config, chain_info, isoform, chain.id, pdb_code)
                cd1_chains[chain.id] = chain_type
                if 'b2m' in chain_type:
                    b2m_chain_id = chain.id
                if isoform in chain_type:
                    cd1_chain_id = chain.id 
                if len(chain_info['ligands']) > 0:
                    lipid_types = assign_lipid_chain(config, chain_info, lipid_names, lipid_roles, chain.id, pdb_code)
            else:   
                lipid_types = assign_lipid_chain(config, chain_info, lipid_names, lipid_roles, chain.id, pdb_code)
            if lipid_types is not None:
                lipid_chains[chain.id] = lipid_types
        else:
            lipid_types = assign_lipid_chain(config, chain_info, lipid_names, lipid_roles, chain.id, pdb_code)
            if lipid_types is not None:
                lipid_chains[chain.id] = lipid_types
        

        chains_info[chain.id] = chain_info

    # Now we need to curate the information we have about the structure

    for chain in cd1_chains:
        if cd1_chains[chain] == [isoform]:
            cd1_chains[chain] = [{'roles': ['cd1'], 'molecule': isoform}]
        elif cd1_chains[chain] == ['b2m']:
            cd1_chains[chain] = [{'roles': ['b2m'], 'molecule': 'b2m'}]
        elif cd1_chains[chain] == []:
            errors.append(f"Unable to assign role for chain {chain} in {pdb_code}/{isoform}")
        elif len(cd1_chains[chain]) > 1:
            cd1_chains[chain] = [{'roles': cd1_chains[chain], 'molecule': isoform}]
            errors.append(f"Multiple roles for chain {chain} found in {pdb_code}/{isoform}: {cd1_chains[chain]}")
            

    structure_info = {
        'cd1': cd1_chains,
        'lipids': lipid_chains,
        'receptors': receptor_chains,
    }



    
    # Save the chain information to a JSON file
    
    filename = f"{generate_facet_folder_path(config, 'information', isoform)}/{pdb_code}_chain_info.json"
    with open(filename, 'w') as f:
        json.dump(chains_info, f, indent=4)

    # Save the structure information to a JSON file
    filename = f"{generate_facet_folder_path(config, 'information', isoform)}/{pdb_code}_structure_info.json"
    with open(filename, 'w') as f:
        json.dump(structure_info, f, indent=4)

    return structure_info, errors



def process_structure(config: Dict[str, Union[str, List[str]]], pdb_code:str, isoform:str, lipid_names: List[str], lipid_roles: List[str]) -> Tuple[Dict[str, Union[str, List[str]]], List[str], str, str]:
    '''
    This function processes the structure by removing hetatoms, water, and relettering chains.
    It also stores information about the chains in a JSON file.
 
    Arguments:
        config (Dict[str, Union[str, List[str]]]): Configuration dictionary
        pdb_code (str): The PDB code of the structure
        isoform (str): The isoform of the structure
        lipid_names (List[str]): List of lipid names

    Returns:
        Tuple[Dict[str, Union[str, List[str]]], List[str], str, str]: A tuple containing the chain information, errors, CD1 chain ID, and B2M chain ID
    '''
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, f"output/structures/aligned/{isoform}/{pdb_code}.pdb")

    structure_info = None
    errors = []
    

    print (f"Processing {pdb_code} for isoform {isoform}")
    structure_info, errors = generate_structure_information(config, structure, pdb_code, isoform, lipid_names, lipid_roles)

    return structure_info, errors



def categorise_structure_chains():
    '''
    This function categorises the chains in the structures by assigning roles to them and whether they need remapping or relettering.

    Arguments:
        None
    
    Returns:
        None
    '''
    log = {}

    # Load the configuration file
    config, errors, success = load_config()

    # Create directories for storing structure information
    for isoform in config['isoforms']:
        if not os.path.exists(generate_facet_folder_path(config, 'information', isoform)):
            os.makedirs(generate_facet_folder_path(config, 'information', isoform))

    structures = load_structures_data()

    mapped_chains = config['chain_letters']


    # Process each structure
    for structure in structures:

        pdb_code = structure['pdb_code'].lower()
        isoform = structure['isoform']
        lipid_names = structure['lipid_names']
        lipid_roles = structure['lipid_roles']

        # Process the structure
        structure_info, errors = process_structure(config, pdb_code, isoform, lipid_names, lipid_roles)
        
        # Look at remapping and relettering the chains/lipids
        
        mappings = {}
        current_roles = {}
        for complex_element in structure_info:
            for chain in structure_info[complex_element]:
                if chain not in current_roles:
                    current_roles[chain] = []
                for role in structure_info[complex_element][chain]:
                    for component in role['roles']:
                        current_roles[chain].append(component)
                        # Check if the current roles match the expected roles
                        
                        if chain != mapped_chains[component]:
                            mappings[component] = {"from": chain, "to": mapped_chains[component]}
        
        print (f"Mappings needed: {json.dumps(mappings, indent=4)} \n")
        
        mapping_file = f"{generate_facet_folder_path(config, 'chain_mappings', None)}/{pdb_code}.json"
        with open(mapping_file, 'w') as f:
            json.dump(mappings, f, indent=4)    
        

        # Add the errors to the log for the step
        if len(errors) > 0:
            log[pdb_code] = errors
        

    print (f"Log: {json.dumps(log, indent=4)}")



if __name__ == "__main__":
    categorise_structure_chains()




# Steps 

# Analysing chains - done

# Storing information about the chains (including previous lettering)

