from typing import Dict, List, Union, Tuple

import os
import json

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import three_to_one

from Levenshtein import distance

from shared import load_config, load_structure_information


class WaterSelect(Select):
    '''
    Select class to select water molecules in a PDB structure.
    '''
    def accept_residue(self, residue):
        if residue.id[0] == 'W':
            return 1
        else:
            return 0

class ResidueSelect(Select):
    '''
    Select class to select polypeptide residues in a PDB structure.
    '''
    def accept_residue(self, residue):
        if residue.id[0] == ' ':
            return 1
        else:
            return 0

class LipidSelect(Select):
    '''
    Select class to select lipid molecules in a PDB structure.
    '''
    def __init__(self, lipid_code):
        self.lipid_code = lipid_code
    def accept_residue(self, residue):
        if residue.id[0].startswith('H_') and residue.get_resname() == self.lipid_code:
            return 1
        else:
            return 0


class AllLipidSelect(Select):
    '''
    Select class to select all named lipid molecules in a PDB structure.
    '''
    def __init__(self, lipid_codes):
        self.lipid_codes = lipid_codes
    def accept_residue(self, residue):
        if residue.id[0].startswith('H_') and residue.get_resname() in self.lipid_codes:
            return 1
        else:
            return 0

class AntigenBindingDomainSelect(Select):
    '''
    Select class to select the antigen binding domain (ABD) in a PDB structure.
    '''
    def __init__(self, chain_id, end_residue_id):
        self.chain_id = chain_id
        self.end_residue_id = end_residue_id
    def accept_residue(self, residue):
        if residue.id[0] == ' ' and residue.id[1] < self.end_residue_id and residue.get_full_id()[2] == self.chain_id:
            return 1
        else:
            return 0



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
    

def create_directories(config: Dict[str, Union[str, List[str]]]) -> None:
    """
    Create directories for storing structure facets

    Args:
        config (Dict[str, Union[str, List[str]]]): Configuration dictionary
    
    Returns:
        None
    """
    structure_facets = config['structure_facets']
    directories = [config['paths'][f"{facet}_structures"] for facet in structure_facets]

    # Create directories for each facet
    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
        # Create sub-directories for each isoform
        for isoform in config['isoforms']:
            isoform_directory = f"{directory}/{isoform}"
            if not os.path.exists(isoform_directory):
                os.makedirs(isoform_directory) 
    pass   



def generate_chain_information(config: Dict[str, Union[str, List[str]]], structure:Structure, pdb_code:str, isoform:str, lipid_names:List[str]) -> Tuple[Dict[str, Union[str, List[str]]], List[str]]:
    chain_info = {}
    errors = []

    b2m_found = False
    cd1_found = False
    cutoff = 20

    cd1_chain_id = None

    for chain in structure.get_chains():
        chain_info[chain.id] = {}
        chain_info[chain.id]['residues'] = [{'aa':three_to_one(residue.get_resname()), 'position':residue.id[1]} for residue in chain.get_residues() if residue.id[0] == ' ' and residue.get_resname() != 'UNK']
        chain_info[chain.id]['waters'] = [{'number': residue.id[1]} for residue in chain.get_residues() if residue.id[0] == 'W']
        chain_info[chain.id]['ligands'] = [residue.get_resname() for residue in chain.get_residues() if residue.id[0].startswith('H_')]
        if len(chain_info[chain.id]['residues']) > 0: 
            chain_info[chain.id]['sequence'] = ''.join([residue['aa'] for residue in chain_info[chain.id]['residues']])
            chain_info[chain.id]['length'] = len(chain_info[chain.id]['residues'])
        else:
            chain_info[chain.id]['sequence'] = ''
            chain_info[chain.id]['length'] = 0
        
        if chain_info[chain.id]['length'] > 0:
            if chain_info[chain.id]['length'] <= 200:
                sequence_distance = distance(chain_info[chain.id]['sequence'], config['sequences']['b2m']['sequence'])
                if sequence_distance < cutoff:
                    chain_info[chain.id]['type'] = 'b2m'
                    b2m_found = True
            else:
                sequence_distance = distance(chain_info[chain.id]['sequence'], config['sequences'][isoform]['sequence'])
                if sequence_distance < cutoff:
                    chain_info[chain.id]['type'] = isoform
                    cd1_found = True
                    cd1_chain_id = chain.id
                else:
                    if config['sequences'][isoform]['sequence'][5:200] in chain_info[chain.id]['sequence']:
                        chain_info[chain.id]['type'] = isoform
                        cd1_found = True
                        cd1_chain_id = chain.id
    if not b2m_found:
        errors.append(f"No b2m chain found in {pdb_code}.")
    
    if not cd1_found:
        errors.append(f"No CD1 chain found in {pdb_code}, should be {isoform}.")
    
    print (f"Cd1 chain id: {cd1_chain_id}")
    # Save the chain information to a JSON file
    if not os.path.exists(generate_facet_folder_path(config, 'information', isoform)):
        os.makedirs(generate_facet_folder_path(config, 'information', isoform))
    filename = f"{generate_facet_folder_path(config, 'information', isoform)}/{pdb_code}_chain_info.json"
    with open(filename, 'w') as f:
        json.dump(chain_info, f, indent=4)

    return chain_info, errors, cd1_chain_id


def process_structure(config: Dict[str, Union[str, List[str]]], pdb_code:str, isoform:str, lipid_names: List[str]) -> List[str]:
    
    # WORK IN PROGRESS

    errors = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, f"output/structures/aligned/{isoform}/{pdb_code}.pdb")

    chain_info, errors, cd1_chain_id = generate_chain_information(config, structure, pdb_code, isoform, lipid_names)

    # Split the structure into separate components
    io = PDBIO()
    io.set_structure(structure)

    # Save the extracellular chains
    filename =  f"{generate_facet_folder_path(config, 'extracellular_domain', isoform)}/{pdb_code}.pdb"
    io.save(filename, select=ResidueSelect())
    
    # Save the water positions
    filename =  f"{generate_facet_folder_path(config, 'water', isoform)}/{pdb_code}.pdb"
    io.save(filename, select=WaterSelect())

    # Save the lipid molecules
    # TODO: maybe add in the ability to save lipids as antigenic or spacer
    for lipid_name in lipid_names:
        filename =  f"{generate_facet_folder_path(config, 'lipid', isoform)}/{pdb_code}_{lipid_name.lower()}.pdb"
        io.save(filename, select=LipidSelect(lipid_name))

    # Save all lipid molecules
    filename =  f"{generate_facet_folder_path(config, 'lipid', isoform)}/{pdb_code}_all.pdb"
    io.save(filename, select=AllLipidSelect(lipid_names))

    #Save the ABDs only
    # TODO: make the abd selection more robust, it will strugggle on single chain constructs
    if cd1_chain_id is not None:
    
        filename =  f"{generate_facet_folder_path(config, 'antigen_binding_domain', isoform)}/{pdb_code}.pdb"
        # get the sequence for the start of the alpha3 linker for the particular isoform
        abd_end = config['abd_breakpoints'][isoform]['alpha3start']
        # locate the residue id for the end of the ABD
        # the end of the ABD is 3 residues after the start of the alpha3 linker
        abd_end_residue_id = chain_info[cd1_chain_id]['sequence'].index(abd_end) + 3



        io.save(filename, select=AntigenBindingDomainSelect(cd1_chain_id, abd_end_residue_id))
    
    #TODO: save the receptors only

    #TODO: make the extracellular domain CD1 selection more robust, removing the receptors

    #TODO: save the hetatm only (no lipids or waters)
    
    return errors






config, success, errors = load_config()


create_directories(config)



raw_structures = load_structure_information()

structures = []

for structure in raw_structures:
    processed_structure = {}
    if not '_' in structure['pdb_code']:
        processed_structure['pdb_code'] = structure['pdb_code'].lower()
        processed_structure['isoform'] = structure['isoform']
        processed_structure['lipid_names'] = []
        if structure['antigenic_pdb_ligand_code'] != '' or structure['antigenic_pdb_lipid_code'] != 'ENDOG':
            processed_structure['lipid_names'].append(structure['antigenic_pdb_ligand_code'])
        if structure['spacer_pdb_ligand_code1'] != '':
            processed_structure['lipid_names'].append(structure['spacer_pdb_ligand_code1'])
        if structure['spacer_pdb_ligand_code2'] != '':
            processed_structure['lipid_names'].append(structure['spacer_pdb_ligand_code2'])
        if structure['receptors'] != '':
            processed_structure['receptors'] = structure['receptors']
        else:
            processed_structure['receptors'] = None

        structures.append(processed_structure)


for structure in structures:
    errors = process_structure(config, structure['pdb_code'], isoform=structure['isoform'], lipid_names=structure['lipid_names'])

    if len(errors) > 0:
        print('')
        print(f"Errors processing structure {structure['pdb_code']}:")
        for error in errors:
            print(error)
        print ('')
    else:
        print(f"Processed structure {structure['pdb_code']} successfully.")

