from typing import Dict, List, Union

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import three_to_one

import os
from shared import load_config


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


def process_structure(config, pdb_code, isoform, lipid_names):
    
    # WORK IN PROGRESS

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('1onq', 'output/structures/aligned/cd1a/1onq.pdb')


    chain_info = {}

    for chain in structure.get_chains():
        chain_info[chain.id] = {}
        chain_info[chain.id]['residues'] = [{'aa':three_to_one(residue.get_resname()), 'position':residue.id[1]} for residue in chain.get_residues() if residue.id[0] == ' ']
        chain_info[chain.id]['waters'] = [{'number': residue.id[1]} for residue in chain.get_residues() if residue.id[0] == 'W']
        chain_info[chain.id]['ligands'] = [residue.get_resname() for residue in chain.get_residues() if residue.id[0].startswith('H_')]


    io = PDBIO()
    io.set_structure(structure)

    # Save the extracellular chains
    filename =  f"{generate_facet_folder_path(config, 'extracellular_domain', isoform)}/{pdb_code}.pdb"
    io.save(filename, select=ResidueSelect())
    
    # Save the water positions
    filename =  f"{generate_facet_folder_path(config, 'water', isoform)}/{pdb_code}.pdb"
    io.save(filename, select=WaterSelect())

    # Save the lipid molecules
    # TODO: make this iterative and driven by the lipids in the structure file
    lipid_code = 'SLF'
    filename =  f"{generate_facet_folder_path(config, 'lipid', isoform)}/{pdb_code}_{lipid_code.lower()}.pdb"
    io.save(filename, select=LipidSelect('SLF'))

    # Save all lipid molecules
    lipid_codes = ['SLF']
    filename =  f"{generate_facet_folder_path(config, 'lipid', isoform)}/{pdb_code}_all.pdb"
    io.save(filename, select=AllLipidSelect(lipid_codes))

    #TODO: save the ABDs only
    #TODO: save the receptor only
    #TODO: make the extracellular domain selection more robust, removing the receptors
    #TODO: save the hetatm only (no lipids or waters)
    







config, success, errors = load_config()


create_directories(config)


structures = ['1onq']

for structure in structures:
    process_structure(config, structure, isoform='cd1a', lipid_names=['SLF'])
    print (f"Processed structure {structure}")

