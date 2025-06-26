from typing import Dict, List, Union, Tuple

import os
import json



from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import three_to_one

from Levenshtein import distance

from shared import load_config, load_structures_data


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






def process_structure(config: Dict[str, Union[str, List[str]]], pdb_code:str, isoform:str, lipid_names: List[str]) -> List[str]:
    
    # WORK IN PROGRESS

    errors = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_code, f"output/structures/aligned/{isoform}/{pdb_code}.pdb")

    # TODO: load in the chain information from the json file
    #chain_info, errors, cd1_chain_id = generate_chain_information(config, structure, pdb_code, isoform, lipid_names)

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
        abd_ends = config['abd_breakpoints'][isoform]['alpha3start']
        # locate the residue id for the end of the ABD
        # the end of the ABD is 3 residues after the start of the alpha3 linker

        abd_end_residue_id = None

        for abd_end in abd_ends:
            if abd_end in chain_info[cd1_chain_id]['sequence']:
                abd_end_residue_id = chain_info[cd1_chain_id]['sequence'].index(abd_end) + 3
                break
        
        if abd_end_residue_id is not None:
            io.save(filename, select=AntigenBindingDomainSelect(cd1_chain_id, abd_end_residue_id))
        else:
            errors.append(f"Could not find the end of the ABD for {pdb_code}.")
    
    #TODO: save the receptors only

    #TODO: make the extracellular domain CD1 selection more robust, removing the receptors

    #TODO: save the hetatm only (no lipids or waters)
    
    return errors






config, success, errors = load_config()


create_directories(config)



structures = load_structures_data()


for structure in structures:
    errors = process_structure(config, structure['pdb_code'], isoform=structure['isoform'], lipid_names=structure['lipid_names'])

    if len(errors) > 0:

        print(f"\nErrors processing structure {structure['pdb_code']}:")
        for error in errors:
            print(error)
        print ('')
    else:
        print(f"Processed structure {structure['pdb_code']} successfully.\n")

