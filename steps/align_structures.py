from typing import Dict, List, Union, Tuple
import os


from pymol import cmd

from shared import load_config, load_structures_data


# Constants for labels
# These labels are used to identify the loaded structures in PyMOL
canonical_label = 'canonical'
mobile_label = 'mobile_structure'


def align_structure(config:Dict, canonical_structure:str, mobile_structure:str, canonical:bool=False, isoform:str=None) -> Dict:
    '''
    Aligns a mobile structure to a canonical structure using PyMOL.

    Args:
        config (Dict): Configuration dictionary containing paths and other settings.
        canonical_structure (str): The name of the canonical structure to align to.
        mobile_structure (str): The name of the mobile structure to be aligned.
        canonical (bool): If True, saves the aligned structure in the canonical directory. Defaults to False.
    
    Returns:
        Dict: A dictionary containing alignment information such as RMSD and atom count.
    '''
    
    print (f"Aligning {mobile_structure} to {canonical_structure}")
    
    # Load the canonical and mobile structures
    canonical_structure_path = config['paths']['canonical_structures']
    mobile_structure_path = config['paths']['raw_structures']

    canonical_structure_file = f"{canonical_structure_path}/{canonical_structure}.pdb"
    mobile_structure_file = f"{mobile_structure_path}/{mobile_structure}-assembly1.cif"

    cmd.load(canonical_structure_file, canonical_label, quiet=1)
    cmd.load(mobile_structure_file, mobile_label, quiet=1)

    # Align the mobile structure to the canonical structure
    align = cmd.align(mobile_label, canonical_label)

    # Remove the canonical structure from the PyMOL session
    # This is done to avoid cluttering the session with multiple loaded structures
    # and to ensure that only the aligned mobile structure is saved
    cmd.delete(canonical_label)

    # Save the aligned mobile structure
    if canonical:
        output_directory = config['paths']['canonical_structures']
        output_formats = ['pdb']
    else:
    
        output_directory = config['paths']['aligned_structures']
        if isoform:
            output_directory = f"{output_directory}/{isoform}"
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
        output_formats = ['pdb', 'cif']

    for file_type in output_formats:
        file_name = f"{output_directory}/{mobile_structure}.{file_type}"
        cmd.save(file_name)

    # Remove the mobile structure from the PyMOL session
    cmd.delete('all')

    # Extract alignment information
    # The alignment information includes RMSD, atom count, cycle count, starting RMSD,
    # starting atom count, match alignment score, and aligned residue count
    alignment_data = {
        'alignment_information':dict(zip(['rmsd','atom_count','cycle_count','starting_rmsd','starting_atom_count','match_align_score','aligned_residue_count'], list(align))),
    }

    
    print (f"Alignment complete. RMSD: {alignment_data['alignment_information']['rmsd']}")
    return alignment_data


def align_canonical_structures(config: Dict[str, Union[str, List[str]]]) -> None:
    canonical_structure_type = 'alignment'
    
    canonical_structure = config['canonical_structures'][canonical_structure_type]

    mobile_structures = [config['canonical_structures'][structure] for structure in config['canonical_structures'] if structure != canonical_structure_type]

    canonical_structure_path = config['paths']['canonical_structures']

    for structure in mobile_structures:
        structure_file_path = f"{canonical_structure_path}/{structure}.pdb"
        if not os.path.exists(structure_file_path):
            alignment_data = align_structure(config, canonical_structure, structure, canonical=True)
        else:
            print(f"Aligned isoform structure {structure} already exists at {structure_file_path}. Skipping alignment.")
            continue
    pass


def align_isoform_structures(config: Dict[str, Union[str, List[str]]]) -> None:

    structures = load_structures_data(parsed=False)
    isoforms = config['isoforms']

    isoform_structures = {}

    for isoform in isoforms:
        isoform_structures[isoform] = list(set([structure['pdb_code'].split('_')[0].lower() for structure in structures if structure['isoform'] == isoform]))

        canonicical_structure = config['canonical_structures'][isoform]
        mobile_structures = isoform_structures[isoform]

        for structure in mobile_structures:
            print(f"Aligning isoform structure {structure} to canonical structure {canonicical_structure}")
            structure_file_path = f"{config['paths']['aligned_structures']}/{structure}.pdb"
            if not os.path.exists(structure_file_path):
                alignment_data = align_structure(config, canonicical_structure, structure, isoform=isoform)
            else:
                print(f"Aligned isoform structure {structure} already exists at {structure_file_path}. Skipping alignment.")
                continue


    pass




if __name__ == "__main__":
    # Load the configuration
    config, success, errors = load_config()
    if not success:
        print("Error loading configuration:", errors)
        exit(1)

    # Align structures
    align_canonical_structures(config)
    align_isoform_structures(config)