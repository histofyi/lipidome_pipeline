from typing import List, Tuple, Dict, Union

from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio import SeqIO

from collections import Counter
from scipy import stats





input_file = open("input/ebi_mafft_output.fasta", "r")
isoform_fasta_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

canonical = 'CD1a'

sequence_dict = {}

gap_positions = []
conserved_positions = []
most_variable_positions = []

position_offset = 7
max_position = 280

for isoform in isoform_fasta_dict:
    for i, residue in enumerate(isoform_fasta_dict[isoform].seq):
        position = i + position_offset
        
        if position <= max_position:
            if position not in sequence_dict:
                sequence_dict[position] = {'residues': [], 'min_max_entropy': 0.0}
            
            sequence_dict[position]['residues'].append(residue)



max_entropy = 2.0
min_entropy = 0.0

entropy_mapping = {
    0.0: 0.0,
    0.41: 0.2,
    0.5: 0.4,
    0.75: 0.6,
    1.0: 1.0
}

entropy_values = []
for position in sequence_dict:
    residues = sequence_dict[position]['residues']
    if '-' in residues:
        gap_positions.append(position)
        

    residue_counter = Counter(residues)
    entropy = stats.entropy(list(Counter(residues).values()), base=2)

    min_max_entropy = round((entropy - min_entropy) / (max_entropy - min_entropy), 2)
    if min_max_entropy not in entropy_values:
        entropy_values.append(min_max_entropy)

    if min_max_entropy == 0.0:
        conserved_positions.append(position)

    if min_max_entropy == 1.0:
        most_variable_positions.append(position)


    sequence_dict[position]['min_max_entropy'] = entropy_mapping[min_max_entropy]
    
    print (f"Position {position}: {residues} - Entropy: {min_max_entropy}")

print (sequence_dict)
print (conserved_positions)
print (most_variable_positions)
print (gap_positions)
print (sorted(entropy_values))




parser = PDBParser(QUIET=True)
structure = parser.get_structure('isoform_differences', 'input/isoform_difference_structure_7mxf.pdb')


for residue in structure.get_residues():
    position = residue.get_id()[1]
    for atom in residue.get_atoms():
        atom.set_bfactor(sequence_dict[position]['min_max_entropy'])


io = PDBIO()
io.set_structure(structure)
filename =  "output/structures/isoform_differences/test.pdb"
io.save(filename)