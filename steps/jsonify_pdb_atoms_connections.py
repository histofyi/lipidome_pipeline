import json

input_code = "1onq_slf"


with open(f"output/structures/lipids/cd1a/{input_code}.pdb", 'r') as file:
    lines = file.readlines()

output = {
    'atoms': [],
    'connections': [],
}

for line in lines:
    if line.startswith('HETATM'):
        split_line = [item for item in line.split() if item != 'HETATM']
        atom = {
            'number': int(split_line[0]),
            'pdb_atom_label': split_line[1],
            'atom_type': split_line[-1],
        }
        output['atoms'].append(atom)
        print (atom)
    elif line.startswith('CONECT'):
        split_line = [item for item in line.split() if item != 'CONECT']
        output['connections'].append(split_line)
        print (split_line)

with open(f"output/structures/lipid_json/{input_code}_pdb_atom_labels.json", 'w') as file:
    json.dump(output, file, indent=4)

for atom in output['atoms']:
    del atom['pdb_atom_label']

with open(f"output/structures/lipid_json/{input_code}.json", 'w') as file:
    json.dump(output, file, indent=4)