import json

interaction_file = 'output/structures/interaction_information/1onq_slf.json'

pdb_atom_number_mapping_file = 'output/structures/lipid_json/1onq_slf_pdb_atom_labels.json'

lipid_atom_description_file = 'output/structures/lipid_atom_descriptions/1onq_slf.json'


lipid_atom_interactions = {}

with open(interaction_file, 'r') as file:
    interaction_data = json.load(file)

with open(pdb_atom_number_mapping_file, 'r') as file:   
    pdb_atom_number_mapping = json.load(file)

with open(lipid_atom_description_file, 'r') as file:   
    lipid_atom_descriptions = json.load(file)

atom_mappings = {}
for atom in pdb_atom_number_mapping['atoms']:
    atom_mappings[atom['pdb_atom_label']] = atom['number']


interactions = interaction_data['1onq'][0]['interactions']

print (f"Number of interactions: {len(interactions)}")

hydrophobic_interaction_count = 0
hbond_interaction_count = 0


for interaction in interactions:
    for atom in interaction['ligand_atoms']:
        atom_number = atom_mappings[atom]
        if not atom_number in lipid_atom_interactions:
            lipid_atom_interactions[atom_number] = {'interaction_count': 0, 'partner_count':0, 'partners':{}}
        print (interaction)
        partner_key = f"{interaction['end']['chem_comp_id'].lower()}_{interaction['end']['author_residue_number']}"
        if not partner_key in lipid_atom_interactions[atom_number]['partners']:
            lipid_atom_interactions[atom_number]['partners'][partner_key] = {
                'position': interaction['end']['author_residue_number'],
                'residue_name': interaction['end']['chem_comp_id'],
                'interactions': []
            }
            lipid_atom_interactions[atom_number]['partner_count'] += 1
    
        interaction_details = {
            'type': interaction['interaction_type'].lower(),
            'details': [detail.lower() for detail in interaction['interaction_details']],
            'distance': interaction['distance'],
            'atoms': interaction['end']['atom_names']
        }
        lipid_atom_interactions[atom_number]['partners'][partner_key]['interactions'].append(interaction_details)
        lipid_atom_interactions[atom_number]['interaction_count'] += 1
    
    print (interaction)


print (f"Number of lipid atoms: {len(lipid_atom_interactions)}")

print (json.dumps(lipid_atom_interactions, indent=4))

for atom_group in lipid_atom_descriptions['atom_groups']:
    print (f"Processing atom group: {atom_group}")
    for atom in lipid_atom_descriptions['atom_groups'][atom_group]['atoms']:

        if atom['number'] in lipid_atom_interactions:
            atom['interactions'] = lipid_atom_interactions[atom['number']]
            
        else:
            atom['interactions'] = {
                'interaction_count': 0,
                'partner_count': 0,
                'partners': {}
            }
        print (json.dumps(atom, indent=4))


