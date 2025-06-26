import os
import json

input_folder = 'output/structures/chain_mappings'

mappings = {}

keys = ['tcr_alpha', 'tcr_beta', 'antigenic_lipid', 'spacer_lipid1', 'spacer_lipid2', 'cd1', 'b2m']

for key in keys:
    mappings[key] = []

structure_count = 0

for filename in os.listdir(input_folder):
    if filename.endswith('.json'):
        structure_count += 1
        with open(os.path.join(input_folder, filename), 'r') as f:
            data = json.load(f)
            print(f"File: {filename}")
            for key in data:
                if key not in mappings:
                    mappings[key] = []
                mappings[key].append(data[key]['from'])
            print()

print (json.dumps(mappings, indent=4))

remapped_mappings = {}

for key in mappings:
    if key == 'tcr_gamma':
        remapped_key = 'tcr_alpha'
    elif key == 'tcr_delta':    
        remapped_key = 'tcr_beta'
    else:
        remapped_key = key
    if remapped_key not in remapped_mappings:
        remapped_mappings[remapped_key] = {'values': [], 'count': 0, 'options':{}}
    remapped_mappings[remapped_key]['values'].extend(mappings[key])
    remapped_mappings[remapped_key]['count'] += len(mappings[key])
    for chain in mappings[key]:
        if chain not in remapped_mappings[remapped_key]['options']:
            remapped_mappings[remapped_key]['options'][chain] = 0
        remapped_mappings[remapped_key]['options'][chain] += 1

print (json.dumps(remapped_mappings, indent=4))

print (f"Total number of structures: {structure_count}")