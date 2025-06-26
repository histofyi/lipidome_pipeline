from typing import Dict, List, Union

import anthropic
import json
import csv
import os

with open('secrets.json') as f:
    secrets = json.load(f)

client = anthropic.Anthropic(
    api_key=secrets['claude_api_key'],
)




tmp_dir = 'tmp/claude'
if not os.path.exists(tmp_dir):
    os.makedirs(tmp_dir)


def prompt_claude_for_lipid(lipid:Dict, verbose:bool=False) -> Dict:
    '''
    This function prompts Claude for lipid information from a SMILES string.

    Args:
        lipid (Dict): A dictionary containing lipid information, including 'code', 'molecular_formula', 'molecular_weight', and 'smiles'.
        verbose (bool): If True, prints additional information during processing.
    
    Returns:
        Dict: A dictionary containing the lipid analysis results.
    '''

    tmp_file = f"{tmp_dir}/{lipid['lipid_pdb_code']}.txt"

    if not os.path.exists(tmp_file):


        prompt = "You are an expert biochemist specializing in lipid analysis. Your task is to analyze a lipid molecule based on its SMILES string representation and provide a detailed description along with a structured JSON object containing key information about the lipid.\n\nHere is the SMILES string of the lipid molecule you need to analyze:\n\n<SMILES_STRING>\n{{SMILES_STRING}}\n</SMILES_STRING>\n\nPlease follow these steps to analyze the lipid molecule:\n\n1. Parse the SMILES string and identify the key components of the lipid molecule.\n2. Determine the lipid type and class.\n3. Calculate the molecular formula and molecular weight.\n4. Analyze the headgroup, identifying its type and chemical name.\n5. Analyze each tail, providing the following information:\n   - Type of tail (e.g., fatty acid, sphingoid base)\n   - Name of the tail\n   - Length (number of carbon atoms)\n   - Saturation status\n   - For unsaturated tails:\n     - Position of unsaturations from the terminal end (opposite to the headgroup)\n     - Stereochemistry of the unsaturations\n\nBefore providing your final output, perform a thorough analysis of the lipid molecule within <lipid_analysis> tags. Follow these steps in your analysis:\n\n1. Break down the SMILES string into its components\n2. Identify the headgroup and tails\n3. Analyze each component separately\n4. Determine the lipid type and class based on the analysis\n5. Calculate the molecular formula and weight\n\n\nThis will help ensure a comprehensive interpretation of the molecular structure.\n\nYour final output should consist of three parts:\n\n1. A detailed text description of the lipid molecule, including all the information mentioned above.\n\n2. A JSON object with the following structure:\n\n```json\n{\n  \"lipid_type\": \"\",\n  \"lipid_class\": \"\",\n  \"molecular_formula\": \"\",\n  \"molecular_weight\": \"\",\n  \"molecule_name\":,\n  \"headgroup\": {\n    \"type\": \"\",\n    \"chemical_name\": \"\"\n  },\n  \"tails\": [\n    {\n      \"type\": \"\",\n      \"name\": \"\",\n      \"length\": 0,\n      \"saturation\": \"\",\n      \"unsaturations\": [\n        {\n          \"position_from_terminal_end\": \"\",\n          \"configuration\": \"\"\n        }\n      ]\n    }\n  ]\n}\n```\n\nPlease ensure that your analysis is thorough and accurate, and that both the text description and JSON object contain all the required information."
        prompt = prompt.replace("{{SMILES_STRING}}", lipid['smiles'])
        message = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=20000,
            temperature=0,
            messages=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": "<examples>\n<example>\n<SMILES_STRING>\nCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O[S](O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC\n</SMILES_STRING>\n<ideal_output>\n<lipid_analysis>\n1. Breaking down the SMILES string:\n   - CCCCCCCCCCCCCCCCCC(=O)N: Long carbon chain with an amide group\n   - [C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O[S](O)(=O)=O)[C@H]1O): Sugar moiety with a sulfate group\n   - [C@H](O)\\C=C\\CCCCCCCCCCCCC: Another long carbon chain with a double bond\n\n2. Identifying components:\n   - Headgroup: Sulfated sugar moiety\n   - Tails: Two long carbon chains\n\n3. Analyzing components:\n   - Headgroup: The sugar appears to be a galactose with a sulfate group at position 3\n   - Tail 1: 18-carbon chain with an amide linkage (fatty acid)\n   - Tail 2: 18-carbon chain with a double bond (sphingoid base)\n\n4. Determining lipid type and class:\n   - The presence of a sugar headgroup and two tails (one fatty acid, one sphingoid base) indicates this is a glycosphingolipid\n   - The sulfate group on the sugar classifies this as a sulfatide (cerebroside sulfate)\n\n5. Calculating molecular formula and weight:\n   - Counting atoms: C42H81NO11S\n   - Molecular weight: 808.17 g/mol\n</lipid_analysis>\n\n{\n  \"lipid_type\": \"Glycosphingolipid\",\n  \"lipid_class\": \"Sulfatide (Cerebroside sulfate)\",\n  \"molecular_formula\": \"C42H81NO11S\",\n  \"molecular_weight\": \"808.17\",\n  \"headgroup\": {\n    \"type\": \"Sulfated monosaccharide\",\n    \"chemical_name\": \"Galactose-3-sulfate\"\n  },\n  \"tails\": [\n    {\n      \"type\": \"Fatty acid\",\n      \"name\": \"Stearic acid\",\n      \"length\": 18,\n      \"saturation\": \"Saturated\",\n      \"unsaturations\": [],\n      \"additional_features\": [\"Amide linkage\"]\n    },\n    {\n      \"type\": \"Sphingoid base\",\n      \"name\": \"Sphingosine\",\n      \"length\": 18,\n      \"saturation\": \"Monounsaturated\",\n      \"unsaturations\": [\n        {\n          \"position_from_terminal_end\": \"14-15\",\n          \"configuration\": \"cis\"\n        }\n      ]\n    }\n  ]\n}\n</ideal_output>\n</example>\n</examples>\n\n"
                        },
                        {
                            "type": "text",
                            "text": prompt
                        }
                    ]
                },
                {
                    "role": "assistant",
                    "content": [
                        {
                            "type": "text",
                            "text": "<lipid_analysis>"
                        }
                    ]
                }
            ]
        )

        output = message.content[0].text
        if verbose:
            print("Claude's response:")
            print("========================================")
            print(output)

        with open(tmp_file, 'w') as f:
            f.write(output)

    else:
        with open(tmp_file, 'r') as f:
            output = f.read()


    json_output = output.split('json')[1].replace('```', '').strip()

    print ('')
    print (json_output)

    json_data = json.loads(json_output)


    output_folder = 'output/structures/lipid_descriptions'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    output_file = f"{output_folder}/{lipid['lipid_pdb_code']}.json"
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=4)

    print(f"Output saved to {output_file}")
    print(f"Lipid analysis for {lipid['lipid_pdb_code']} completed successfully.") 

    return json_data



# import the information on the lipids from the tsv file

filename = 'input/structure_lipids.tsv'

with open(filename, 'r') as f:
    reader = csv.DictReader(f, delimiter='\t')
    lipids = [row for row in reader]


for lipid in lipids:
    print(f"Processing lipid: {lipid['lipid_pdb_code']}")

    json_data = prompt_claude_for_lipid(lipid)

