import os

import pandas as pd

structure_dir = 'DALI-benchmark/pdb70'

multiple_models_pdbs = 0
for filename in os.listdir(structure_dir):
    lines = []
    with open(os.path.join(structure_dir, filename), 'r', encoding='utf-8') as pdb_file:
        atom_number = 0
        for line in pdb_file:
             if('ATOM' in line):
                parts = ' '.join(line.split()).split(' ')
                if(int(parts[1]) < atom_number):
                    multiple_models_pdbs += 1
                    break
                else:
                    atom_number = int(parts[1])
             lines.append(line)
    with open(os.path.join(''.join([structure_dir, 'b']), filename), 'w', encoding='utf-8') as pdb_file:
        pdb_file.write(''.join(lines))
print(multiple_models_pdbs)
