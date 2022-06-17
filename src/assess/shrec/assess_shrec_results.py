from collections import defaultdict
import pandas as pd
import os
import io
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.size': 8})
import seaborn as sns

results_path = os.path.join('../../../output_shrec_new')

conformations = defaultdict(list)
with open('SHREC2018_ref.cla', 'r', encoding='utf-8') as groundtruth_file:
    line_index = 0
    current_protein = ''
    for line in groundtruth_file:
        if (line_index > 1):
            parts = line.strip()
            if (len(parts) > 0):
                parts = parts.split(' ')
                if(len(parts) == 1 and len(current_protein) > 0):
                    conformations[current_protein].append(parts[0])
                else:
                    current_protein = parts[0]
        line_index += 1

# with open('shrec.log', 'w', encoding='utf-8') as log_file:
#     log_file.write('\t'.join(['filename', 'found', 'in_top', 'missed', 'median-3D', 'in_top20']))
#     log_file.write('\n')
#     for protein in conformations:
#         if (len(conformations[protein]) > 1):
#             for filename in conformations[protein]:
#                 file_path = os.path.join(results_path, '_'.join([filename, 'A', 'whole']), 'candidates', ''.join([filename, '_A-merged-noenriched_shrec.csv']))
#                 results = pd.read_csv(file_path, sep='\t')
#                 results['filename'] = results['filename'].apply(str)
#                 found = 0
#                 lines_read = 0
#                 found_in_top = 0
#                 for index, row in results.iterrows():
#                     if(index == 20):
#                         found_in_top = (found / 20) * 100
#                     if row['filename'] in conformations[protein]:
#                         found += 1
#                         lines_read = index
#
#                 missed = len(conformations[protein]) - 1 - found
#                 log_file.write('\t'.join([filename, repr(found), repr(lines_read + 1), repr(missed), repr(results['3D'].median()), repr(found_in_top)]))
#                 log_file.write('\n')

logged = pd.read_csv('shrec.log', sep='\t')
print('\nmean values:\n')
for key in ['in_top', 'missed', 'median-3D', 'in_top20']:
    print(key, ': ', repr(logged[key].mean()))

print('\nmedian values:\n')
for key in ['in_top', 'missed', 'median-3D', 'in_top20']:
    print(key, ': ', repr(logged[key].median()))


print('\nstd values:\n')
for key in ['in_top', 'missed', 'median-3D', 'in_top20']:
    print(key, ': ', repr(logged[key].std()))


logged[['found', 'missed']].rename({'found': 'Present', \
        'missed': 'Excluded'}, axis=1).plot(kind='area', colormap= ListedColormap(sns.color_palette('colorblind', 2).as_hex()))
plt.ylabel('Total conformations')
plt.xlabel('Reference proteins')
plt.title('Conformations of reference in the results')
plt.xlim([0, len(logged.index)])
plt.savefig('SHREC-searches.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('SHREC-searches.eps', format='eps', dpi=300)


logged[['median-3D']].rename({'median-3D': 'TM-Score', \
        'missed': 'Excluded'}, axis=1).plot(kind='area', colormap= ListedColormap(sns.color_palette('colorblind', 1).as_hex()))
plt.ylabel('Median 3D similarity (%)')
plt.xlabel('Reference proteins')
plt.title('3D similarity in the results')
plt.xlim([0, len(logged.index)])
plt.savefig('SHREC-searches-3D.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('SHREC-searches-3D.eps', format='eps', dpi=300)



print('\nworst case:\n')

# worst_case = logged[logged['missed'] == logged['missed'].max()]
# print(worst_case)
# print('\n')
#
# protein = 'CaM'
# filename = '1737'
# file_path = os.path.join(results_path, '_'.join([filename, 'A', 'segment']), 'candidates', ''.join([filename, '_A_site0-metrics-merged-notenriched.csv']))
# results = pd.read_csv(file_path, sep='\t')
# results['filename'] = [x.split('_')[0] for x in results['structureId']]
# found = 0
# lines_read = 0
# found_in_top = 0
# for index, row in results.iterrows():
#     if(index == 20):
#         found_in_top = (found / 20) * 100
#     if row['filename'] in conformations[protein]:
#         found += 1
#         lines_read = index
#
# missed = len(conformations[protein]) - found
# result = '\n'.join(['\t'.join(['filename', 'found', 'in_top', 'missed', 'in_top20']),
#       '\t'.join([filename, repr(found), repr(lines_read + 1), repr(missed), repr(found_in_top)])])
#
# result = pd.read_csv((io.StringIO(result)), sep='\t')
# print(result)