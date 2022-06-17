import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.size': 8})

logged = pd.read_csv('cath.log', sep='\t')

print('median')
for key in ['c_found_top', 'c_found_all', 'ca_found_top', 'ca_found_all', 'cath_found_top', 'cath_found_all', 'tm_found_top', 'tm_found_all',
            'tm_found_top_cand', 'tm_found_all_cand', ]:
    print(key, ': ', repr(logged[key].median()))

print('mean')
for key in ['c_found_top', 'c_found_all', 'ca_found_top', 'ca_found_all', 'cath_found_top', 'cath_found_all', 'tm_found_top', 'tm_found_all',
            'tm_found_top_cand', 'tm_found_all_cand', ]:
    print(key, ': ', repr(logged[key].mean()))

print('std')
for key in ['c_found_top', 'c_found_all', 'ca_found_top', 'ca_found_all', 'cath_found_top', 'cath_found_all', 'tm_found_top', 'tm_found_all',
            'tm_found_top_cand', 'tm_found_all_cand', ]:
    print(key, ': ', repr(logged[key].std()))

# logged[['filename', 'c_found_all', 'ca_found_all', 'cath_found_all']].rename({'c_found_all': 'Same C family occurrences', \
#         'ca_found_all': 'Same CA family occurrences', 'cath_found_all': 'Same CATH family occurrences'}, axis=1).plot(x='filename', kind='bar', stacked=True, \
#         colormap= ListedColormap(sns.color_palette('colorblind', 3).as_hex()), width=0.9)

logged[['filename', 'cath_found_all']].rename({'cath_found_all': 'Same CATH family detection rate'}, axis=1).plot(x='filename', \
                                            kind='bar', colormap= ListedColormap(sns.color_palette('colorblind', 1).as_hex()), width=0.9)

plt.ylabel('Percentage (%)')
plt.xlabel('Reference proteins')
plt.tick_params(axis='x', labelsize=3)
plt.title('Detection of same CATH family proteins in the results')
plt.legend(loc = "lower right")
plt.savefig('CATH-searches-results.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('CATH-searches-results.eps', format='eps', dpi=300)

logged[['filename', 'tm_found_all']].rename({'tm_found_all': 'TM-Score'}, axis=1).plot(x='filename', kind='bar', width=0.9)
plt.ylabel('Median 3D similarity (%)')
plt.xlabel('Reference proteins')
plt.tick_params(axis='x', labelsize=3)
plt.title('3D similarity of the proteins in the results')
plt.savefig('CATH-searches-median-3D.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('CATH-searches-median-3D.eps', format='eps', dpi=300)


logged[['filename', 'c_found_top', 'ca_found_top', 'cath_found_top']].rename({'c_found_top': 'Same C family occurrences', \
        'ca_found_top': 'Same CA family occurrences', 'cath_found_top': 'Same CATH family occurrences'}, axis=1).plot(x='filename', kind='bar', stacked=True, \
        colormap= ListedColormap(sns.color_palette('colorblind', 3).as_hex()), width=0.9)
plt.yscale('log')
plt.ylabel('Sum of percentages (log scale)')
plt.xlabel('Reference proteins')
plt.tick_params(axis='x', labelsize=3)
plt.legend(loc = "lower right")
plt.title('CATH classification of the proteins in the top 20 results')
plt.savefig('CATH-searches-top-results.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('CATH-searches-top-results.eps', format='eps', dpi=300)

logged[['filename', 'tm_found_top']].rename({'tm_found_top': 'TM-Score'}, axis=1).plot(x='filename', kind='bar', width=0.9)
plt.ylabel('Median 3D similarity (%)')
plt.xlabel('Reference proteins')
plt.tick_params(axis='x', labelsize=3)
plt.title('3D similarity of the proteins in the top 20 results')
plt.savefig('CATH-searches-top-median-3D.png', format='png', dpi=300)
ax = plt.gca()
ax.set_rasterized(True)
plt.savefig('CATH-searches-top-median-3D.eps', format='eps', dpi=300)