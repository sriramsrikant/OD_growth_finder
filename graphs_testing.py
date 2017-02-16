__author__ = 'nmcollin'

def histogram(metric, save=False):

    for m in metric:
        sns.barplot(x=keys[0], y=m, hue=keys[1], data=data) #hue_order=['YPD', '0.5 M NaCl'])
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

grouped_data = data[metric].groupby([data[x] for x in keys])
data_calc = grouped_data.agg([np.mean, np.std])

# functions from notebook for data analysis:

def get_clone(row):
    strain = row['strain']
    try: clone = strain[3]
    except: clone = np.nan
    return clone

def strip_clone(row):
    strain = row['strain']
    try:
        clone = strain[3]
        strain = strain[:3]
    except: clone = np.nan
    return strain

# to assign value to a column based on contents (and manipulation) of another column
data003['strain'] = data003.apply(lambda row: strip_clone(row), axis=1 )

# to get subset of dataframe based on column value
data003_IMP2 = data003_genes.loc[data003_genes['gene'] == 'IMP2']

expt1228_means.reset_index()
expt1228_means.columns = [' '.join(col).strip() for col in expt1228_means.columns.values]


pre0218_gr = pre0218.results.loc[:,['strain', 'growth rate', 'media']]
pre0218_gr['media'] = 'pre-incubation'
growthrates0218 = pd.concat([pre0218_gr, expt0218_results], join='inner')
growthrates0218

# barplot
sns.set_context('poster', font_scale=1)
sns.set_style('whitegrid')
sns.barplot(x='strain', y='growth rate', hue='media', data=growthrates0218,
               hue_order=['pre-incubation', 'YPD', '0.5 M NaCl'])
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.title('Growth Rates 02-18-16')
plt.savefig('/Users/nmcollin/Desktop/20160218/all0218_barplot.png')

#add multiple columns to dataframe
expt_df.columns

# create function for plotting quantitative data by any given characteristic
# designed for graphing HU survival data by strain

input_data = pd.read_excel('/Users/nmcollin/Desktop/HU_assay_quantification.xlsx',
                        converters ={'strain': lambda x: str(x), 'pop_or_clone': lambda x: str(x)})
ref_file = pd.read_excel('/Users/nmcollin/Desktop/strain-genotype.xlsx',
                                converters ={'strain': lambda x: str(x)})
input_data_matched = pd.merge(input_data, ref_file, how='left', on='strain', suffixes=('_expt', '_ref'))

key_subset = input_data_matched[input_data_matched[keycol] == key]
subset_expts = pd.unique(data_subset['expt_date'])
subset_expt_rows = input_data_matched[input_data_matched['expt_date'].isin(subset_expts)]
subset_ctrls = subset_expt_rows[subset_expt_rows['strain'] == '006']
subset_evolved = subset_expt_rows[subset_expt_rows['strain'] == '033']
subset_data = pd.concat((key_subset, subset_ctrls, subset_evolved))

gene_data = gene_data[gene_data['media'].isin(['YPD', '0.5 M NaCl'])]
grouped_data = gene_data['fold_death'].groupby([gene_data.strain, gene_data.pop_or_clone, gene_data.media])
data_calc = grouped_data.agg(['count', np.mean, np.std])
data_calc.reset_index(inplace=True)
