#!/usr/bin/env python
# -*- coding: utf-8 -*-
# plot-groc.py

# created to plot results from growth_curve_analysis.py
# created from jupyter notebook analyze_gc_results.ipynb

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, argparse, sys, datetime, re

import seaborn as sns
sns.set(context='paper', font_scale=1.3, style='whitegrid')
sns.set_style({'axes.edgecolor': 'black', 'grid.color': 'black'})
plt.rc('text', usetex=False)

def transform(x):
    if x > correction[0]:
        y = np.around(correction[1] * np.exp(correction[2] * x), 3)
    else: y = x
    return y


def reformat_time(x):
    if type(x) not in ['int64', 'float64']:  # if numeric, assume minutes
        if type(x) is datetime.time:
            time_in_seconds = (60. * 60. * x.hour + 60 * x.minute + x.second)
        elif type(x) is datetime.datetime:
            time_in_seconds = (24 * 3600 * x.day + 60. * 60. * x.hour + 60 * x.minute + x.second)
        elif type(x) is str:
            x = datetime.datetime.strptime(x, '%H:%M:%S')
            time_in_seconds = (60. * 60. * x.hour + 60 * x.minute + x.second)
        else:
            print 'Unrecognized time format for ' + str(x)
            return None
        elapsed_time = time_in_seconds / 60.
    else: elapsed_time = x
    time_in_hours = elapsed_time/60
    return time_in_hours


def get_long_form(data, val):
    lf_data = pd.melt(data, id_vars = ['name', 'media', 'expt_date', 'well'], var_name='time', value_name=val)
    lf_data['expt_date'] = lf_data['expt_date'].apply(lambda x: x.strftime('%Y-%m-%d'))
    lf_data['sample'] = lf_data['expt_date'] + ' ' + lf_data['well']
    return lf_data


def import_raw_data(dates, input_dir, blank_file):  # for agg plots
    data_dict = {}
    layout_dict = {}
    for d in dates:
        data_file = os.path.join(input_dir, d, d + ' ypd nacl data.txt')
        data = pd.read_csv(data_file, sep='\t', converters={'strain': lambda x: str(x), 'name': lambda x: str(x)})
        no_time = data.iloc[:, 1:]
        blanked_data = pd.DataFrame(no_time.values - blank_file.values, columns=no_time.columns)
        adj_data = blanked_data.applymap(transform)
        adj_data['Kinetic read'] = data['Kinetic read']
        data_dict[d] = adj_data
        layout_file = os.path.join(input_dir, d, d + ' ypd nacl layout.xlsx')
        layout_dict[d] = pd.read_excel(layout_file, converters={'strain': lambda x: str(x), 'name': lambda x: str(x)})
    # calculate log-transformed, normalized data and compile dataframe
    all_log_norm_data = []
    for d in dates:
        orig_data = data_dict[d]
        orig_data['time'] = orig_data['Kinetic read'].apply(lambda x: reformat_time(x))
        times = orig_data['time']
        only_ODs = orig_data.drop(['Kinetic read', 'time'], axis=1)
        norm_data = only_ODs / only_ODs.iloc[0]
        log_norm_data = np.log(norm_data)
        log_norm_data['time'] = times
        flipped_data = log_norm_data.T
        flipped_data.columns = flipped_data.loc['time']
        flipped_data.drop('time', inplace=True)
        layout = layout_dict[d]
        log_norm_keyed_data = pd.merge(flipped_data, layout, right_on='well', left_index=True)
        all_log_norm_data.append(log_norm_keyed_data)
    log_norm_data_df = pd.concat(all_log_norm_data)
    trimmed_log_norm_data = log_norm_data_df.drop(['run', 'strain', 'clone', 'replicate'], axis=1)
    return trimmed_log_norm_data


def plot_agg_strain(strain, data, val='ln(OD600/init OD)'): # for one strain, all media types
    subset = data[data.name == strain]
    lf_data = get_long_form(subset, val)
    fig = plt.figure()
    ax = sns.tsplot(data=lf_data, time='time', value=val, unit='sample', condition='media', err_style='unit_traces')
    fig.savefig('/Users/nwespe/Desktop/'+strain+'_traces.png', dpi=200, bbox_inches='tight')


first_six_hrs = pd.concat([trimmed_log_norm_data.iloc[:, 0:37], trimmed_log_norm_data.iloc[:, -4:]], axis=1)
# two_six_hrs = pd.concat([trimmed_log_norm_data.iloc[:, 12:37], trimmed_log_norm_data.iloc[:, -4:]], axis=1)
data = first_six_hrs
strains = ['003', '006A', '214', '213', '129', '130', '033B', '188']
group='Hal5 Trk1'


def plot_agg_mult_strains(data, strains, group):    # for multiple strains on one plot, with media split into 2 plots
    fig = plt.figure()
    palette = sns.color_palette()
    media = ['YPD', '0.5 M NaCl']
    media1_data = data[data.media == media[0]]
    media2_data = data[data.media == media[1]]

    med1 = fig.add_subplot(121, adjustable='box')
    x = 0
    for s in strains:
        subset = media1_data[media1_data.name == s]
        lf_data = get_long_form(subset, val)
        ax1 = sns.tsplot(data=lf_data, time='time', value=val, unit='sample', condition='media',
                         color=palette[x], legend=False, ci=95)
        x += 1
    plt.title(media[0])

    med2 = fig.add_subplot(122, adjustable='box')
    x = 0
    for s in strains:
        subset = media2_data[media2_data.name == s]
        lf_data = get_long_form(subset, val)
        ax2 = sns.tsplot(data=lf_data, time='time', value=val, unit='sample', condition='media',
                         color=palette[x], legend=False, ci=95)
        x += 1
    plt.title(media[1])

    med2.legend(strains, title='Strains', bbox_to_anchor=(1, 1), loc=2)
    plt.tight_layout()
    fig.savefig('/Users/nwespe/Desktop/'+group+'_agg.png', dpi=200, bbox_inches='tight')


input = '/Users/nwespe/Dropbox/Research - Data/PG GC Experiments/EffGR 2-6hrs/'
def compile_gc_results(input_dir=input):  # to compile multiple data output files into one results file
    all_results = []
    file_list = os.listdir(input_dir)
    file_list = [f for f in file_list if f[0] != '~']
    for f in file_list:
        f = os.path.join(input_dir, f)
        data = pd.read_excel(f, converters={'strain': lambda x: str(x), 'name': lambda x: str(x),
                                            'expt_date': lambda x: x.strftime('%Y-%m-%d')})
        data['sample'] = data['expt_date'] + ' ' + data['well']
        results = data[['expt_date', 'sample', 'name', 'media', 'growth rate 120-360',
                        'r-squared 120-360', 'doubling time 120-360', 'saturation time',
                        'strain', 'clone']]
        all_results.append(results)
    gc = pd.concat(all_results)
    gc.to_excel('/Users/nwespe/Desktop/all_effgr_2-6hr_results.xlsx')


def plot_together(param=p, gcdata=gc, strain_list = strains, strain_labels=labels, group, pltype='violin'):
    g=sns.factorplot(x='name', y=param + ' 120-360', hue='media', data=gcdata, order=strain_list,
                     kind=pltype, aspect=2, legend=False,
                     hue_order=['YPD', '0.5 M NaCl'], width=1)  # , scale='width'
    sns.despine(offset=20, trim=False)
    g.set_xticklabels(strain_labels, rotation=45)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    plt.title('Growth 2-6 hrs post-NaCl addition')
    if p == 'growth rate':
        g.set(ylim=(0, 0.008), xlabel='Strain', ylabel='Growth rate (1/min)')
        plt.savefig('/Users/nwespe/Desktop/growth plots/' + group + '_gr_' + pltype + '.svg', format='svg',
                    bbox_inches='tight')
    elif p == 'doubling time':
        g.set(ylim=(0, 500), xlabel='Strain', ylabel='Doubling time (min)')  # 0,500
        plt.savefig('/Users/nwespe/Desktop/growth plots/' + group + '_dt_' + pltype + '.svg', format='svg',
                    bbox_inches='tight')