#!/usr/bin/env python
# -*- coding: utf-8 -*-
# OD_growth_finder.py

import pandas as pd
import numpy as np
from scipy import interpolate, signal, stats
import bisect
import matplotlib.pyplot as plt
import datetime as datetime
import os, argparse, sys

import seaborn as sns
sns.set_context('poster', font_scale=1.25)
sns.set_style('ticks')


def reformat_time(x):

    """Assumes x is a datetime.time object."""
    time_in_seconds = (60.*60.*x.hour + 60*x.minute + x.second)
    return time_in_seconds/60.  # return time in minutes for computation of doubling time in minutes


class OD_growth_experiment(object):
    """Each experiment consists of one data file, one plate_layout file (optional), and some other arguments specific
    to the experiment: blank well(s) or value, organism (determines window size for now), and output directory
    """

    def __init__(self, path_to_data, plate_layout=None, blank=0, organism='yeast', out_dir='./'):
        self.path_to_data = path_to_data
        self.data = pd.read_excel(path_to_data)
        self.data.dropna(inplace=True, axis=1)  # Drop the rows that have NAN's, usually at the end
        self.sample_list = []
        self.plate_layout = plate_layout

        self.organism = organism
        self.out_dir = out_dir
        self.name = os.path.splitext(os.path.basename(path_to_data))[0]

        # Get the times from the data - check for format and reformat if necessary
        self.times = self.data.iloc[:,0].values  # assume first column is time data, returns numpy array
        if self.times.dtype not in ['int64', 'float64']:  # if numeric, assume minutes
            self.times = self.times.astype(datetime.time)
            self.elapsed_time = reformat_time(self.times)
        else: self.elapsed_time = self.times
        self.data.drop(self.data.columns[[0]], inplace=True, axis=1)  # remove time column from data to be analyzed

        #if 'Temperature' in self.data.columns: # TODO: not finding temperature column...need fnmatch?
        #    print "found temperature column"
        #    self.data.drop['Temperature']  # remove Temperature column from data to be analyzed

        # check blank
        if type(blank) == 'str':  # TODO: accommodate input of multiple blank wells, average value at each time point
            self.blank = self.data.mean()[blank]  # if blank input is well, blank is average value of that well
        else: self.blank = float(blank)

    def create_sample_list(self):  # creates list of Sample objects

        for well_str in self.data.columns: # well_str is name of column
            raw_data = self.data.loc[:, well_str].values # returns numpy array
            self.sample_list.append(Sample(self, well_str, raw_data))

        return self.sample_list

    def analyze_sample_data(self, method='smooth_n_slide', sample_plots=False, s=0.05):

        for sample in self.sample_list:
            if method == 'smooth_n_slide':
                sample.smooth_n_slide(s)
            elif method == 'spline':
                sample.spline_max_growth_rate(s)
            elif method == 'sliding_window':
                sample.sliding_window()

            if sample_plots:
                folder_name = self.name + '_plots'
                plot_folder = os.path.join(self.out_dir, folder_name)
                if not os.path.exists(plot_folder): os.makedirs(plot_folder)
                sample.plot_growth_parameters(show=False, save=True, folder=plot_folder) # creates one plot per sample
                # default show/save assumes that saving files is desired if function is called for an entire experiment

        return self.sample_list

    def get_window_size(self):
        # determine a good window size - start with 1-1.5*(wt doubling time in minutes)/(time interval in minutes)
        # determine time interval from self.elapsed_time
        interval = self.elapsed_time[1] - self.elapsed_time[0]
        if self.organism == 'yeast':
            window_size = int(90/interval)
        elif self.organism == 'bacteria':
            window_size = int(1.5*20/interval)
        else:
            window_size = 5
            print "I don't know what organism you are using. Here's some data anyway."
        return window_size, interval

    def output_data(self, save=False):

        self.results = pd.DataFrame(
            [[  sample.name,
                sample.lag_time,
                sample.growth_rate,
                sample.doubling_time,
                sample.time_of_max_rate,
                sample.sat_time,
                sample.max_OD,
                sample.time_of_max_OD]
            for sample in self.sample_list],
            columns=("well", "lag time", "growth rate", "doubling time", "time of max growth rate",
                     "saturation time", "max OD", "time of max OD"))

        if self.plate_layout is not None:
            plate_info = pd.read_excel(self.plate_layout)
            self.results = pd.merge(self.results, plate_info, how='outer', on='well')

        if save:
            output_name = self.name + '_output.xlsx'
            output_file = os.path.join(self.out_dir, output_name)
            self.results.to_excel(output_file)

        return self.results

    def plot_histogram(self, show=True, save=False, metric='doubling time', unit='minutes'):
        # plot histogram of growth rates for entire experiment
        sns.distplot(self.results[metric], kde=False, rug=True)
        plt.xlabel(metric + ' (' + unit +')')
        plt.ylabel('number of samples')

        if save:
            hist_name = self.name +'_'+ metric.replace(' ', '_') +'_histogram.png'
            hist_file = os.path.join(self.out_dir, hist_name)
            plt.savefig(hist_file, dpi=300, bbox_inches='tight')

        if show:
            return plt.show()

    def plot_heatmap(self, show=True, save=False, metric='growth rate'):
        # only works if sample names are well IDs
        self.results['row'] = self.results['well'].apply(lambda x: ord(x[0]) - 65)
        self.results['column'] = self.results['well'].apply(lambda x: int(x[1:]) - 1)

        if self.results['row'].max() > 7 or self.results['column'].max() > 11:
            results_arr = np.empty((16, 24))  # assume 384-well plate
            results_arr.fill(np.nan)
        else:
            results_arr = np.empty((8, 12))  # assume 96-well plate
            results_arr.fill(np.nan)

        results_arr[self.results['row'], self.results['column']] = self.results[metric]

        # only show if not saved?
        data = pd.DataFrame(results_arr, index=['A','B','C','D','E','F','G','H'], columns=range(1,13))
        sns.set_style('darkgrid')
        fig = sns.heatmap(data, cmap='spring_r', linewidths=0.01)
        plt.yticks(rotation=0)
        fig.xaxis.set_ticks_position('top')
        plt.title(metric, y=1.1)

        if save:
            heat_name = self.name +'_'+ metric.replace(' ', '_') +'_heatmap.png'
            heat_file = os.path.join(self.out_dir, heat_name)
            plt.savefig(heat_file, dpi=300, bbox_inches='tight')

        if show:
            return plt.show()


class Sample(object):  # create Sample class for attributes collected for each column of data

    def __init__(self, experiment, name, data):
        self.experiment = experiment # now ref attributes as self.experiment.attr
        self.elapsed_time = self.experiment.elapsed_time
        self.out_dir = self.experiment.out_dir
        self.name = name
        self.raw_data = data
        self.cal_data = self.raw_data - self.experiment.blank # subtract blank value, either single value or vector
        self.log_data = np.log(self.cal_data) # Log of the calibrated OD
        # TODO: measure OD values of serial dilution, remove or underweight very low and very high values if not linear

        # get max OD and time of max_OD
        self.max_OD_index = np.argmax(self.raw_data)
        self.max_OD = self.raw_data[self.max_OD_index]
        self.time_of_max_OD = self.elapsed_time[self.max_OD_index]

        # get min OD for lag time calculation
        self.min_logOD = np.amin(self.log_data)

    def spline_max_growth_rate(self, s):

        interpolator = interpolate.UnivariateSpline(self.elapsed_time, self.log_data, k=4, s=s) #k can be 3-5
        der = interpolator.derivative()

        # Get the approximation of the derivative at all points
        der_approx = der(self.elapsed_time)

        # Get the maximum
        self.maximum_index = np.argmax(der_approx)
        self.growth_rate = der_approx[self.maximum_index]
        self.doubling_time = np.log(2)/self.growth_rate
        self.time_of_max_rate = self.elapsed_time[self.maximum_index]

        # Get estimates of lag time and saturation time from 2nd derivative
        der2 = der.derivative()
        der2_approx = der2(self.elapsed_time)
        try: self.lag_index = signal.argrelmax(der2_approx)[0][0]  # find first max
        except: self.lag_index = 0
        if self.lag_index > self.maximum_index: self.lag_index = 0
        self.lag_time = self.elapsed_time[self.lag_index]
        minima = signal.argrelmin(der2_approx)[0]  # find first min after maximum_index
        which_min = bisect.bisect(minima, self.maximum_index)
        try: self.sat_index = minima[which_min]
        except: self.sat_index = len(self.elapsed_time) - 1
        self.sat_time = self.elapsed_time[self.sat_index]

        self.fit_y_values = interpolator(self.elapsed_time)
        self.intercept = self.log_data[self.maximum_index] - (self.growth_rate*self.time_of_max_rate) # b = y - ax

    def sliding_window(self, data_to_use=None):

        if data_to_use is None: data_to_use = self.log_data
        window_size, interval = self.experiment.get_window_size()

        rates = []
        intercepts = []
        r_values = []
        num_windows = len(data_to_use)-window_size+1
        for i in range(0, num_windows):
            window_times = self.elapsed_time[i:i+window_size]
            sub_data = data_to_use[i:i+window_size]
            results = stats.linregress(window_times, sub_data)  # get fit parameters - this takes a while...
            rates.append(results[0])
            intercepts.append(results[1])
            r_values.append(results[2])

        maximum_rate = max(rates)
        # find other slopes at 95% of max rate, use all points to calculate new rate
        max95_rates = [rates.index(i) for i in rates if i >= 0.95*maximum_rate]
        start_point = max95_rates[0]
        end_point = max95_rates[-1] + window_size
        # num_points = len(max95_rates) TODO: include start, end, and number of points in output

        window_times = self.elapsed_time[start_point:end_point]
        sub_data = data_to_use[start_point:end_point]
        results = stats.linregress(window_times, sub_data)
        self.growth_rate = results[0]
        self.intercept = results[1]
        self.r_squared = results[2]**2

        self.maximum_index = int((end_point + start_point)/2) # returns midpoint of time window used to calc rate
        self.time_of_max_rate = self.elapsed_time[self.maximum_index]
        self.doubling_time = np.log(2)/self.growth_rate

        # get lag_time and saturation_time parameters from slope-intercept of max slope
        self.lag_time = (np.amin(data_to_use) - self.intercept) / self.growth_rate # find x where y = min logOD
        self.lag_index = int(self.lag_time/interval)
        if self.lag_index > self.maximum_index: self.lag_index = 0
        self.sat_time = (np.amax(data_to_use) - self.intercept) / self.growth_rate # find x where y = max logOD
        self.sat_index = int(self.sat_time/interval)
        if self.sat_index > len(self.elapsed_time): self.sat_index = len(self.elapsed_time) - 1
        self.fit_y_values = [((self.growth_rate * x) + self.intercept) for x in self.elapsed_time]

    def smooth_n_slide(self, s):

        # first get a spline:
        interpolator = interpolate.UnivariateSpline(self.elapsed_time, self.log_data, k=4, s=s) #k can be 3-5
        self.smooth_log_data = interpolator(self.elapsed_time)

        # now compute rate from sliding window using smoothed data points
        self.sliding_window(data_to_use=self.smooth_log_data)

    def plot_growth_parameters(self, show=True, save=False, folder='./'):
        """Default show/save assumes that showing plot is desired if function is called for one sample
        """
        fit_on_orig = [(np.exp(self.growth_rate * x) * np.exp(self.intercept) + self.experiment.blank) \
                            for x in self.elapsed_time]  # TODO: enable a vector of values for blank

        fig = plt.figure()
        # orig data with points marked for max growth, lag time, saturation time, and max OD
        orig = fig.add_subplot(211)
        orig.plot(self.elapsed_time, self.raw_data, ls='', marker='.', label='raw data')
        orig.plot(self.time_of_max_rate, self.raw_data[self.maximum_index], 'ro', label='max growth')
        orig.plot(self.elapsed_time[self.lag_index], self.raw_data[self.lag_index], 'bo', label='lag time')
        orig.plot(self.elapsed_time[self.sat_index], self.raw_data[self.sat_index], 'go', label='saturation time')
        orig.plot(self.time_of_max_OD, self.raw_data[self.max_OD_index], 'ko', label='max OD')
        orig.plot(self.elapsed_time, fit_on_orig, 'r-', label='fit')

        orig.set_ylabel('OD600')
        orig.set_ylim(0,self.max_OD + 0.05)
        orig.legend(loc='best')

        # log(OD) data with fit line
        logOD = fig.add_subplot(212)
        logOD.plot(self.elapsed_time, self.log_data, ls='', marker='.', label='ln(OD)')
        logOD.autoscale(False) # don't want plot rescaled for fit line
        logOD.plot(self.elapsed_time, self.fit_y_values, 'r-', label='fit')
        logOD.set_xlabel('elapsed time (minutes)')
        logOD.set_ylabel('ln(OD600)')
        logOD.legend(loc='best')

        if show: plt.show(fig)  # TODO: check if this works - don't want plots shown unless show=True

        if save:
            self.plot_file = self.name + '_plot.png'
            plot_path = os.path.join(folder, self.plot_file)
            fig.savefig(plot_path, dpi=200, bbox_inches='tight')
            fig.clf()


def analyze_experiment(
        data_file, plate_layout=None, blank=0, method='smooth_n_slide', organism='yeast',
        sample_plots=False, out_dir='./', s=0.05):
    experiment = OD_growth_experiment(data_file, plate_layout, blank, organism, out_dir)
    print "created experiment"
    experiment.create_sample_list()
    print "input samples"
    experiment.analyze_sample_data(method, sample_plots, s)  # these arguments only apply to analysis
    print "analyzed samples"
    experiment.output_data(save=True)
    print "created output data file"
    return experiment


def make_plots(experiment, show=True, save=False, metric1='doubling time', unit='minutes', metric2='growth rate'):
    experiment.plot_histogram(show, save, metric1, unit)
    experiment.plot_heatmap(show, save, metric2)


def compute_means(experiment, metric=['growth rate'], save=False, keys=None):
    data = experiment.results
    if keys is None:
        print 'Please indicate one or more ways to group samples with key list.'
        sys.exit(1)

    means={}
    for m in metric:
        data_mean = data[m].groupby([data[x] for x in keys]).mean()
        means[m] = data_mean

    data_means = pd.DataFrame(means)

    if save:
        if len(metric) > 1: output_name = experiment.name + '_means.csv'
        else: output_name = experiment.name +'_'+ metric.replace(' ', '_') +'_means.csv'
        output_file = os.path.join(experiment.out_dir, output_name)
        data_means.to_excel(output_file)  # can I save list of

    return data_means


def main():  # Defaults include making all plots and saving all files
    parser = argparse.ArgumentParser(description='Specify a data file to be analyzed.')
    parser.add_argument('-f', '--data_file', required=True, help='Full path to data file.')
    parser.add_argument('-p', '--plate_layout_file', default=None, help='Full path to file with plate layout information.')
    parser.add_argument('-b', '--blank', default=0, help='Either a blank value or well to calculate blank value.')
    parser.add_argument('-m', '--method', default='smooth_n_slide', help='Choose method used for data analysis.')
    parser.add_argument('-o', '--output_directory', default='./', help='Full path to output directory.')
    parser.add_argument('--organism', default='yeast', help='Organism in experiment [yeast, bacteria]; '
                                                            'used in sliding window methods')
    args = parser.parse_args()

    data_file = args.data_file
    plate_layout = args.plate_layout_file # include expt date and run description here
    blank = args.blank
    out_dir = args.output_directory
    method = args.method
    organism = args.organism

    experiment = OD_growth_experiment(data_file, plate_layout, blank, organism, out_dir)
    print "created an experiment"

    experiment.create_sample_list()
    print "input samples"

    experiment.analyze_sample_data(method=method, sample_plots=True)
    print "analyzed samples"

    # create pandas dataframe of sample info and write to output
    results = experiment.output_data(save=True)
    print "created output data file"

    # plots!
    experiment.plot_histogram(results, save=True)
    experiment.plot_heatmap(results, save=True)
    print "saved plots for experiment data"


if __name__ == "__main__":
    main()