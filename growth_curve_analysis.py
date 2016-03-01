#!/usr/bin/env python
# -*- coding: utf-8 -*-
# growth_curve_analysis.py

# adapted from OD_growth_finder.py created by Bryan Weinstein
# renamed to reflect additional calculation methods and distinguish from original program

import pandas as pd
import numpy as np
from scipy import interpolate, signal, stats
import bisect
import matplotlib.pyplot as plt
import os, argparse, sys, datetime, re

import seaborn as sns
sns.set_context('poster', font_scale=1.25)
sns.set_style('darkgrid')


def reformat_time(x):

    """Assumes x is a datetime.time object."""
    time_in_seconds = (60.*60.*x.hour + 60*x.minute + x.second)
    return time_in_seconds/60.  # return time in minutes for computation of doubling time in minutes


class OD_growth_experiment(object):
    """Each experiment consists of one data file, one plate_layout file (optional), and some other arguments specific
    to the experiment: blank well(s) or value, organism (determines window size for now), and output directory
    """

    def __init__(self, path_to_data, plate_layout=None, blank=None, organism='yeast', window_size=None, out_dir='./'):
        self.path_to_data = path_to_data
        self.data = pd.read_excel(path_to_data)
        self.data.dropna(inplace=True, axis=1)  # Drop the rows that have NAN's, usually at the end
        self.samples = {}  # make this a dictionary instead
        self.plate_layout = plate_layout
        self.results = None

        self.organism = organism
        self.window_size = window_size
        self.out_dir = out_dir
        self.name = os.path.splitext(os.path.basename(path_to_data))[0]

        # Get the times from the data - check for format and reformat if necessary
        self.times = self.data.iloc[:, 0].values  # assume first column is time data, returns numpy array
        if self.times.dtype not in ['int64', 'float64']:  # if numeric, assume minutes
            self.elapsed_time = [reformat_time(x) for x in self.times.tolist()]
        else: self.elapsed_time = self.times.tolist()
        self.data.drop(self.data.columns[[0]], inplace=True, axis=1)  # remove time column from data to be analyzed

        # remove Temperature column from data
        temp = [c for c in self.data.columns if re.match('Temp.*', c)]
        if temp: self.data.drop(temp, inplace=True, axis=1)

        # check blank TODO: accommodate input of multiple blank wells, average value at each time point
        if type(blank) is str:
            self.blank = self.data.mean()[blank]  # if blank input is well, blank is average value of that well
        elif blank is None:
            self.blank = 0
        else:
            self.blank = blank  # assumes input is a number

        self.create_sample_list()

    def create_sample_list(self):  # creates dictionary of Sample objects

        for well_str in self.data.columns:  # well_str is name of column
            raw_data = self.data.loc[:, well_str].values  # returns numpy array
            self.samples[well_str] = Sample(self, well_str, raw_data)

        return self.samples

    def analyze_sample_data(self, method='sliding_window', sample_plots=False, s=0.05, droplow=False, start=0, end=None):

        # TODO: put common calculations before sample iterator, e.g. get_window_size()
        self.interval = self.elapsed_time[1] - self.elapsed_time[0]
        if self.window_size is None:
            self.window_size = self.get_window_size(self.interval)

        self.method = method
        if method == 'sliding_window':
            for sample in self.samples.itervalues():
                sample.get_max_growth_parameters(droplow=droplow)
                sample.get_lag_sat_parameters()
        elif method == 'smooth_n_slide':
            for sample in self.samples.itervalues():
                sample.smooth_n_slide(s, droplow=droplow)
        elif method == 'spline':
            for sample in self.samples.itervalues():
                sample.spline_max_growth_rate(s, droplow=droplow)
        elif method == 'effective_growth_rate':
            if end is None:
                end = self.elapsed_time[-1]
            if type(start) not in [int, float]:  # if numeric, assume minutes
                start = reformat_time(datetime.datetime.strptime(start, '%H:%M:%S').time())
            if type(end) not in [int, float]:  # if numeric, assume minutes
                end = reformat_time(datetime.datetime.strptime(end, '%H:%M:%S').time())
            self.start = start
            self.end = end
            for sample in self.samples.itervalues():
                sample.effective_growth_rate(start=start, end=end)

        if sample_plots:
            for sample in self.samples.itervalues():
                folder_name = self.name + '_plots'
                plot_folder = os.path.join(self.out_dir, folder_name)
                if not os.path.exists(plot_folder): os.makedirs(plot_folder)
                sample.plot_growth_parameters(show=False, save=True, folder=plot_folder) # creates one plot per sample
                # default show/save assumes that saving files is desired if function is called for an entire experiment

        return self

    def get_window_size(self, interval):
        # determine a good window size - start with 1-1.5*(wt doubling time in minutes)/(time interval in minutes)
        # determine time interval from self.elapsed_time
        if self.organism == 'yeast':
            window_size = int(90/interval)
        elif self.organism == 'bacteria':
            window_size = int(30/interval)
        else:
            window_size = 5
        if window_size > len(self.elapsed_time)/2:  # recalculate window size for small number of data points
            window_size = int(len(self.elapsed_time)/2)  # max window size possible is half the number of data points
        return window_size

    def output_data(self, method=None, save=False):

        method = self.method
        if method in ['sliding_window', 'smooth_n_slide', 'spline']:
            self.results = pd.DataFrame(
                [[  sample.name,
                    sample.growth_rate,
                    sample.doubling_time,
                    sample.time_of_max_rate,
                    sample.lag_time,
                    sample.lag_OD,
                    sample.sat_time,
                    sample.sat_OD,
                    sample.max_OD,
                    sample.time_of_max_OD]
                for sample in self.samples.itervalues()],
                columns=("well", "growth rate", "doubling time", "time of max growth rate", "lag time", "OD at end of lag",
                         "saturation time", "OD at saturation", "max OD", "time of max OD"))

            # only works if sample names are well IDs - add check?
            self.results['row'] = self.results['well'].apply(lambda x: ord(x[0]) - 64)
            self.results['column'] = self.results['well'].apply(lambda x: int(x[1:]))
            self.results.sort_values(['row','column'], inplace=True)

            if self.plate_layout is not None:
                plate_info = pd.read_excel(self.plate_layout)
                self.results = pd.merge(self.results, plate_info, how='inner', on='well', sort=False)

            self.results.set_index(np.arange(1, len(self.results.index)+1), inplace=True)

        elif method == 'effective_growth_rate':

            start_to_end = str(self.start) + '-' + str(self.end)
            eff_data = pd.DataFrame(
                [[sample.name, sample.effective_gr, sample.effective_dt]
                 for sample in self.samples.itervalues()],
                columns=('well', 'growth rate' + start_to_end, 'doubling time' + start_to_end))
            if self.results:
                self.results = pd.merge(self.results, eff_data, how='left', on='well', sort=False)
            else:
                self.results = eff_data
                self.results['row'] = self.results['well'].apply(lambda x: ord(x[0]) - 64)
                self.results['column'] = self.results['well'].apply(lambda x: int(x[1:]))
                self.results.sort_values(['row','column'], inplace=True)
                if self.plate_layout is not None:
                    plate_info = pd.read_excel(self.plate_layout)
                    self.results = pd.merge(self.results, plate_info, how='inner', on='well', sort=False)
                self.results.set_index(np.arange(1, len(self.results.index)+1), inplace=True)
        else:
            return None

        if save:  # this is either the first results file or a new one with more data
            output_name = os.path.join(self.out_dir, (self.name + '_output.xlsx'))
            i = 0
            while os.path.exists(output_name):
                i += 1
                output_name = os.path.join(self.out_dir, (self.name + '_output' + str(i) + '.xlsx'))
            output_file = pd.ExcelWriter(output_name, engine='xlsxwriter', datetime_format='mmm d yyyy')
            self.results.to_excel(output_file)
            output_file.close()

        return self.results

    def plot_histogram(self, show=True, save=False, metric='doubling time', unit='minutes'):
        # plot histogram of growth rates for entire experiment
        sns.distplot(self.results[metric].dropna(), kde=False, rug=True)
        if unit: plt.xlabel(metric + ' (' + unit +')')
        else: plt.xlabel(metric)
        plt.ylabel('number of samples')

        if save:
            hist_name = self.name +'_'+ metric.replace(' ', '_') +'_histogram.png'
            hist_file = os.path.join(self.out_dir, hist_name)
            plt.savefig(hist_file, dpi=300, bbox_inches='tight')

        if show:
            return plt.show()

    def plot_heatmap(self, show=True, save=False, metric='growth rate', unit=None):

        if self.results['row'].max() > 8 or self.results['column'].max() > 12:
            results_arr = np.empty((16, 24))  # assume 384-well plate
            results_arr.fill(np.nan)
            indices=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
            columns=range(1,25)
        else:
            results_arr = np.empty((8, 12))  # assume 96-well plate
            results_arr.fill(np.nan)
            indices=['A','B','C','D','E','F','G','H']
            columns=range(1,13)

        results_arr[(self.results['row']-1), (self.results['column']-1)] = self.results[metric]

        # only show if not saved?
        data = pd.DataFrame(results_arr, index=indices, columns=columns)
        # TODO: optional input for min/max values of heatmap? outliers can greatly distort scale
        fig = sns.heatmap(data, vmin = 0, cmap='spring_r', linewidths=0.01)
        plt.yticks(rotation=0)
        fig.xaxis.set_ticks_position('top')
        if unit: plt.title(metric + ' (' + unit +')', y=1.1)
        else: plt.title(metric, y=1.1)

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

        self.spline = None

    def spline_max_growth_rate(self, s, droplow=False):

        if droplow: data = np.where(self.log_data < -2.3, 'nan', self.log_data)
        else: data = self.log_data
        interpolator = interpolate.UnivariateSpline(self.elapsed_time, data, k=4, s=s)  #k can be 3-5
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
        self.lag_OD = self.raw_data[self.lag_index]
        minima = signal.argrelmin(der2_approx)[0]  # find first min after maximum_index
        which_min = bisect.bisect(minima, self.maximum_index)
        try: self.sat_index = minima[which_min]
        except: self.sat_index = len(self.elapsed_time) - 1
        self.sat_time = self.elapsed_time[self.sat_index]
        self.sat_OD = self.raw_data[self.sat_index]

        self.spline = interpolator(self.elapsed_time)
        self.intercept = self.log_data[self.maximum_index] - (self.growth_rate*self.time_of_max_rate) # b = y - ax
        self.fit_y_values = [((self.growth_rate * x) + self.intercept) for x in self.elapsed_time]  # for plotting

    def sliding_window(self, data=None, masked=False):

        if data is None: data = self.log_data
        window_size = self.experiment.window_size

        rates = []
        intercepts = []
        num_windows = len(data)-window_size+1
        for i in range(0, num_windows):
            window_times = self.elapsed_time[i:i+window_size]
            sub_data = data[i:i+window_size]
            if masked and sub_data.count() < window_size:  # exclude windows with masked values
                results = [np.nan, np.nan]
            else:
                results = stats.mstats.linregress(window_times, sub_data)  # get fit parameters - this takes a while...
            rates.append(results[0])
            intercepts.append(results[1])
        return rates, intercepts, window_size

    def get_max_growth_parameters(self, data=None, droplow=False):  # run sliding window for max growth rate

        if data is None: data = self.log_data
        if droplow:
            masked_data = np.ma.masked_less(np.copy(data), -2.3)
            data = masked_data  # require more than two points after masking?
        rates, intercepts, window_size = self.sliding_window(data=data, masked=droplow)
        maximum_rate = np.nanmax(rates)
        # find other slopes within 5% of max rate, use all points to calculate new rate
        if maximum_rate <= 0. or np.isnan(maximum_rate):
            self.growth_rate, self.intercept, self.maximum_index, self.time_of_max_rate, self.doubling_time, \
            self.fit_y_values = None, None, None, None, None, None
        else:
            max95_rates = [rates.index(i) for i in rates if i >= 0.95*maximum_rate]
            start_point = max95_rates[0]
            try: end_point = max95_rates[-1] + window_size
            except: end_point = -1
            # num_points = len(max95_rates) TODO: include start, end, and number of points in output
            window_times = self.elapsed_time[start_point:end_point]
            sub_data = data[start_point:end_point]
            results = stats.mstats.linregress(window_times, sub_data)
            new_rate = results[0]
            if new_rate <=0 or np.isnan(new_rate):
                self.growth_rate, self.intercept, self.maximum_index, self.time_of_max_rate, self.doubling_time, \
                    self.fit_y_values = None, None, None, None, None, None
            else:
                self.growth_rate = new_rate
                self.intercept = results[1]
                # self.r_squared = results[2]**2
                self.maximum_index = int((end_point + start_point)/2)  # returns midpoint of time window used to calc rate
                self.time_of_max_rate = self.elapsed_time[self.maximum_index]
                self.doubling_time = np.log(2)/self.growth_rate
                self.fit_y_values = [((self.growth_rate * x) + self.intercept) for x in self.elapsed_time]  # for plotting

    def get_lag_sat_parameters(self):  # run sliding window for raw data rates

        log_data = list(self.log_data)
        first_logOD = log_data[0]
        # find points within 0.2 on logOD scale:
        low_ODs = []
        low_times = []
        for i, p in enumerate(log_data):
            if first_logOD-0.05 <= p <= first_logOD+0.05:
                low_ODs.append(p)
                low_times.append(self.elapsed_time[i])
        lag_line = stats.linregress(low_times, low_ODs)
        if self.growth_rate > 0 and lag_line[0] < self.growth_rate/4:
            self.lag_time = (lag_line[1] - self.intercept)/(self.growth_rate - lag_line[0])
            #self.lag_OD = self.growth_rate*self.lag_time + self.intercept
            self.lag_index = int(self.lag_time/self.experiment.interval)
            if self.lag_index < len(self.raw_data):
                self.lag_OD = self.raw_data[self.lag_index]
            else: self.lag_index, self.lag_time, self.lag_OD = None, None, None
        else: self.lag_index, self.lag_time, self.lag_OD = None, None, None

        rates, intercepts, window_size = self.sliding_window(data=self.raw_data)
        last_rate = rates[-1]
        if self.growth_rate > 0 and last_rate < self.growth_rate/4:
            self.sat_index = self.get_flex_point(rates, 'saturation', window_size)
            if self.sat_index > self.maximum_index:
                self.sat_time = self.elapsed_time[self.sat_index]
                self.sat_OD = self.raw_data[self.sat_index]
            else: self.sat_index, self.sat_time, self.sat_OD = None, None, None
        else: self.sat_index, self.sat_time, self.sat_OD = None, None, None

    def get_flex_point(self, rates, point, window_size):  # poor man's 2nd der, find where rate is changing the fastest
        rate_diffs = [(rates[x]-rates[x+1]) for x in (range(0,len(rates)-1))]
        np.array(rate_diffs)
        sum_rate_diffs = np.add(rate_diffs[:-2], rate_diffs[1:-1])
        summed_rate_diffs = np.add(sum_rate_diffs, rate_diffs[2:])  # locally sum to amplify signal
        if point == 'lag':
            flex_point = np.argmin(summed_rate_diffs) + int(window_size/2) + 1
        elif point == 'saturation':
            flex_point = np.argmax(summed_rate_diffs) + int(window_size/2) + 1
        else: flex_point = None
        return flex_point

    def smooth_n_slide(self, s, droplow):

        if droplow: data = np.where(self.log_data < -2.3, 'nan', self.log_data)
        else: data = self.log_data
        # first get a spline:
        interpolator = interpolate.UnivariateSpline(self.elapsed_time, data, k=4, s=s) #k can be 3-5
        self.smooth_log_data = interpolator(self.elapsed_time)

        # now compute rate from sliding window using smoothed data points
        self.get_max_growth_parameters(self.smooth_log_data, droplow)
        self.get_lag_sat_parameters()

    def effective_growth_rate(self, start=0, end=None):  # calculate for passed sample
        if end is None:
            end = self.elapsed_time[-1]
        if type(start) not in [int, float]:  # if numeric, assume minutes
            start = reformat_time(datetime.datetime.strptime(start, '%H:%M:%S').time())
        if type(end) not in [int, float]:  # if numeric, assume minutes
            end = reformat_time(datetime.datetime.strptime(end, '%H:%M:%S').time())

        # get indices of start and end times in self.elapsed_time list:
        start_index = self.elapsed_time.index(start)
        end_index = self.elapsed_time.index(end)
        time_change = end - start

        # now get ODs and calculate for each sample
        start_logOD = self.log_data[start_index]
        end_logOD = self.log_data[end_index]
        logOD_change = end_logOD - start_logOD
        self.effective_gr = logOD_change/time_change
        self.effective_dt = np.log(2)/self.effective_gr

    def plot_growth_parameters(self, show=True, save=False, folder='./'):
        #  Default show/save assumes that showing plot is desired if function is called for one sample

        fig = plt.figure()
        # orig data with points marked for max growth, lag time, saturation time, and max OD
        orig = fig.add_subplot(211)
        orig.plot(self.elapsed_time, self.raw_data, ls='', marker='.', label='raw data')
        orig.autoscale(False)
        if self.growth_rate is not None:
            fit_on_orig = [(np.exp(self.growth_rate * x) * np.exp(self.intercept) + self.experiment.blank) for x in
                           self.elapsed_time]  # TODO: enable a vector of values for blank
            orig.plot(self.elapsed_time, fit_on_orig, 'r-', label='fit')
            orig.plot(self.time_of_max_rate, self.raw_data[self.maximum_index], 'ro', label='max growth')
        if self.lag_index is not None:
            orig.plot(self.lag_time, self.lag_OD, 'bo', label='lag time')
        if self.sat_index is not None:
            orig.plot(self.sat_time, self.sat_OD, 'go', label='saturation time')
        orig.plot(self.time_of_max_OD, self.max_OD, 'ko', label='max OD')

        orig.set_ylabel('OD600')
        # orig.set_ylim(0,self.max_OD + 0.05)
        orig.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # log(OD) data with fit line
        logOD = fig.add_subplot(212)
        logOD.plot(self.elapsed_time, self.log_data, ls='', marker='.', label='ln(OD)')
        logOD.autoscale(False)  # don't want plot rescaled for fit line
        if self.spline is not None:
            logOD.plot(self.elapsed_time, self.spline, 'k-', label='spline')
        if self.growth_rate is not None:
            logOD.plot(self.elapsed_time, self.fit_y_values, 'r-', label='fit')
        logOD.set_xlabel('elapsed time (minutes)')
        logOD.set_ylabel('ln(OD600)')
        logOD.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if show: plt.show(fig)  # TODO: check if this works - don't want plots shown unless show=True

        if save:
            self.plot_file = self.name + '_plot.png'
            plot_path = os.path.join(folder, self.plot_file)
            fig.savefig(plot_path, dpi=200, bbox_inches='tight')
            fig.clf()
            plt.close()


def analyze_experiment(
        data_file, plate_layout=None, blank=0, method='sliding_window', organism='yeast',
        sample_plots=False, out_dir='./', s=0.1, droplow=False):
    experiment = OD_growth_experiment(data_file, plate_layout, blank, organism, out_dir)
    print "initialized experiment"
    experiment.analyze_sample_data(method, sample_plots, s, droplow)  # these arguments only apply to analysis
    print "analyzed samples"
    experiment.output_data(save=True)
    print "created output data table"
    return experiment


def make_plots(experiment, show=True, save=False, metric1='doubling time', unit1='minutes', metric2='growth rate'):
    experiment.plot_histogram(show, save, metric1, unit1)
    experiment.plot_heatmap(show, save, metric2)


def compute_means(experiment, metric=['growth rate'], save=False, keys=None):  # TODO: add other calculations
    data = experiment.results
    if keys is None:
        print 'Please indicate one or more ways to group samples with key list.'
        sys.exit(1)

    grouped_data = data[metric].groupby([data[x] for x in keys])
    data_calc = grouped_data.agg([np.mean, np.std])

    if save:
        if len(metric) > 1: output_name = experiment.name + '_means.xlsx'
        else: output_name = experiment.name +'_'+ metric[0].replace(' ', '_') +'_means.xlsx'
        output_file = os.path.join(experiment.out_dir, output_name)
        data_calc.to_excel(output_file)

    return data_calc

# TODO: rewrite main function to reflect all inputs into above functions
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
