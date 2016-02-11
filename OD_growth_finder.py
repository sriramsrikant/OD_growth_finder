#!/usr/bin/env python
# -*- coding: utf-8 -*-
# OD_growth_finder.py

import pandas as pd
import numpy as np
from scipy import interpolate, signal, stats
import bisect
import matplotlib.pyplot as plt
import datetime as datetime
import os, argparse


def reformat_time(x):
    """Assumes x is a datetime.time object."""
    time_in_seconds = (60.*60.*x.hour + 60*x.minute + x.second)
    return time_in_seconds/(60.) # return time in minutes for computation of doubling time in minutes


class OD_growth_experiment(object):

    def __init__(self, path_to_data, blank, method, organism, out_dir, s = 0.05):
        self.path_to_data = path_to_data
        self.data = pd.read_excel(path_to_data) # not sure if read_csv might be better
        # Drop the rows that have NAN's, usually at the end
        self.data.dropna(inplace=True, axis=1)
        self.sample_list = []
        self.method = method
        self.organism = organism
        self.out_dir = out_dir
        self.name = os.path.splitext(os.path.basename(path_to_data))[0]

        # Get the times from the data - check for format and reformat if necessary
        self.times = self.data.iloc[:,0].values # just assume first column is time data, returns numpy array
        if self.times.dtype != 'int64': # if numeric, assumes minutes
            self.times = self.times.astype(datetime.time)
            self.elapsed_time = self.times.apply(reformat_time)
        else: self.elapsed_time = self.times

        self.data.drop(self.data.columns[[0]], inplace=True, axis=1) # remove time column from data to be analyzed

        # check blank
        if type(blank) == 'str':
            self.blank = self.data.mean()[blank] # if blank input is well, blank is average value of that well
        else: self.blank = float(blank)

        # Set the default s for fitting...deals with how close the fit is to the points
        self.s = s

    def get_sample_data(self):
        for well_str in self.data.columns: # well_str is name of column
            # create Sample object Sample(experiment, data) and add to sample list
            raw_data = self.data.loc[:, well_str].values # returns numpy array
            self.sample_list.append(Sample(self, well_str, raw_data))
        return self.sample_list

    def analyze_sample_data(self):

        for sample in self.sample_list:
            if self.method == 'spline':
                sample.spline_max_growth_rate()
            elif self.method == "sliding_window":
                sample.sliding_window(sample.log_data)
            elif self.method == "smooth_n_slide":
                sample.smooth_n_slide()
            # create list for each sample with results from any of the above methods
            sample.plot_growth_parameters()

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

    def create_dataframe(self):
        expt_df = pd.DataFrame(
            [   [sample.name,
                sample.lag_time,
                sample.growth_rate,
                sample.doubling_time,
                sample.time_of_max_rate,
                sample.sat_time,
                sample.max_OD,
                sample.time_of_max_OD]
            for sample in self.sample_list],
            columns = ("well", "lag time", "growth rate", "doubling time", "time of max growth rate",
                        "saturation time", "max OD", "time of max OD"))
        return expt_df

    def plot_histogram(self, expt_df):
        # plot histogram of growth rates for entire experiment
        plt.hist(expt_df["doubling time"])
        plt.xlabel("number of samples")
        plt.ylabel("doubling time")
        plt.savefig(self.name + '_histogram.png', dpi=300, bbox_inches='tight')

    def plot_heat_map(self, expt_df):
        # check if sample names are well IDs:
        if expt_df["well"][0] == "A01":
        expt_array = expt_df["growth rates"]
        # plot heat map of growth rates for entire experiment
        plt.pcolor(expt_array)

# create Sample class for attributes collected for each column of data
class Sample(object):

    def __init__(self, experiment, name, data):
        self.experiment = experiment # now ref attributes as self.experiment.attr
        self.elapsed_time = self.experiment.elapsed_time
        self.out_dir = self.experiment.out_dir
        self.name = name
        self.raw_data = data
        self.cal_data = self.raw_data - self.experiment.blank # subtract blank value, shouldn't be necessary...
        self.log_data = np.log(self.cal_data) # Log of the calibrated OD

        # get max OD and time of max_OD
        self.max_OD_index = np.argmax(self.raw_data)
        self.max_OD = self.raw_data[self.max_OD_index]
        self.time_of_max_OD = self.elapsed_time[self.max_OD_index]

        # get min OD for lag time calculation
        self.min_logOD = np.amin(self.log_data)

    def spline_max_growth_rate(self):

        interpolator = interpolate.UnivariateSpline(self.elapsed_time, self.log_data, k=4, s=0.05) #k can be 3-5
        der = interpolator.derivative()

        # Get the approximation of the derivative at all points
        der_approx = der(self.elapsed_time)

        # Get the maximum
        self.maximum_index = np.argmax(der_approx)
        self.growth_rate = der_approx[self.maximum_index]
        self.doubling_time = np.log(2)/self.growth_rate
        self.time_of_max_rate = self.elapsed_time[self.maximum_index]
        print self.name, self.growth_rate, self.time_of_max_rate

        # Get estimates of lag time and saturation time from 2nd derivative
        der2 = der.derivative()
        der2_approx = der2(self.elapsed_time)
        try: self.lag_index = signal.argrelmax(der2_approx)[0][0] # find first max
        except: self.lag_index = 0
        if self.lag_index > self.maximum_index: self.lag_index = 0
        self.lag_time = self.elapsed_time[self.lag_index]
        minima = signal.argrelmin(der2_approx)[0] # find first min after maximum_index
        which_min = bisect.bisect(minima,self.maximum_index)
        try: self.sat_index = minima[which_min]
        except: self.sat_index = -1
        self.sat_time = self.elapsed_time[self.sat_index]

        self.fit_y_values = interpolator(self.elapsed_time)
        # we have a slope (self.growth_rate) and a point on the line (self.time_of_max_rate, self.log_data[self.maximum_index])
        # we can draw a line representing the fit derivative:
        self.intercept = self.log_data[self.maximum_index] - (self.growth_rate*self.time_of_max_rate) # b = y - mx
        self.fit_on_orig = [(np.exp(self.growth_rate * x) * np.exp(self.intercept) + self.experiment.blank) \
                            for x in self.elapsed_time]

    def sliding_window(self, data_to_use):

        window_size, interval = self.experiment.get_window_size()

        rates = []
        intercepts = []
        r_values = []
        num_windows = len(data_to_use)-window_size+1
        for i in range(0, num_windows):
            window_times = self.elapsed_time[i:i+window_size]
            sub_data = data_to_use[i:i+window_size] # now use this sub_data to fit a line (not exp since we have the log2 data)
            results = stats.linregress(window_times, sub_data) # get fit parameters
            # save parameters
            rates.append(results[0]) #first two items in results are slope and intercept
            intercepts.append(results[1])
            r_values.append(results[2])

        maximum_rate = max(rates) # get max rate

        # extension: find other slopes at 95% of max rate, use all points to calculate new rate
        max95_rates = [rates.index(i) for i in rates if i >= 0.95*maximum_rate]
        start_point = max95_rates[0]
        end_point = max95_rates[-1] + window_size

        window_times = self.elapsed_time[start_point:end_point]
        sub_data = data_to_use[start_point:end_point]
        results = stats.linregress(window_times, sub_data)
        self.growth_rate = results[0]
        self.intercept = results[1]
        self.r_squared = results[2]**2

        self.maximum_index = int((end_point - start_point)/2) # returns midpoint of time window used to calc rate
        self.time_of_max_rate = self.elapsed_time[self.maximum_index]
        self.doubling_time = np.log(2)/self.growth_rate

        # get lag_time and saturation_time parameters - how? from slope-intercept of max slope
        print 'fit line is y = ' + str(self.growth_rate) + '*x + ' + str(self.intercept)
        print 'r-squared is ' + str(self.r_squared)
        self.lag_time = (np.amin(data_to_use) - self.intercept) / self.growth_rate # find x where y = min logOD
        self.lag_index = int(self.lag_time/interval)
        if self.lag_index > self.maximum_index: self.lag_index = 0
        self.sat_time = (np.amax(data_to_use) - self.intercept) / self.growth_rate # find x where y = max logOD
        self.sat_index = int(self.sat_time/interval)
        if self.sat_index > len(self.elapsed_time): self.sat_index = len(self.elapsed_time) - 1
        self.fit_y_values = [((self.growth_rate * x) + self.intercept) for x in self.elapsed_time]
        self.fit_on_orig = [(np.exp(self.growth_rate * x) * np.exp(self.intercept) + self.experiment.blank) \
                            for x in self.elapsed_time]

    def smooth_n_slide(self):

        # first get a spline:
        interpolator = interpolate.UnivariateSpline(self.elapsed_time, self.log_data, k=4, s=0.05) #k can be 3-5
        self.smooth_log_data = interpolator(self.elapsed_time)

        # Get the approximation of the derivative at all points
        der = interpolator.derivative()
        der_approx = der(self.elapsed_time)

        # now compute rate from sliding window using smoothed data points
        self.sliding_window(self.smooth_log_data)

    def plot_growth_parameters(self):

        fig = plt.figure()
        # orig data with points marked for max growth, lag time, saturation time, and max OD
        orig = fig.add_subplot(211)
        orig.plot(self.elapsed_time, self.raw_data, ls='', marker='.', label='raw data')
        orig.plot(self.time_of_max_rate, self.raw_data[self.maximum_index], 'ro', label='max growth')
        orig.plot(self.elapsed_time[self.lag_index], self.raw_data[self.lag_index], 'bo', label='lag time')
        orig.plot(self.elapsed_time[self.sat_index], self.raw_data[self.sat_index], 'go', label='saturation time')
        orig.plot(self.time_of_max_OD, self.raw_data[self.max_OD_index], 'ko', label='max OD')
        #if self.experiment.method == 'sliding_window':
        orig.plot(self.elapsed_time, self.fit_on_orig, 'r-', label='fit')

        #orig.set_xlabel('elapsed time (minutes)')
        orig.set_ylabel('OD600')
        orig.set_ylim(0,self.max_OD + 0.1)
        orig.legend(loc='best')

        # log(OD) data with fit line
        logOD = fig.add_subplot(212)
        logOD.plot(self.elapsed_time, self.log_data, ls='', marker='.', label='ln(OD)')
        logOD.autoscale(False) # don't want plot rescaled for fit line
        logOD.plot(self.elapsed_time, self.fit_y_values, 'r-', label='fit')
        logOD.set_xlabel('elapsed time (minutes)')
        logOD.set_ylabel('ln(OD600)')
        logOD.legend(loc='best')

        self.plot_file = self.name + '_plot.png'
        plot_path = os.path.join(self.out_dir, self.name)
        fig.savefig(plot_path, dpi=200, bbox_inches='tight')
        fig.clf()

    # put functions of all actual data analysis and plotting here - create list of Sample objects in Experiment


def main():
    parser = argparse.ArgumentParser(description = "Specify a data file to be analyzed.")
    parser.add_argument('-f', '--data_file', required=True, help='Full path to data file.')
    parser.add_argument('-p', '--plate_layout_file', help='Full path to file with plate layout information.')
    parser.add_argument('-b', '--blank', default=0, help='Either a blank value or well to calculate blank value.')
    parser.add_argument('-m', '--method', default='spline', help='Choose method used for data analysis.')
    parser.add_argument('-o', '--output_directory', help='Full path to output directory.')
    parser.add_argument('--organism', default='yeast', help='Organism in experiment [yeast, bacteria]; used in sliding_window method')
    args = parser.parse_args()

    data_file = args.data_file
    plate_layout = args.plate_layout_file # include expt date and run description here
    blank = args.blank
    out_dir = args.output_directory
    method = args.method
    organism = args.organism

    print "got the inputs"

    experiment = OD_growth_experiment(data_file, blank, method, organism, out_dir)
    print "created an experiment"

    experiment.get_sample_data()
    print "created samples"

    sample_output = experiment.analyze_sample_data()
    print "analyzed samples"

    # create pandas dataframe of sample info:
    expt_df = experiment.create_dataframe()

    # plots!
    experiment.plot_histogram(expt_df)
    experiment.plot_heat_map(expt_df)

    # write output to file!
    output_name = os.path.basename(data_file) + '_output.txt'
    new_file = os.path.join(out_dir, output_name)
    output_file = open(new_file,'w')

    output_file.write("well \t lag time \t max growth rate \t doubling time \t time of max rate \t saturation time " \
                      "\t max OD \t time of max OD \n") #\t plot_file
    for sample in sample_output:

        output_file.write("%s \t %i \t %f \t %f \t %i \t %i \t %f \t %i \n" % ( # \t %s
            sample.name, sample.lag_time, sample.growth_rate, doubling_time, sample.time_of_max_rate, sample.sat_time, \
            sample.max_OD, sample.time_of_max_OD #, sample.plot_file
        ))
    print "wrote the output"

if __name__ == "__main__":
    main()