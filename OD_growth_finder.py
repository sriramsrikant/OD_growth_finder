import pandas as pd
import numpy as np
from scipy import interpolate, signal, stats
import bisect
#from scipy.interpolate import UnivariateSpline
# import scipy as sp
import matplotlib.pyplot as plt
import datetime as datetime

def reformat_time(x):
    """Assumes x is a datetime.time object."""
    time_in_seconds = (60.*60.*x.hour + 60*x.minute + x.second)
    return time_in_seconds/(60.) # return time in minutes for computation of doubling time in minutes

class OD_growth_experiment(object):

    def __init__(self, path_to_data, blank, method, output_path = './', organism):
        self.path_to_data = path_to_data
        self.data = pd.read_excel(path_to_data) # not sure if read_csv might be better
        # Drop the rows that have NAN's, usually at the end
        self.data.dropna(inplace=True, axis=1)
        self.sample_list = []
        self.method = method
        self.organism = organism

        # Get the times from the data - check for format and reformat if necessary
        self.times = self.data.iloc[:,0] # just assume first column is time data
        if self.times.dtype != 'int64': # if numeric, assumes minutes
            self.times = self.times.astype(datetime.time)
            self.elapsed_time = self.times.apply(reformat_time)
        else: self.elapsed_time = self.times

        self.data.drop(self.data.columns[[0]], inplace=True, axis=1) # remove time column from data to be analyzed

        # check blank
        if type(blank) == 'str':
            self.blank = self.data.mean()[blank] # if blank input is well, blank is average value of that well
        else: self.blank = blank

        # Set the output path
        self.output_path = output_path

        # Set the default s for fitting...deals with how close the fit is to the points
        #self.s = s

    def get_sample_data(self):
        for well_str in self.data.columns: # well_str is name of column
            # create Sample object Sample(experiment, data) and add to sample list
            raw_data = self.data.loc[:, well_str]
            self.sample_list.append(Sample(self, well_str, raw_data))
        return self.sample_list

    def analyze_sample_data(self):

        for sample in self.sample_list:
            if self.method == 'spline':
                sample.spline_max_growth_rate()
#            elif self.method == "lowess":
            elif self.method == "sliding_window":
                sample.sliding_window()
            # create list for each sample with results from any of the above methods
            sample.plot_growth_parameters()

        return self.sample_list


# create Sample class for attributes collected for each column of data
class Sample(object):

    def __init__(self, experiment, name, data):
        self.experiment = experiment # now ref attributes as self.experiment.attr
        self.elapsed_time = self.experiment.elapsed_time
        self.name = name
        self.raw_data = data
        # self.cal_data = self.raw_data - self.experiment.blank # subtract blank value, shouldn't be necessary...
        self.log_data = np.log(self.raw_data) # Log of the calibrated OD

        # get max OD and time of max_OD
        max_OD_index = np.argmax(self.raw_data)
        self.max_OD = self.raw_data[max_OD_index]
        self.time_of_max_OD = self.elapsed_time[max_OD_index]

    def spline_max_growth_rate(self):

        interpolator = interpolate.UnivariateSpline(self.elapsed_time, self.log2_data, k=4, s=0.05) #k can be 3-5
        der = interpolator.derivative()

        # Get the approximation of the derivative at all points
        der_approx = der(self.elapsed_time)

        # Get the maximum
        maximum_index = np.argmax(der_approx)
        self.maximum_rate = der_approx[maximum_index]
        self.time_of_max_rate = self.elapsed_time[maximum_index]

        # Get estimates of lag time and saturation time from 2nd derivative
        der2 = der.derivative()
        der2_approx = der2(self.elapsed_time)
        lag_index = signal.argrelmax(der2_approx)[0][0] # find first max
        self.lag_time = self.elapsed_time[lag_index]
        minima = signal.argrelmin(der2_approx)[0]
        sat_index = minima[bisect.bisect(minima,lag_index)] # find first min after lag_index
        self.sat_time = self.elapsed_time[sat_index]

        self.fit_y_values = interpolator(self.elapsed_time)

    def sliding_window(self):
        #determine a good window size - start with 1.5*(wt doubling time in minutes)/time interval in minutes
        # determine time interval from self.elapsed_time
        interval = self.elapsed_time[1] - self.elapsed_time[0]
        if self.experiment.organism == 'yeast':
            window_size = int(1.5*90/interval)
        elif self.experiment.organism == 'bacteria':
            window_size = int(1.5*20/interval)
        else: window_size = 5 # print "I don't know what organism you are using. Here's some data anyway."
        window_times = self.elapsed_time[:window_size] # always use first time points

        rates = []
        intercepts = []
        num_windows = len(self.log_data)-window_size+1
        for i in range(0, num_windows):
            sub_data = self.log_data[i:i+window_size] # now use this sub_data to fit a line (not exp since we have the log2 data)
            results = stats.linregress(window_times, sub_data) # get fit parameters
            # save b parameter
            rates.append(results[0]) #first two items in results are slope and intercept
            intercepts.append(results[1])

        self.maximum_rate = max(rates) # get max rate
        maximum_index = rates.index(self.maximum_rate) + int(window_size/2) # returns midpoint of time window for max rate
        self.time_of_max_rate = self.elapsed_time[maximum_index]

        # get lag_time and saturation_time parameters - how? from slope-intercept of max slope
        intercept_max = intercepts[maximum_index]
        self.lag_time = -(intercept_max) / self.maximum_rate # find x where y = 0
        self.sat_time = (self.max_OD - intercept_max) / self.maximum_rate # find x where y = max_OD

        self.fit_y_values = [(self.maximum_rate*x + intercept_max) for x in self.elapsed_time]

    def plot_growth_parameters(self):

        fig = plt.figure()
        # orig data with points marked for max growth, lag time, saturation time, and max OD
        orig = fig.add_subplot(211)
        orig.plot(self.elapsed_time, self.raw_data, ls='', marker='.', label='raw data')
        orig.plot(self.time_of_max_rate, self.raw_data[self.time_of_max_rate], 'ro', label='max growth', alpha=0.5)
        orig.plot(self.lag_time, self.raw_data[self.lag_time], 'bo', label='lag time')
        orig.plot(self.sat_time, self.raw_data[self.sat_time], 'go', label='saturation time')
        orig.plot(self.time_of_max_OD, self.raw_data[self.time_of_max_OD], 'ko', label='max OD')

        orig.xlabel('elapsed time (minutes)')
        orig.ylabel('OD600')
        orig.legend(loc='best')

        # log(OD) data with fit line
        logOD = fig.add_subplot(212)
        logOD.plot(self.elapsed_time, self.log_data, ls='', marker='.', label='ln(OD)')
        logOD.plot(self.elapsed_time, self.fit_y_values, 'r-', label='fit')
        logOD.xlabel('elapsed time (minutes)')
        logOD.ylabel('ln(OD600)')
        logOD.legend(loc='best')

        fig.savefig(self.name + '_plot.png', dpi=200, bbox_inches='tight')
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

    experiment = OD_growth_experiment(data_file, blank, method, out_dir, organism)
    experiment.get_sample_data()
    sample_output = experiment.analyze_sample_data()

    return pd.DataFrame(growth_rate_data, columns=['well', 'lag_time', 'max_growth_rate', 'time_of_max_rate', \
                                                       'saturation_time', 'max_OD', 'time_of_max_OD', 'plot_file'])
