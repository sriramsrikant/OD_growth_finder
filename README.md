# OD_growth_finder

Global functions:

reformat_time(x)  called by OD_growth_experiment.__init__ for time input

analyze_experiment(data_file, plate_layout=None, blank=0, method='smooth_n_slide', organism='yeast', sample_plots=False,
                   out_dir='./', s=0.05)
        
make_plots(experiment, show=True, save=False, metric1='doubling time', unit='minutes', metric2='growth rate')

compute_means(experiment, metric=['growth rate'], save=False, keys=None)

main()
    flags:
    '-f', '--data_file', required=True, help='Full path to data file.'
    '-p', '--plate_layout_file', default=None, help='Full path to file with plate layout information.'
    '-b', '--blank', default=0, help='Either a blank value or well to calculate blank value.'
    '-m', '--method', default='smooth_n_slide', help='Choose method used for data analysis.'
    '-o', '--output_directory', default='./', help='Full path to output directory.'
    '--organism', default='yeast', help='Organism in experiment [yeast, bacteria]; used in sliding window methods'


class functions

class OD_growth_experiment(path_to_data, plate_layout=None, blank=None, organism='yeast', out_dir='./'):
    create_sample_list(self)  creates Sample objects
    analyze_sample_data(self, method='smooth_n_slide', sample_plots=False, s=0.05)
    get_window_size(self) # called by other methods
    output_data(self, save=False)
    plot_histogram(self, show=True, save=False, metric='doubling time', unit='minutes')
    plot_heatmap(self, show=True, save=False, metric='growth rate', unit=None)
    
class Sample(experiment, name, data):
    spline_max_growth_rate(self, s)
    sliding_window(self, data=None)  -> if data is None: data = self.log_data
    get_max_growth_parameters(self, data=None)  -> if data is None: data = self.log_data, calls sliding_window
    get_lag_sat_parameters(self)  calls sliding_window, get_flex_point
    get_flex_point(self, rates, point, window_size)
    smooth_n_slide(self, s)  calls get_max_growth_parameters, get_lag_sat_parameters
    plot_growth_parameters(self, show=True, save=False, folder='./')