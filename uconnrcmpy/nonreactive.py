"""
Nonreactive pressure trace processor
"""
# System imports
from pathlib import Path

# Third-party imports
import numpy as np
import matplotlib.pyplot as plt
import yaml

# Local imports
from .pressure_traces import NonReactivePressureTrace, ReactivePressureTrace
from .utilities import copy


class NonReactiveExperiments(object):
    """
    Class storing the non-reactive experiments corresponding to a
    particular reactive experiment.
    """

    def __init__(self):
        """
        When the class is initialized, load the reactive pressure trace
        and then add one non-reactive pressure trace.
        """
        with open('volume-trace.yaml') as yaml_file:
            self.yaml_data = yaml.load(yaml_file)
            self.reacfile = Path(self.yaml_data['reacfile'])

        self.reactive_case = ReactivePressureTrace(self.reacfile)
        """Associated reactive experiment."""

        self.fig = plt.figure('Non-Reactive Comparison')
        """The figure for plotting."""
        self.ax = self.fig.add_subplot(1, 1, 1)
        """The axes for plotting."""
        self.ax.plot(self.reactive_case.ztim, self.reactive_case.pressure)

        self.nonreactive_cases = {}
        """Dictionary of the cases of non-reactive experiments."""
        self.add_nonreactive_experiment()
        m = plt.get_current_fig_manager()
        m.window.showMaximized()

    def add_nonreactive_experiment(self):
        """
        Add a non-reactive experiment to this case. First process the
        pressure trace, then add the time/pressure array to the
        dictionary of all cases.
        """
        nr_trace = NonReactivePressureTrace()
        case_name = '{}'.format(nr_trace.date)
        self.nonreactive_cases[case_name] = np.vstack(
            (nr_trace.ztim, nr_trace.pressure)).transpose()

        copy('\t'.join(map(str,
                           [nr_trace.time_of_day, nr_trace.pin, nr_trace.Tin,
                            nr_trace.p_EOC, 'NR', 'NR', '', nr_trace.spacers,
                            nr_trace.shims])
                       )
             )

        self.ax.plot(nr_trace.ztim, nr_trace.pressure)

    def __call__(self):
        """
        When the instance of the class is called, add a non-reactive
        experiment.
        """
        self.add_nonreactive_experiment()
