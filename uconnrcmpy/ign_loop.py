# System imports
from pathlib import Path

# Third-party imports
import matplotlib.pyplot as plt

# Local imports
from .dataprocessing import Experiment
from .utilities import copy


def main(path='.', plot=False):
    p = Path(path)
    result = []

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('[0-3]*.txt'):
        print(f)
        case = Experiment(f.resolve())
        result.append('\t'.join(map(str, [
            case.experiment_parameters['time_of_day'], case.experiment_parameters['pin'],
            case.experiment_parameters['Tin'], case.pressure_trace.p_EOC, case.ignition_delay,
            case.first_stage, case.T_EOC, case.experiment_parameters['spacers'],
            case.experiment_parameters['shims']])))
        if plot:
            ax.plot(case.zeroed_time, case.pressure, label=case.experiment_parameters['date'])

    copy('\n'.join(sorted(result)))
    print('Finished')

if __name__ == '__main__':
    main()
