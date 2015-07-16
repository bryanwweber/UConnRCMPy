# System imports
from pathlib import Path

# Third-party imports
import matplotlib.pyplot as plt

# Local imports
from .ignitiondelayexp import ExperimentalIgnitionDelay
from .utilities import copy


class LoopIgnitionDelay(ExperimentalIgnitionDelay):

    def __init__(self, filename):
        super(ExperimentalIgnitionDelay, self).__init__(filename)
        self.calculate_ignition_delay()
        self.calculate_EOC_temperature()


def main(path='.', plot=False):
    p = Path(path)
    result = []

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

    for f in p.glob('[0-3]*.txt'):
        print(f)
        case = LoopIgnitionDelay(f.resolve())
        result.append('\t'.join(map(str, [
            case.time_of_day, case.pin, case.Tin, case.p_EOC,
            case.ignition_delay, case.first_stage, case.T_EOC,
            case.spacers, case.shims])))
        if plot:
            ax.plot(case.ztim, case.pressure, label=case.date)

    copy('\n'.join(sorted(result)))
    print('Finished')

if __name__ == '__main__':
    main()
