from .ignitiondelayexp import ExperimentalIgnitionDelay
from .pressure_traces import ReactivePressureTrace
from .utilities import copy
from pathlib import Path


class LoopIgnitionDelay(ExperimentalIgnitionDelay):

    def __init__(self, filename):
        super(ExperimentalIgnitionDelay, self).__init__(filename)
        self.calculate_ignition_delay()
        self.calculate_EOC_temperature()


def main(path='.'):
    p = Path(path)
    result = []

    for f in p.glob('[0-3]*.txt'):
        print(f)
        case = LoopIgnitionDelay(f.resolve())
        result.append('\t'.join(map(str, [
            case.time_of_day, case.pin, case.Tin, case.p_EOC,
            case.ignition_delay, case.first_stage, case.T_EOC,
            case.spacers, case.shims])))

    copy('\n'.join(sorted(result)))
    print('Finished')

if __name__ == '__main__':
    main()
