from .ignitiondelayexp import ExperimentalIgnitionDelay
from .pressure_traces import ReactivePressureTrace
from .utilities import copy
import os


class LoopIgnitionDelay(ExperimentalIgnitionDelay):

    def __init__(self, filename):
        super(ReactivePressureTrace, self).__init__(filename)
        self.process_pressure_trace(filename)
        self.calculate_ignition_delay()
        self.calculate_EOC_temperature()


def main():
    pth = os.listdir('.')
    result = []

    for filename in pth:
        if filename[:1].isdigit():
            print(filename)
            case = LoopIgnitionDelay(filename)
            result.append('\t'.join(map(str, [
                case.time_of_day, case.pin, case.Tin, case.p_EOC,
                case.ignition_delay, case.first_stage, case.T_EOC,
                case.spacers, case.shims])))

    copy('\n'.join(['\t'.join(map(str, res)) for res in sorted(result)]))
    print('Finished')

if __name__ == '__main__':
    main()
