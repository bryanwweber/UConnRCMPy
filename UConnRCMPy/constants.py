import cantera as ct

one_atm_in_torr = 760.0
one_atm_in_bar = 1.01325
one_bar_in_pa = 1E5

try:
    cantera_version = map(int, ct.__version__.split('.'))
except ValueError:
    cantera_version = map(int, ct.__version__.split('.')[0:2])
    cantera_version.append(int(ct.__version__.split('.')[2][0]))
