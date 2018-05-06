
from itertools import product
import numpy as np
from string import ascii_uppercase


_clades = {}

alphabet = list(ascii_uppercase)

for i in range(26*3+26):
    used_ids = list(_clades.keys())
    potential_clade_name = np.setdiff1d(alphabet, used_ids)

    size = 1

    while len(potential_clade_name) == 0:
        a = product(ascii_uppercase, repeat=size)
        a = [''.join(s) for s in a]
        potential_clade_name = np.setdiff1d(a, used_ids)
        size += 1

    _clades[potential_clade_name[0]] = -1

print(_clades)
