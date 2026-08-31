# SPDX-License-Identifier: MIT
#
# Generate response KAT digest from SageMath SNOVA
#
# Copyright (c) 2026 SNOVA TEAM

import math
import os
import subprocess
import sys

from hashlib import shake_256


recommended = [
    ['SNOVA_I_K', 29, 3, 16, 4, 8, 5],
    ['SNOVA_I_B', 27, 4, 16, 4, 6, 5],
    ['SNOVA_I_S', 27, 5, 16, 4, 4, 5],

    ['SNOVA_III_K', 38, 4, 16, 4, 8, 7],
    ['SNOVA_III_B', 38, 5, 16, 4, 6, 7],
    ['SNOVA_III_S', 38, 6, 16, 4, 5, 7],

    ['SNOVA_V_K', 40, 4, 16, 5, 8, 6],
    ['SNOVA_V_B', 40, 5, 16, 5, 6, 6],
    ['SNOVA_V_S', 40, 6, 16, 5, 5, 6],
]


for paramset in recommended:
    pname = paramset[0]

    if paramset[1] == 'Sym':
        var = paramset[2:]
        symmetric = True
    else:
        var = paramset[1:]
        symmetric = False

    v = var[0]
    o = var[1]
    q = var[2]
    l = var[3]

    if len(var) > 4:
        r = var[4]
    else:
        r = l

    if len(var) > 5:
        m1 = var[5]
    else:
        m1 = math.floor(o * r / l)

    for aes in [True, False]:
        name = pname + '_AES' if aes else pname

        with open('snova.sage') as infile:
            data = infile.read()

        # Create sage file from parameters
        param_dict = {
            "print('# SNOVA', v, o, q, l, r, m1, m2, 'AES' if aes else 'SHAKE', n_alpha)":
            f"print('# {name}')",
            "v = 27": f"v = {v}",
            "o = 5": f"o = {o}",
            "q = 16": f"q = {q}",
            "l = 4": f"l = {l}",
            "r = l": f"r = {r}",
            "m1 = 5": f"m1 = {m1}",
            "aes = False": f"aes = {aes}",
            "range(1)": "range(100)",
        }
        for key in param_dict.keys():
            data = data.replace(key, param_dict[key])
        if symmetric:
            data = data.replace('ASYMMETRIC_PUBMAT = True', 'ASYMMETRIC_PUBMAT = False')

        with open(f'_{name}.sage', 'w') as outfile:
            outfile.write(data)

        # Run and print digest

        result = subprocess.run(['sage', f'_{name}.sage'], capture_output=True, text=True)
        digest = shake_256(result.stdout.encode()).digest(24).hex()
        print(f'{name}.rsp  {digest}')
        sys.stdout.flush()

        os.remove(f'_{name}.sage')
