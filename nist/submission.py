# SPDX-License-Identifier: MIT
#
# Script to create the recommended SNOVA instances
#
# Copyright (c) 2026 SNOVA TEAM

import math
import os
import shutil


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

alternatives = [
    # SL 1
    ['snova1a', 27, 3, 16, 4, 7],
    ['snova1b', 27, 4, 16, 4, 5],
    ['snova1c', 27, 4, 19, 4, 6, 5],
    ['snova1d', 27, 5, 13, 4, 4, 5],
    ['snova1x', 50, 17, 16, 2],

    # SL 3
    ['snova3a', 27, 4, 16, 5, 5],
    ['snova3b', 33, 5, 16, 4, 5],
    ['snova3c', 39, 5, 11, 4, 6],
    ['snova3d', 38, 5, 19, 4, 6, 7],
    ['snova3x', 76, 25, 16, 2],

    # SL 5
    ['snova5a', 40, 5, 19, 5, 6, 6],
    ['snova5b', 47, 6, 19, 4, 6, 9],
    ['snova5x', 104, 33, 16, 2, 2, 33],
]

# recommended += alternatives

aes_list = [False, True]

gen_sources = [
    'LICENSE',
    'PQCgenKAT_sign.c',
    'aes.c',
    'api.h',
    'keccak_opt64.h',
    'rng.c',
    'rng.h',
    'sign.c',
    'snova.h',
    'speed.c',
    'symmetric.h',
]

ref_sources = gen_sources + [
    'snova_ref.c',
    'symmetric_ref.c',
]

opt_sources = gen_sources + [
    'keccak_avx2.h',
    'keccak_avx512.h',
    'symmetric.c',
]

source_dir = 'generic/'

for target in ['ref', 'opt', 'avx2', 'gfni']:
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)

    shutil.copyfile('digests.sh', target + '/digests.sh')

    mf = open(target + '/Makefile', 'w')
    print('MAKEFLAGS += --no-print-directory\n\nall: kat\n\ndigest:\n\t@sh digests.sh\n', file=mf)

    print('\nkat:', file=mf)
    for paramset in recommended:
        pname = paramset[0]
        if paramset[1] == 'Sym':
            param = paramset[2:]
            symmetric = True
        else:
            param = paramset[1:]
            symmetric = False

        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        if len(param) > 4:
            r = param[4]
        else:
            r = l

        avxsource = target in ['avx2', 'gfni'] and r != 2
        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            print(f'\tmake -C {name} nistkat' if avxsource else f'\tmake -C {name} kat', file=mf)

    print('\nspeed:', file=mf)
    for paramset in recommended:
        pname = paramset[0]
        if paramset[1] == 'Sym':
            param = paramset[2:]
            symmetric = True
        else:
            param = paramset[1:]
            symmetric = False

        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        if len(param) > 4:
            r = param[4]
        else:
            r = l

        avxsource = target in ['avx2', 'gfni'] and r != 2
        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            print(f'\t@make -C {name} {"bench" if avxsource else "speed"}', file=mf)

    print('\nclean:', file=mf)
    for paramset in recommended:
        pname = paramset[0]
        if paramset[1] == 'Sym':
            param = paramset[2:]
            symmetric = True
        else:
            param = paramset[1:]
            symmetric = False

        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        if len(param) > 4:
            r = param[4]
        else:
            r = l
        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            print(f'\tmake -C {name} clean', file=mf)

    mf.close()

    for paramset in recommended:
        pname = paramset[0]
        if paramset[1] == 'Sym':
            param = paramset[2:]
            symmetric = True
        else:
            param = paramset[1:]
            symmetric = False

        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        if len(param) > 4:
            r = param[4]
            if len(param) > 5:
                m1 = param[5]
            else:
                m1 = math.floor((o * r) / l)
        else:
            r = l
            m1 = math.floor((o * r) / l)

        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            dirname = target + '/' + name + '/'

            if target in ['avx2', 'gfni'] and r != 2:
                coredir = 'core/'
                shutil.copytree(coredir, dirname)

                with open(dirname + 'snova_config.mk', 'w') as outfile:
                    print('ARCH=x86_avx2', file=outfile)
                    print(f'GFNI={1 if target == 'gfni' else 0}', file=outfile)
                    print(f'SNOVA_V={v}', file=outfile)
                    print(f'SNOVA_O={o}', file=outfile)
                    print(f'SNOVA_Q={q}', file=outfile)
                    print(f'SNOVA_L={l}', file=outfile)
                    print(f'SNOVA_R={r}', file=outfile)
                    print(f'SNOVA_M1={m1}', file=outfile)
                    print(f'AES={1 if aes else 0}', file=outfile)
                    print(f'SNOVA_NAME={name}', file=outfile)

                continue

            if target == 'ref':
                os.makedirs(dirname)
                shutil.copyfile('Makefile.ref', dirname + 'Makefile')
                for file in ref_sources:
                    shutil.copyfile(source_dir + file, dirname + file)
            else:
                if l == 2:
                    snova_src = 'snova_opt_16_2' if q == 16 else 'snova_opt_q_r'
                elif l == 5:
                    snova_src = 'snova_opt_16_5' if q == 16 else 'snova_opt_q_r'
                elif r != l:
                    snova_src = 'snova_opt_16' if q == 16 else 'snova_opt_q_r'
                else:
                    snova_src = 'snova_opt_16' if q == 16 else 'snova_opt_q_s'

                os.makedirs(dirname)
                with open('Makefile.opt') as infile:
                    data = infile.read()
                with open(dirname + 'Makefile', 'w') as outfile:
                    outfile.write(data.replace('snova_opt_q', snova_src))
                for file in opt_sources:
                    shutil.copyfile(source_dir + file, dirname + file)
                shutil.copyfile(source_dir + snova_src + '.c', dirname + snova_src + '.c')

            sp = open(dirname + 'snova_params.h', 'w')
            print('#define SNOVA_NAME', name, file=sp)
            print('#define SNOVA_v', v, file=sp)
            print('#define SNOVA_o', o, file=sp)
            print('#define SNOVA_q', q, file=sp)
            print('#define SNOVA_l', l, file=sp)
            print('#define SNOVA_r', r, file=sp)
            print('#define SNOVA_m1', m1, file=sp)
            if len(param) > 6:
                print('#define SNOVA_alpha', param[6], file=sp)
            if aes:
                print('#define AESCTR', file=sp)
            if symmetric:
                print('#define SYMMETRIC', file=sp)
            sp.close()
