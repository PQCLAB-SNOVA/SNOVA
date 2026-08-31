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

short_signatures = [
    ['SNOVA_I_X', 50, 17, 16, 2],
    ['SNOVA_III_X', 76, 25, 16, 2],
    ['SNOVA_V_X', 104, 33, 16, 2, 2, 33],
]

recommended += short_signatures

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

source_dir = '../src/'

for target in ['ref', 'opt', 'mem', 'avx2']:
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

        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            print(f'\tmake -C {name} kat', file=mf)

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

        for aes in aes_list:
            name = pname + '_AES' if aes else pname
            print(f'\t@make -C {name} speed', file=mf)

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
            os.makedirs(dirname)

            with open('Makefile.ref') as infile:
                data = infile.read()

            if target == 'ref':
                for file in ref_sources:
                    shutil.copyfile(source_dir + file, dirname + file)
            else:
                sources = opt_sources

                if target == 'mem':
                    snova_src = 'snova_memopt_16'
                    opt = 5
                    if o == 3 and l == 4 and r == 8:
                        sources += [
                            'abq_snova_i_k.h',
                        ]
                    elif o == 4 and l == 4 and r == 6:
                        sources += [
                            'abq_snova_i_b.h',
                        ]
                    elif o == 5 and l == 4 and r == 4:
                        sources += [
                            'abq_snova_i_s.h',
                        ]
                    elif o == 4 and l == 4 and r == 8:
                        sources += [
                            'abq_snova_iii_k.h',
                        ]
                    elif o == 5 and l == 4 and r == 6:
                        sources += [
                            'abq_snova_iii_b.h',
                        ]
                    elif o == 6 and l == 4 and r == 5:
                        sources += [
                            'abq_snova_iii_s.h',
                        ]
                    elif o == 4 and l == 5 and r == 8:
                        sources += [
                            'abq_snova_v_k.h',
                        ]
                    elif o == 5 and l == 5 and r == 6:
                        sources += [
                            'abq_snova_v_b.h',
                        ]
                    elif o == 6 and l == 5 and r == 5:
                        sources += [
                            'abq_snova_v_s.h',
                        ]
                    elif o == 17 and l == 2 and r == 2:
                        sources += [
                            'abq_snova_i_x.h',
                        ]
                    elif o == 25 and l == 2 and r == 2:
                        sources += [
                            'abq_snova_iii_x.h',
                        ]
                    elif o == 33 and l == 2 and r == 2:
                        sources += [
                            'abq_snova_v_x.h',
                        ]
                elif l == 2:
                    snova_src = 'snova_opt_16_2'
                    opt = 21
                elif target == 'avx2':
                    snova_src = 'snova_avx2_16'
                    opt = 20
                else:
                    snova_src = 'snova_opt_16'
                    opt = 10

                data = data.replace('symmetric_ref.o', 'symmetric.o')
                data = data.replace('snova_ref', snova_src)
                data = data.replace('-DSNOVA_OPT=0', f'-DSNOVA_OPT={opt}')

                for file in sources:
                    shutil.copyfile(source_dir + file, dirname + file)
                shutil.copyfile(source_dir + snova_src + '.c', dirname + snova_src + '.c')

            with open(dirname + 'Makefile', 'w') as outfile:
                outfile.write(data)

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
