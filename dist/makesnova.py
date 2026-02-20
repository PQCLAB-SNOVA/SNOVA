import shutil
import os

params = [
    [26, 5, 19, 4],
    [24, 5, 16, 4],
    [48, 17, 16, 2],

    [37, 8, 19, 4],
    [37, 8, 16, 4],
    [72, 25, 16, 2],

    [60, 10, 19, 4],
    [60, 10, 16, 4],
    [97, 33, 16, 2],
]


mf = open('Makefile.root', 'w')
print('''MAKEFLAGS += --no-print-directory

all: kat

digest:
	@sh digest.sh */PQCsign*''', file=mf)

print('\nkat:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    for aes in [True, False]:
        print('''\tmake -C SNOVA_{}_{}_{}_{}{} kat'''.format(
            v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nspeed:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    for aes in [True, False]:
        print('''\t@make -C SNOVA_{}_{}_{}_{}{} speed'''.format(
            v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nclean:', file=mf)
for param in params:
    v = param[0]
    o = param[1]
    q = param[2]
    l = param[3]
    for aes in [True, False]:
        print('''\tmake -C SNOVA_{}_{}_{}_{}{} clean'''.format(
            v, o, q, l, '_AES' if aes else ''), file=mf)

mf.close()

gen_sources = [
    '.gitignore',
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

for target in ['ref', 'opt', 'avx2']:
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)
    shutil.copy('Makefile.root', target + '/Makefile')
    shutil.copy('digest.sh', target)

    for param in params:
        v = param[0]
        o = param[1]
        q = param[2]
        l = param[3]
        if len(param) > 4:
            r = param[4]
        else:
            r = l
        for aes in [True, False]:

            dirname = target + '/SNOVA_{}_{}_{}_{}{}/'.format(
                v, o, q, l, '_AES' if aes else '')

            if target == 'ref':
                os.makedirs(dirname)
                shutil.copy('Makefile.ref', dirname + 'Makefile')
                for file in ref_sources:
                    shutil.copy(source_dir + file, dirname)
            else:
                if r != l:
                    snova_src = 'snova_opt_16' if q == 16 else 'snova_rect_q'
                elif q == 16:
                    snova_src = 'snova_opt_16' if target == 'opt' or l == 4 else 'snova_avx2_16'
                else:
                    snova_src = 'snova_opt_q'

                os.makedirs(dirname)
                with open('Makefile.opt') as infile:
                    data = infile.read()
                with open(dirname + 'Makefile', 'w') as outfile:
                    outfile.write(data.replace('snova_opt_q', snova_src))
                for file in opt_sources:
                    shutil.copy(source_dir + file, dirname)
                shutil.copy(source_dir + snova_src + '.c', dirname)

            sp = open(dirname + 'snova_params.h', 'w')
            print('#define SNOVA_v', v, file=sp)
            print('#define SNOVA_o', o, file=sp)
            print('#define SNOVA_q', q, file=sp)
            print('#define SNOVA_l', l, file=sp)
            if aes:
                print('#define AESCTR', file=sp)
            if len(param) > 4:
                print('#define SNOVA_r', param[4], file=sp)
            if len(param) > 5:
                print('#define SNOVA_m1', param[5], file=sp)
            if len(param) > 6:
                print('#define SNOVA_alpha', param[6], file=sp)
            sp.close()

os.remove('Makefile.root')
