from subprocess import call
import shutil
import os

params = [
    [43, 17, 16, 2, True],
    [43, 17, 16, 2, False],
    [24, 5, 16, 4, True],
    [24, 5, 16, 4, False],
    [24, 5, 23, 4, True],
    [24, 5, 23, 4, False],

    [69, 25, 16, 2, True],
    [69, 25, 16, 2, False],
    [37, 8, 16, 4, True],
    [37, 8, 16, 4, False],
    [24, 5, 16, 5, True],
    [24, 5, 16, 5, False],
    [37, 8, 19, 4, True],
    [37, 8, 19, 4, False],

    [99, 33, 16, 2, True],
    [99, 33, 16, 2, False],
    [60, 10, 16, 4, True],
    [60, 10, 16, 4, False],
    [29, 6, 16, 5, True],
    [29, 6, 16, 5, False],
    [60, 10, 23, 4, True],
    [60, 10, 23, 4, False],
]


mf = open('Makefile.root', 'w')
print('''MAKEFLAGS += --no-print-directory

all: kat

digest:
	@sh digest.sh */PQCsign*''', file=mf)

print('\nkat:', file=mf)
for param in params:
    v, o, q, l, aes = param
    print('''\tmake -C SNOVA_{}_{}_{}_{}{} kat'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nspeed:', file=mf)
for param in params:
    v, o, q, l, aes = param
    print('''\t@make -C SNOVA_{}_{}_{}_{}{} speed'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

print('\nclean:', file=mf)
for param in params:
    v, o, q, l, aes = param
    print('''\tmake -C SNOVA_{}_{}_{}_{}{} clean'''.format(
        v, o, q, l, '_AES' if aes else ''), file=mf)

mf.close()

ref_sources = [
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
    'snova_ref.c',
    'speed.c',
    'symmetric.h',
    'symmetric_ref.c',
]

source_dir = '../src/'

for target in ['ref', 'opt', 'avx2']:
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)
    shutil.copy('Makefile.root', target + '/Makefile')
    shutil.copy('digest.sh', target)

    for param in params:
        v, o, q, l, aes = param

        dirname = target + '/SNOVA_{}_{}_{}_{}{}/'.format(
            v, o, q, l, '_AES' if aes else '')

        if target == 'ref':
            os.makedirs(dirname)
            shutil.copy('README.ref', dirname + 'README.md')
            shutil.copy('Makefile.ref', dirname + 'Makefile')
            for file in ref_sources:
                shutil.copy(source_dir + file, dirname)
        else:
            shutil.copytree(source_dir, dirname)

        if target == 'opt':
            os.remove(dirname + 'snova_avx2.c')
            os.remove(dirname + 'snova_avx2_q.c')
            os.remove(dirname + 'snova_avx2_16.c')
            shutil.copy('README.opt', dirname + 'README.md')
            shutil.copy('Makefile.opt', dirname + 'Makefile')

        sp = open(dirname + 'snova_params.h', 'w')
        print('#define SNOVA_v', v, file=sp)
        print('#define SNOVA_o', o, file=sp)
        print('#define SNOVA_q', q, file=sp)
        print('#define SNOVA_l', l, file=sp)
        if aes:
            print('#define AESCTR', file=sp)
        sp.close()

os.remove('Makefile.root')
