from subprocess import call
import shutil
import os

params = [
    [24, 5, 23, 4, False],
    [37, 8, 19, 4, False],
    [60, 10, 23, 4, False],

    [24, 5, 16, 4, False],
    [37, 8, 16, 4, False],
    [60, 10, 16, 4, False],

    [24, 5, 23, 4, True],
    [37, 8, 19, 4, True],
    [60, 10, 23, 4, True],

    [24, 5, 16, 4, True],
    [37, 8, 16, 4, True],
    [60, 10, 16, 4, True],
]

ref_sources = [
    '.gitignore',
    'LICENSE',
    'PQCgenKAT_sign.c',
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

source_dir = '../oddqsrc/'

for target in ['ref', 'opt', 'avx2']:
    shutil.rmtree(target, ignore_errors=True)
    os.makedirs(target)
    shutil.copy('Makefile.root', target + '/Makefile')
    shutil.copy('../KAT/digest.sh', target)

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
            os.remove(dirname + 'snova_avx2_q.c')
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
