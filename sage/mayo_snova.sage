# SPDX-License-Identifier: MIT
#
# SNOVA Round 2 implementation in SageMath
#
# This version uses the formulation in terms of an E matrix
# which is derived from the ABQ matrices or from a seed
#
# This version also generates Rund 2 MAYO KATs
#
# Copyright (c) 2025 SNOVA TEAM

import hashlib
import traceback

try:
    import nistrng
except:
    print('Error importing nistrng')
    print('Try: export PYTHONPATH=`pwd`')
    quit()

################################################################

# SNOVA parameters

if True:
    # Round 2 SNOVA_24_5_16_4

    v = 24
    o = 5
    q = 16
    l = 4
    aes = False
    r = l
    m1 = o
    m2 = o * l * r
    RANDOM_E = False

else:
    # Round 2 MAYO-1

    v = 78
    o = 8
    q = 16
    l = 1
    aes = True
    r = 10
    m1 = 78
    m2 = o * l * r
    RANDOM_E = False

################################################################

n_alpha = r * r + r
n = v + o

MAYO = l == 1


# Set GF

GF_q = GF(q, 'x')
x = GF_q.gen()

if q == 16:
    def from_int(x): return GF_q.from_integer(x)
    def to_int(x): return x.to_integer()
else:
    def from_int(x): return x
    def to_int(x): return int(x)


# Set constants

if q == 16:
    PACK_GF = 2
    PACK_BYTES = 1

elif q == 19:
    Q_A = 1
    Q_B = 3
    Q_C = 15
    PACK_GF = 15
    PACK_BYTES = 8

elif q == 23:
    Q_A = 1
    Q_B = 11
    Q_C = 22
    PACK_GF = 7
    PACK_BYTES = 4


# Derived constants

def BYTES_GF(x):
    return (PACK_BYTES * (x) + PACK_GF - 1) // PACK_GF


GF16_HASH = m2
BYTES_HASH = BYTES_GF(GF16_HASH)

NUM_GEN_PUB_GF = m1 * (v * v + 2 * v * o) * l**2 + o * n_alpha * (r * (r + l) + 2 * l)
NUMGF_PK = m1 * o * l * (o * l)

if q == 16:
    NUM_GEN_PUB_BYTES = NUM_GEN_PUB_GF // 2
else:
    NUM_GEN_PUB_BYTES = NUM_GEN_PUB_GF


# Create the S matrix

if q == 16:
    S = matrix(GF_q, l, l, lambda i, j: from_int(abs(8 - (i + j))))
    if l == 5:
        S[4, 4] = from_int(9)
else:
    S = matrix(GF_q, l, l)
    for i in range(l):
        for j in range(i, l):
            S[i, j] = (Q_A + i + j) & Q_B
            S[j, i] = S[i, j]
    S[l - 1, l - 1] = Q_C


S_times_o = matrix.block_diagonal([S for _ in range(o)])
S_times_v = matrix.block_diagonal([S for _ in range(v)])
S_times_n = matrix.block_diagonal([S for _ in range(n)])


# Utils

def expand_gf(data, num):
    # Convert bytes to elements of $\mathbb{F}_{q}$
    res = []
    idx = 0
    while idx < len(data):
        sum = 0
        for i in range(min(PACK_BYTES, len(data) - idx)):
            sum += int(data[idx + i]) * int(256**i)
        idx += PACK_BYTES
        if q == 16:
            res.append(from_int(int(sum) % 16))
            res.append(from_int(int(sum / 16) % 16))
            sum = sum // 256
        else:
            for i in range(PACK_GF):
                res.append(from_int(int(sum) % q))
                sum = sum // q
    return res[:num]


def compress_gf(data, num):
    # Convert elements of $\mathbb{F}_{q}$ to bytes
    res = []
    idx = 0
    while idx < len(data):
        sum = 0
        for i in range(min(PACK_GF, len(data) - idx)):
            sum += to_int(data[idx + i]) * int(q**i)
        idx += PACK_GF
        for i in range(PACK_BYTES):
            res.append(int(sum) % 256)
            sum = sum // 256
    return bytes(res[:BYTES_GF(num)])


def hash_combined(msg, pk_seed, salt):
    # Get message hash in $\mathbb{F}_{q}$
    state = hashlib.shake_256()
    if MAYO:
        dgst = hashlib.shake_256(msg).digest(32)
        state.update(dgst)
    else:
        state.update(pk_seed)
        dgst = hashlib.shake_256(msg).digest(64)
        state.update(dgst)
    state.update(salt)
    res = state.digest(BYTES_HASH)
    res_gf = expand_gf(res, GF16_HASH)

    if MAYO:
        msg_hash = res_gf[:m1]
    else:
        # Necessary to be compliant to KATs from C-Reference
        msg_hash = [res_gf[mi * l * r + j1 * l + i1] for mi in range(o) for i1 in range(l) for j1 in range(r)]

    return msg_hash


# XOF

def snova_xof(seed):
    # $\texttt{SNOVA{\_}SHAKE}$ public key expansion
    if aes:
        return nistrng.aesctr(seed, NUM_GEN_PUB_BYTES)
    else:
        # snova_shake
        blocks = (NUM_GEN_PUB_BYTES + 167) // 168
        res = bytearray()
        for i in range(blocks):
            blockseed = bytearray(seed)
            for j in range(8):
                blockseed.append((i >> (8 * j)) % 256)
            res += hashlib.shake_128(blockseed).digest(168)
        return bytes(res[:NUM_GEN_PUB_BYTES])


# Expand secret

def expand_T12(sk_data):
    # Generate the secret map $T_{12}$

    def gen_a_FqS(coefs):
        # Generate elements of $\mathbb{F}_{q}[S]$
        if not MAYO and coefs[l - 1] == 0:
            coefs[l - 1] = q - (coefs[0] if coefs[0] != 0 else 1)
        F = matrix(GF_q, l, l)
        for i in range(l):
            F += S**i * from_int(coefs[i])
        return F

    coef = []
    idx = 0
    i = 0
    while i < o * v * l:
        b = sk_data[idx]
        if q == 16:
            coef.append(b % 16)
            coef.append(b // 16)
            i += 2
        else:
            if b < (256 // q) * q:
                coef.append(b % q)
                i += 1
        idx += 1
    T12 = [gen_a_FqS(coef[l * i:]) for i in range(o * v)]

    # Convert to a single matrix
    T12m = matrix(GF_q, v * l, o * l)
    for ni in range(v):
        for nj in range(o):
            for i1 in range(l):
                for j1 in range(l):
                    T12m[ni * l + i1, nj * l + j1] = T12[ni * o + nj][i1, j1]
    return T12m


# Expand public

def convert_bytes_to_GF(data):
    # Expand public XOF data
    if q == 16:
        res = []
        for item in data:
            res.append(from_int(item % 16))
            res.append(from_int(item // 16))
        return res
    else:
        return [item % q for item in data]


def fixed_abq():
    NUM_ABQ = o * n_alpha * (r * (r + l) + 2 * l)
    abqdata = hashlib.shake_256(b'SNOVA_ABQ').digest(NUM_ABQ)
    return convert_bytes_to_GF(abqdata)


def expand_public_snova(seed):
    # Generate the random part of public key for $q=16$
    bindata = snova_xof(seed)
    data = convert_bytes_to_GF(bindata)

    idx = 0
    Pm11 = []
    for _ in range(m1):
        p11 = matrix(GF_q, v * l, v * l)
        for ni in range(v):
            for nj in range(v):
                for i1 in range(l):
                    for j1 in range(l):
                        p11[ni * l + i1, nj * l + j1] = data[idx]
                        idx += 1
        Pm11.append(p11)
    Pm12 = []
    for _ in range(m1):
        p12 = matrix(GF_q, v * l, o * l)
        for ni in range(v):
            for nj in range(o):
                for i1 in range(l):
                    for j1 in range(l):
                        p12[ni * l + i1, nj * l + j1] = data[idx]
                        idx += 1
        Pm12.append(p12)
    Pm21 = []
    for _ in range(m1):
        p21 = matrix(GF_q, o * l, v * l)
        for ni in range(o):
            for nj in range(v):
                for i1 in range(l):
                    for j1 in range(l):
                        p21[ni * l + i1, nj * l + j1] = data[idx]
                        idx += 1
        Pm21.append(p21)
    return Pm11, Pm12, Pm21, fixed_abq() if l < 4 else data[idx:]


def expand_public_mayo(seed):
    bindata = nistrng.aesctr(seed, NUM_GEN_PUB_BYTES)
    data = convert_bytes_to_GF(bindata)

    idx = 0
    Pm11 = [matrix(GF_q, v * l, v * l) for _ in range(m1)]
    Pm12 = [matrix(GF_q, v * l, o * l) for _ in range(m1)]
    Pm21 = [matrix(GF_q, o * l, v * l) for _ in range(m1)]
    for ni in range(v):
        for nj in range(ni, v):
            for mi in range(m1):
                Pm11[mi][ni, nj] = data[idx]
                idx += 1
    for ni in range(v):
        for nj in range(o):
            for mi in range(m1):
                Pm12[mi][ni, nj] = data[idx]
                idx += 1
    return Pm11, Pm12, Pm21, data[idx:]


def expand_public(seed):
    if MAYO:
        return expand_public_mayo(seed)
    else:
        return expand_public_snova(seed)


def compress_p22(pub22):
    # Pack the generated public key as bytes
    pk = bytearray()
    res = []
    if MAYO:
        for ni in range(o):
            for nj in range(ni, o):
                for mi in range(m1):
                    if ni == nj:
                        res.append(pub22[mi][ni, nj])
                    else:
                        res.append(pub22[mi][ni, nj] + pub22[mi][nj, ni])
    else:
        for mi in range(m1):
            for ni in range(o):
                for nj in range(o):
                    for i1 in range(l):
                        for j1 in range(l):
                            res.append(pub22[mi][ni * l + i1, nj * l + j1])
    pk += compress_gf(res, NUMGF_PK)
    return pk


def expand_p22(p22bytes):
    # Expand public key
    data = expand_gf(p22bytes, NUMGF_PK)
    idx = 0
    if MAYO:
        P22 = [matrix(GF_q, o * l, o * l) for _ in range(m1)]
        for ni in range(o):
            for nj in range(ni, o):
                for mi in range(m1):
                    P22[mi][ni, nj] = data[idx]
                    idx += 1
    else:
        P22 = []
        for _ in range(m1):
            pub22 = matrix(GF_q, o * l, o * l)
            for ni in range(o):
                for nj in range(o):
                    for i1 in range(l):
                        for j1 in range(l):
                            pub22[ni * l + i1, nj * l + j1] = data[idx]
                            idx += 1
            P22.append(pub22)
    return P22


def decode_sig(sig_bytes):
    # Decode sig
    if MAYO:
        salt = sig_bytes[-24:]
    else:
        salt = sig_bytes[-16:]
    gfsig = expand_gf(sig_bytes[:-16], n * l * r)
    if len(gfsig) < n * l * r:
        raise Exception('Verify failed.')
    if MAYO:
        sig = matrix(GF_q, n * l, r, lambda i, j: gfsig[j * n * l + i])
    else:
        sig = matrix(GF_q, n * l, r, lambda i, j: gfsig[i * r + j])

    return sig, salt


# Generate E matrix

def gen_E(abqdata):
    if MAYO:
        import mayo_emat
        emat = [[[matrix(GF_q, m1, r**2, lambda mj, i1: from_int(mayo_emat.mayo1[mi][0][0][mj][i1]))
                  for a in range(l)] for b in range(l)] for mi in range(m1)]
        return emat

    if RANDOM_E:
        NUM_E = m1 * m2 * r**2 * l**2
        e_bytes = hashlib.shake_256(b'SNOVA_E').digest(NUM_E)
        e_data = convert_bytes_to_GF(e_bytes)
        emat = [[[matrix(GF_q, m2, r**2, lambda i, j: e_data[(((a * l + b) * m1 + mi) * m2 + i) * r**2 + j])
                 for a in range(l)] for b in range(l)] for mi in range(m1)]
        return emat

    def create_AB(data, r1, r2):
        # Generate invertible matrices
        M = matrix(GF_q, r1, r2, lambda i, j: data[i * r2 + j])
        if l == r1 and l == r2:
            f1 = 1
            while M.det() == 0 and f1 < q:
                M += from_int(f1) * S
                f1 += 1
            if f1 == q:
                raise Exception('f1 == q')
        return M

    def nonzero_q(data):
        fq = [to_int(data[i]) for i in range(l)]
        if fq[l - 1] == 0:
            fq[l - 1] = q - (fq[0] if fq[0] != 0 else 1)
        return [from_int(fq[i]) for i in range(l)]

    A = [create_AB(abqdata[i * r**2:], r, r) for i in range(o * n_alpha)]
    B = [create_AB(abqdata[o * n_alpha * r**2 + i * l * r:], r, l) for i in range(o * n_alpha)]
    q1 = [nonzero_q(abqdata[o * n_alpha * r * (r + l) + i * l:]) for i in range(o * n_alpha)]
    q2 = [nonzero_q(abqdata[o * n_alpha * r * (r + l) + o * n_alpha * l + i * l:]) for i in range(o * n_alpha)]

    emat = [[[matrix(GF_q, m2, r**2) for a in range(l)] for b in range(l)] for mi in range(m1)]
    for mi in range(o):
        for alpha in range(n_alpha):
            mia = mi * n_alpha + alpha
            mi_prime = (mi + alpha) % m1
            for a in range(l):
                for b in range(l):
                    for i1 in range(l):
                        for i2 in range(r):
                            for j1 in range(r):
                                for j2 in range(r):
                                    emat[mi_prime][a][b][mi * l * r + i1 * r + i2, j1 * r + j2] += \
                                        A[mia][i2, j1] * q1[mia][a] * q2[mia][b] * B[mia][j2, i1]
    return emat

################################################################

# API functions
# Generate keypair from seed


def genkeys(seed):
    # Generate Public key
    if MAYO:
        sk_seed = seed[:24]
        sk = sk_seed
        sk_data = hashlib.shake_256(sk_seed).digest(2 * o * v * l)  # Overdimensioned
        pk_seed = sk_data[:16]
        t12data = sk_data[16:]
    else:
        sk_seed = seed[16:]
        sk = seed
        pk_seed = seed[:16]
        t12data = hashlib.shake_256(sk_seed).digest(2 * o * v * l)  # Overdimensioned

    T12 = expand_T12(t12data)
    P11, P12, P21, _ = expand_public(pk_seed)

    P22 = [-(T12.transpose() * (P11[mi] * T12 + P12[mi]) + P21[mi] * T12) for mi in range(m1)]

    return sk, pk_seed + compress_p22(P22)


# Sign message

def sign(sk, msg, salt):
    if MAYO:
        sk_seed = sk[:24]
        sk_data = hashlib.shake_256(sk_seed).digest(2 * o * v * l)  # Overdimensioned
        pk_seed = sk_data[:16]
        t12data = sk_data[16:]
        m_g = m1
    else:
        sk_seed = sk[16:]
        pk_seed = sk[:16]
        t12data = hashlib.shake_256(sk_seed).digest(2 * o * v * l)  # Overdimensioned
        m_g = m2

    T12 = expand_T12(t12data)

    P11, P12, P21, abqdata = expand_public(pk_seed)
    emat = gen_E(abqdata)

    # Expand private key
    F12 = []
    F21 = []
    for mi in range(m1):
        F12.append(P11[mi] * T12 + P12[mi])
        F21.append(T12.transpose() * P11[mi] + P21[mi])

    num_sign = 0
    while True:
        num_sign += 1
        if num_sign == 255:
            raise Exception('signing failed')

        # Assign values to vinegar variables
        # Vinegar from sk and salt

        if MAYO:
            dgst = hashlib.shake_256(msg).digest(32)
            salt_byte = hashlib.shake_256(dgst + salt + sk_seed).digest(24)
            msg_hash = hash_combined(msg, pk_seed, salt_byte)
            vinegar_byte = hashlib.shake_256(dgst + salt_byte + sk_seed + bytes([num_sign - 1])).digest(BYTES_GF(r * n))
            vinegar_gf = expand_gf(vinegar_byte, n * l * r)
            vinegar = matrix(GF_q, v * l, r, lambda i, j: vinegar_gf[j * v + i])
            rvec = vector([vinegar_gf[v * l * r + j1 * o * l + i1] for j1 in range(r) for i1 in range(o * l)])
        else:
            salt_byte = salt
            msg_hash = hash_combined(msg, pk_seed, salt_byte)
            v_state = hashlib.shake_256()
            v_state.update(sk_seed)
            dgst = hashlib.shake_256(msg).digest(64)
            v_state.update(dgst)
            v_state.update(salt)
            v_state.update(num_sign.to_bytes(1))
            vinegar_byte = v_state.digest(BYTES_GF(v * l * r))
            vinegar_gf = expand_gf(vinegar_byte, v * l * r)
            vinegar = matrix(GF_q, v * l, r, lambda i, j: vinegar_gf[i * r + j])
            rvec = vector(GF_q, m2)

        # Compute the vinegar part of the central map
        p_vin = [[[vinegar.transpose() * S_times_v**b * P11[mi] * S_times_v**a * vinegar
                   for a in range(l)] for b in range(l)] for mi in range(m1)]

        # Apply emulsifier
        res = vector(GF_q, m_g)
        for mi in range(m1):
            for a in range(l):
                for b in range(l):
                    vin_vec = vector([p_vin[mi][a][b][i1, j1] for i1 in range(r) for j1 in range(r)])
                    res += emat[mi][a][b] * vin_vec
        F_vv = list(res)

        # Get msg vinegar part
        msg_vv = vector([msg_hash[idx] - F_vv[idx] for idx in range(m_g)])

        # Compute the coefficient matrix of the oil variable
        # compute the coefficients of Xo and put into gauss matrix
        F21_ab = [[[S_times_o**b * F21[mi] * S_times_v**a * vinegar for a in range(l)] for b in range(l)] for mi in range(m1)]
        F12_ab = [[[vinegar.transpose() * S_times_v**b * F12[mi] * S_times_o**a for a in range(l)] for b in range(l)] for mi in range(m1)]

        A = matrix(GF_q, m_g, m2)
        for ti1 in range(m_g):
            for tj1 in range(o * l):
                for tj2 in range(r):
                    for mi_prime in range(m1):
                        for a in range(l):
                            for b in range(l):
                                for k1 in range(r):
                                    A[ti1, tj2 * o * l + tj1] += \
                                        emat[mi_prime][a][b][ti1, tj2 * r + k1] * F21_ab[mi_prime][a][b][tj1, k1] \
                                        + emat[mi_prime][a][b][ti1, k1 * r + tj2] * F12_ab[mi_prime][a][b][k1, tj1]

        try:
            solution = A.solve_right(msg_vv - A * rvec) + rvec
            solution = [solution[tj2 * o * l + tj1] for tj1 in range(o * l) for tj2 in range(r)]
        except ValueError:
            # print('no solution')
            continue

        sol_mat = matrix(GF_q, o * l, r, lambda i, j: solution[i * r + j])
        vinegar += T12 * sol_mat
        sig_gf = [vinegar[mi * l + i1, j1] for mi in range(v) for i1 in range(l) for j1 in range(r)]
        sig_gf += solution

        break

    if MAYO:
        sig_gf = [sig_gf[i * r + j] for j in range(r) for i in range(n * l)]

    return compress_gf(sig_gf, n * l * r) + salt_byte


# Verify signature of message

def verify(pk, sig_bytes, msg):
    # Verify signature of message

    # Decode sig
    sig, salt = decode_sig(sig_bytes)

    # Expand pubkey
    pk_seed = pk[:16]
    P11, P12, P21, abqdata = expand_public(pk_seed)
    P22 = expand_p22(pk[16:])
    P = [matrix.block([[P11[mi], P12[mi]], [P21[mi], P22[mi]]]) for mi in range(m1)]

    # Whip-up signature
    p_sig = [[[sig.transpose() * S_times_n**b * P[mi] * S_times_n**a * sig
               for a in range(l)] for b in range(l)] for mi in range(m1)]

    # Apply emulsifier
    emat = gen_E(abqdata)

    res = vector(GF_q, emat[0][0][0].nrows())
    for mi in range(m1):
        for a in range(l):
            for b in range(l):
                sig_vec = vector([p_sig[mi][a][b][i1, j1] for i1 in range(r) for j1 in range(r)])
                res += emat[mi][a][b] * sig_vec
    sig_hash = list(res)

    # Check against expected hash
    msg_hash = hash_combined(msg, pk_seed, salt)

    if msg_hash != sig_hash:
        raise Exception('Verify failed')


################################################################


# Generate KATs
entropy_input = bytearray(48)
for i in range(48):
    entropy_input[i] = i
drbg = nistrng.rng(entropy_input)
srbdg = nistrng.rng(entropy_input)

print('# SNOVA', v, o, q, l, aes, r, m1)
print()

for count in range(1):
    seed = drbg.random_bytes(48)
    srbdg = nistrng.rng(seed)

    if MAYO:
        keygen_seed = srbdg.random_bytes(24)
        salt = srbdg.random_bytes(24)
    else:
        keygen_seed = srbdg.random_bytes(48)
        salt = srbdg.random_bytes(16)

    mlen = 33 * (count + 1)
    msg = drbg.random_bytes(mlen)

    try:
        # Keygen
        sk, pk = genkeys(keygen_seed)

        print('count =', count)
        print('seed =', seed.hex().upper())
        print('mlen =', mlen)
        print('msg =', msg.hex().upper())
        print('pk =', pk.hex().upper())
        print('sk =', sk.hex().upper())

        # Sign
        sig = sign(sk, msg, salt)
        sm = sig + msg

        print('smlen =', len(sm))
        print('sm =', sm.hex().upper())
        print()

        # Verify
        verify(pk, sig, msg)

    except Exception as exc:
        traceback.print_exc()
        quit()
