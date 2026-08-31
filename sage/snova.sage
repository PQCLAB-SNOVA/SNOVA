# SPDX-License-Identifier: MIT
#
# SNOVA implementation in SageMath
#
# Copyright (c) 2026 SNOVA TEAM

from hashlib import shake_128, shake_256
import math
import traceback

try:
    import nistrng
except:
    print('Error importing nistrng')
    print('Try: export PYTHONPATH=`pwd`')
    quit()

################################################################

# SNOVA parameters

v = 27
o = 5
q = 16
l = 4
r = l
m1 = 5
aes = False

################################################################

n = v + o
n_alpha = l * r + 2 * r

HASH_PK = l > 2
ABQ_ALG2 = True

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

if q == 11:
    Q_A = 0
    Q_B = 3
    Q_C = 6
    PACK_GF = 16
    PACK_BYTES = 7

elif q == 13:
    Q_A = 2
    Q_B = 11
    Q_C = 3
    PACK_GF = 15
    PACK_BYTES = 7

elif q == 16:
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


GF16_HASH = o * l * r
BYTES_HASH = BYTES_GF(GF16_HASH)

NUM_GEN_PUB_GF = m1 * (v * v + 2 * v * o) * l**2 + o * n_alpha * (r * (r + l) + 2 * l)
NUMGF_PK = m1 * o * l * (o * l)

if q == 16:
    NUM_GEN_PUB_BYTES = math.ceil((NUM_GEN_PUB_GF + 1) / 2)
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

S_times_v = matrix.block_diagonal([S for _ in range(v)])
S_times_n = matrix.block_diagonal([S for _ in range(n)])


# Utils

def expand_gf(data, num, check=False):
    # Convert bytes to elements of $\mathbb{F}_{q}$
    res = []
    idx = 0
    while idx < len(data):
        sum = 0
        for i in range(min(PACK_BYTES, len(data) - idx)):
            sum += int(data[idx + i]) * int(256**i)
        idx += PACK_BYTES
        for i in range(PACK_GF):
            res.append(from_int(int(sum) % q))
            sum = sum // q
        if check and sum:
            raise Exception('expand_gf illegal data')
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
    res = shake_256(pk_seed + shake_256(msg).digest(64) + salt).digest(BYTES_HASH)
    res_gf = expand_gf(res, GF16_HASH)

    # Reorder, necessary to be compliant to KATs from C-Reference
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
            res += shake_128(blockseed).digest(168)
        return bytes(res[:NUM_GEN_PUB_BYTES])


# Expand secret

def expand_T12(seed):
    # Generate the secret map $T_{12}$

    def gen_a_FqS(coefs):
        # Generate elements of $\mathbb{F}_{q}[S]$
        F = matrix(GF_q, l, l)
        for i in range(l):
            F += S**i * from_int(coefs[i])
        return F

    sk_data = shake_256(seed).digest(2 * o * v * l)  # Overdimensioned
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


def expand_public(seed):
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
    return Pm11, Pm12, Pm21


def compress_p22(pub22):
    # Pack the generated public key as bytes
    pk = bytearray()
    res = []
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
    data = expand_gf(p22bytes, NUMGF_PK, check=True)
    P22 = []
    idx = 0
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


def gen_ABQ():
    # Generate public ABQ

    def create_AB(data, r1, r2):
        # Improve public matrices
        M = matrix(GF_q, r1, r2, lambda i, j: data[i * r2 + j])
        if ABQ_ALG2 and l == r1 and l == r2:
            f1 = 1
            while M.det() == 0 and f1 < q:
                M += from_int(f1) * S
                f1 += 1
            if f1 == q:
                raise Exception('f1 == q')
        return M

    NUM_ABQ = o * n_alpha * (r * (r + l) + 2 * l)
    seed = b'SNOVA_ABQ'
    if l == 2:
        if o == 17:
            seed = b'SNOVA_ABQ_2'
        else:
            raise Exception('Unsupported ', l, o)
    abqbytes = shake_256(seed).digest(NUM_ABQ)
    abqdata = convert_bytes_to_GF(abqbytes)

    A = [create_AB(abqdata[i * r**2:], r, r) for i in range(o * n_alpha)]
    B = [create_AB(abqdata[o * n_alpha * r**2 + i * l * r:], r, l) for i in range(o * n_alpha)]
    q1 = [abqdata[o * n_alpha * r * (r + l) + i * l:] for i in range(o * n_alpha)]
    q2 = [abqdata[o * n_alpha * r * (r + l) + o * n_alpha * l + i * l:] for i in range(o * n_alpha)]

    return A, B, q1, q2


################################################################

# API functions
# Generate keypair from seed


def genkeys(seed):
    # Generate Public key
    sk_seed = seed[16:]
    T12 = expand_T12(sk_seed)

    pk_seed = seed[:16]
    P11, P12, P21 = expand_public(pk_seed)

    P22 = [-(T12.transpose() * (P11[mi] * T12 + P12[mi]) + P21[mi] * T12) for mi in range(m1)]

    pk = pk_seed + compress_p22(P22)
    if HASH_PK:
        return seed + shake_256(pk).digest(48), pk
    else:
        return seed, pk


# Sign message

def sign(sk, msg, salt):
    # Sign message
    sk_seed = sk[16:48]
    T12 = expand_T12(sk_seed)

    pk_seed = sk[:16]
    P11, P12, P21 = expand_public(pk_seed)

    # Expand private key
    F12 = []
    F21 = []
    for mi in range(m1):
        F12.append(P11[mi] * T12 + P12[mi])
        F21.append(T12.transpose() * P11[mi] + P21[mi])

    A, B, q1, q2 = gen_ABQ()

    Q1 = []
    Q2 = []

    for idx in range(o * n_alpha):
        q1mat = matrix(GF_q, l, l)
        q2mat = matrix(GF_q, l, l)
        for ab in range(l):
            q1mat += q1[idx][ab] * S**ab
            q2mat += q2[idx][ab] * S**ab
        Q1.append(q1mat)
        Q2.append(q2mat)

    msg_hash = hash_combined(msg, sk[48:] if HASH_PK else pk_seed, salt)

    num_sign = 0
    while True:
        num_sign += 1
        if num_sign == 255:
            raise Exception('signing failed')

        # Assign values to vinegar variables
        # Vinegar from sk and salt

        v_state = shake_256(sk_seed + shake_256(msg).digest(64) + salt + num_sign.to_bytes(1))
        vinegar_byte = v_state.digest(BYTES_GF(v * l * r))
        vinegar_gf = expand_gf(vinegar_byte, v * l * r)
        vinegar = matrix(GF_q, v * l, r, lambda i, j: vinegar_gf[i * r + j])

        # Compute the vinegar part of the central map
        p_vin = [[[vinegar.transpose() * S_times_v**a * P11[mi] * S_times_v**b * vinegar
                   for b in range(l)] for a in range(l)] for mi in range(m1)]

        # Apply emulsifier
        F_vv_mat = [matrix(GF_q, r, l) for mi in range(o)]
        for mi in range(o):
            for alpha in range(n_alpha):
                mia = mi * n_alpha + alpha
                mi_prime = (mi + alpha) % m1
                pqq = matrix(GF_q, r, r)
                for a in range(l):
                    for b in range(l):
                        pqq += q1[mia][a] * p_vin[mi_prime][a][b] * q2[mia][b]
                F_vv_mat[mi] += A[mia] * pqq * B[mia]
        F_vv = [F_vv_mat[mi][i1][j1] for mi in range(o) for j1 in range(l) for i1 in range(r)]

        # Get msg vinegar part
        msg_vv = vector([msg_hash[idx] - F_vv[idx] for idx in range(o * l * r)])

        # Compute the coefficient matrix of the oil variable
        # compute the coefficients of Xo and put into gauss matrix
        gauss = matrix(GF_q, o * l * r, o * l * r)
        for mi in range(o):
            for alpha in range(n_alpha):
                mia = mi * n_alpha + alpha
                mi_prime = (mi + alpha) % m1

                Q1_v = matrix.block_diagonal([Q1[mia] for _ in range(v)])
                Q1_o = matrix.block_diagonal([Q1[mia] for _ in range(o)])
                Q2_v = matrix.block_diagonal([Q2[mia] for _ in range(v)])
                Q2_o = matrix.block_diagonal([Q2[mia] for _ in range(o)])

                H_1 = Q1_o * (F21[mi_prime] * Q2_v * vinegar) * B[mia]
                H_2 = A[mia] * (vinegar.transpose() * Q1_v * F12[mi_prime]) * Q2_o

                for i0 in range(o):
                    for i1 in range(l):
                        for i2 in range(r):
                            for j1 in range(l):
                                for j2 in range(r):
                                    val = H_1[i0 * l + j1, i1] * A[mia][i2, j2] \
                                        + H_2[i2, i0 * l + j1] * B[mia][j2, i1]
                                    gauss[mi * l * r + i1 * r + i2, i0 * l * r + j1 * r + j2] += val

        if gauss.det() == 0:
            # print('Try with another vinegar value', num_sign)
            continue

        solution = gauss.solve_right(msg_vv)

        sol_mat = matrix(GF_q, o * l, r, lambda i, j: solution[i * r + j])
        vinegar += T12 * sol_mat
        sig_gf = [vinegar[mi * l + i1, j1] for mi in range(v) for i1 in range(l) for j1 in range(r)]
        sig_gf += solution

        break

    return compress_gf(sig_gf, n * l * r) + salt


# Verify signature of message

def verify(pk, sig_bytes, msg):
    # Verify signature of message

    # Decode sig
    salt = sig_bytes[-16:]
    gfsig = expand_gf(sig_bytes[:-16], n * l * r, check=True)
    if len(gfsig) < n * l * r:
        raise Exception('Verify failed.')
    sig = matrix(GF_q, n * l, r, lambda i, j: gfsig[i * r + j])

    # Expand pubkey
    pk_seed = pk[:16]
    P11, P12, P21 = expand_public(pk_seed)
    P22 = expand_p22(pk[16:])
    P = [matrix.block([[P11[mi], P12[mi]], [P21[mi], P22[mi]]]) for mi in range(m1)]
    A, B, q1, q2 = gen_ABQ()

    # Whip-up signature
    p_sig = [[[sig.transpose() * S_times_n**b * P[mi] * S_times_n**a * sig
               for a in range(l)] for b in range(l)] for mi in range(m1)]

    # Apply emulsifier
    temp = [matrix(GF_q, r, l) for mi in range(o)]
    for mi in range(o):
        for alpha in range(n_alpha):
            mia = mi * n_alpha + alpha
            mi_prime = (mi + alpha) % m1
            pqq = matrix(GF_q, r, r)
            for a in range(l):
                for b in range(l):
                    pqq += q1[mia][a] * p_sig[mi_prime][a][b] * q2[mia][b]
            temp[mi] += A[mia] * pqq * B[mia]
    sig_hash = [temp[mi][i1][j1] for mi in range(o) for j1 in range(l) for i1 in range(r)]

    # Check against expected hash
    msg_hash = hash_combined(msg, shake_256(pk).digest(48) if HASH_PK else pk_seed, salt)

    if msg_hash != sig_hash:
        raise Exception('Verify failed')


################################################################

# Generate KATs
entropy_input = bytearray(48)
for i in range(48):
    entropy_input[i] = i
drbg = nistrng.rng(entropy_input)

m2 = r * l * o
print('# SNOVA', v, o, q, l, r, m1, m2, 'AES' if aes else 'SHAKE', n_alpha)
print()

for count in range(1):
    seed = drbg.random_bytes(48)
    srbdg = nistrng.rng(seed)

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
