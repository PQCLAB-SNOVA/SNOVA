# SPDX-License-Identifier: MIT
#
# SNOVA implementation in SageMath
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

v = 24
o = 5
q = 23
l = 4
aes = False

# RectSNOVA parameters

r = l
m1 = o

################################################################

n_alpha = r * r + r
n = v + o

ASYMMETRIC_PUBMAT = q == 16
ROUND2_KAT = l == r and q == 16

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


GF16_HASH = o * l * r
BYTES_HASH = BYTES_GF(GF16_HASH)

if ASYMMETRIC_PUBMAT:
    NUM_GEN_PUB_GF = m1 * (v * v + 2 * v * o) * l**2 + o * n_alpha * (r * (r + l) + 2 * l)
    NUMGF_PK = m1 * o * l * (o * l)
else:
    NUM_GEN_PUB_GF = m1 * (v * (v + 1) // 2 + v * o) * l**2 + o * n_alpha * (r * (r + l) + 2 * l)
    NUMGF_PK = m1 * o * l * (o * l + 1) // 2

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
    state.update(pk_seed)
    if ROUND2_KAT:
        dgst = hashlib.shake_256(msg).digest(64)
        state.update(dgst)
    else:
        state.update(msg)
    state.update(salt)
    res = state.digest(BYTES_HASH)
    res_gf = expand_gf(res, GF16_HASH)

    # Necessary to be compliant to KATs from C-Reference
    msg_hash = [res_gf[mi * l * r + j1 * l + i1] for mi in range(o) for i1 in range(l) for j1 in range(r)]

    return msg_hash, res


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

def expand_T12(seed):
    # Generate the secret map $T_{12}$

    def gen_a_FqS(coefs):
        # Generate elements of $\mathbb{F}_{q}[S]$
        if coefs[l - 1] == 0:
            coefs[l - 1] = q - (coefs[0] if coefs[0] != 0 else 1)
        F = matrix(GF_q, l, l)
        for i in range(l):
            F += S**i * from_int(coefs[i])
        return F

    sk_data = hashlib.shake_256(seed).digest(2 * o * v * l)  # Overdimensioned
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


def expand_public_sym(seed):
    # Generate the random part of public key for odd $q$
    bindata = snova_xof(seed)
    data = convert_bytes_to_GF(bindata)

    idx = 0
    Pm11 = []
    Pm12 = []
    Pm21 = []
    for _ in range(m1):
        p11 = matrix(GF_q, v * l, v * l)
        p12 = matrix(GF_q, v * l, o * l)
        p21 = matrix(GF_q, o * l, v * l)
        for ni in range(v):
            for i1 in range(l):
                for j1 in range(i1, l):
                    p11[ni * l + i1, ni * l + j1] = data[idx]
                    p11[ni * l + j1, ni * l + i1] = data[idx]
                    idx += 1
            for nj in range(ni + 1, v):
                for i1 in range(l):
                    for j1 in range(l):
                        p11[ni * l + i1, nj * l + j1] = data[idx]
                        p11[nj * l + j1, ni * l + i1] = data[idx]
                        idx += 1
            for nj in range(o):
                for i1 in range(l):
                    for j1 in range(l):
                        p12[ni * l + i1, nj * l + j1] = data[idx]
                        p21[nj * l + j1, ni * l + i1] = data[idx]
                        idx += 1
        Pm11.append(p11)
        Pm12.append(p12)
        Pm21.append(p21)
    return Pm11, Pm12, Pm21, fixed_abq() if l < 4 else data[idx:]


def expand_public_asym(seed):
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


def expand_public(seed):
    if ASYMMETRIC_PUBMAT:
        return expand_public_asym(seed)
    else:
        return expand_public_sym(seed)


def compress_p22(pub22):
    # Pack the generated public key as bytes
    pk = bytearray()
    res = []
    for mi in range(m1):
        for ni in range(o):
            if ASYMMETRIC_PUBMAT:
                for nj in range(o):
                    for i1 in range(l):
                        for j1 in range(l):
                            res.append(pub22[mi][ni * l + i1, nj * l + j1])
            else:
                for i1 in range(l):
                    for j1 in range(i1, l):
                        res.append(pub22[mi][ni * l + i1, ni * l + j1])
                    for nj in range(ni + 1, o):
                        for j1 in range(l):
                            res.append(pub22[mi][ni * l + i1, nj * l + j1])
    pk += compress_gf(res, NUMGF_PK)
    return pk


def expand_p22(p22bytes):
    # Expand public key
    data = expand_gf(p22bytes, NUMGF_PK)
    P22 = []
    idx = 0
    for _ in range(m1):
        pub22 = matrix(GF_q, o * l, o * l)
        for ni in range(o):
            if ASYMMETRIC_PUBMAT:
                for nj in range(o):
                    for i1 in range(l):
                        for j1 in range(l):
                            pub22[ni * l + i1, nj * l + j1] = data[idx]
                            idx += 1
            else:
                for i1 in range(l):
                    for j1 in range(i1, l):
                        pub22[ni * l + i1, ni * l + j1] = data[idx]
                        pub22[ni * l + j1, ni * l + i1] = data[idx]
                        idx += 1
                    for nj in range(ni + 1, o):
                        for j1 in range(l):
                            pub22[ni * l + i1, nj * l + j1] = data[idx]
                            pub22[nj * l + j1, ni * l + i1] = data[idx]
                            idx += 1
        P22.append(pub22)
    return P22


def gen_ABQ(abqdata):
    # Generate public ABQ from XOF data

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

    return A, B, q1, q2


################################################################

# API functions
# Generate keypair from seed


def genkeys(seed):
    # Generate Public key
    sk_seed = seed[16:]
    T12 = expand_T12(sk_seed)

    pk_seed = seed[:16]
    P11, P12, P21, _ = expand_public(pk_seed)

    P22 = [-(T12.transpose() * (P11[mi] * T12 + P12[mi]) + P21[mi] * T12) for mi in range(m1)]

    return seed, pk_seed + compress_p22(P22)


# Sign message

def sign(sk, msg, salt):
    # Sign message
    sk_seed = sk[16:]
    T12 = expand_T12(sk_seed)

    pk_seed = sk[:16]
    P11, P12, P21, abqdata = expand_public(pk_seed)

    # Expand private key
    F12 = []
    F21 = []
    for mi in range(m1):
        F12.append(P11[mi] * T12 + P12[mi])
        F21.append(T12.transpose() * P11[mi] + P21[mi])

    A, B, q1, q2 = gen_ABQ(abqdata)

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

    msg_hash, msg_bytes = hash_combined(msg, pk_seed, salt)

    num_sign = 0
    while True:
        num_sign += 1
        if num_sign == 255:
            raise Exception('signing failed')

        # Assign values to vinegar variables
        # Vinegar from sk and salt

        v_state = hashlib.shake_256()
        v_state.update(sk_seed)
        if ROUND2_KAT:
            dgst = hashlib.shake_256(msg).digest(64)
            v_state.update(dgst)
            v_state.update(salt)
        else:
            v_state.update(msg_bytes)
        v_state.update(num_sign.to_bytes(1))
        vinegar_byte = v_state.digest(BYTES_GF(v * l * r))
        vinegar_gf = expand_gf(vinegar_byte, v * l * r)
        vinegar = matrix(GF_q, v * l, r, lambda i, j: vinegar_gf[i * r + j])

        # Compute the vinegar part of the central map
        p_vin = [[[vinegar.transpose() * S_times_v**b * P11[mi] * S_times_v**a * vinegar
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
                        pqq += q1[mia][a] * p_vin[mi_prime][a][b] * q2[mia][b]
                temp[mi] += A[mia] * pqq * B[mia]
        F_vv = [temp[mi][i1][j1] for mi in range(o) for j1 in range(l) for i1 in range(r)]

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

                temp1 = Q1_o * (F21[mi_prime] * Q2_v * vinegar) * B[mia]
                temp2 = A[mia] * (vinegar.transpose() * Q1_v * F12[mi_prime]) * Q2_o

                for idx in range(o):
                    for ti1 in range(l):
                        for ti2 in range(r):
                            for tj1 in range(l):
                                for tj2 in range(r):
                                    val = temp1[idx * l + tj1, ti1] * A[mia][ti2, tj2] \
                                        + temp2[ti2, idx * l + tj1] * B[mia][tj2, ti1]
                                    gauss[mi * l * r + ti1 * r + ti2, idx * l * r + tj1 * r + tj2] += val

        try:
            solution = gauss.solve_right(msg_vv)
        except ValueError:
            # print('no solution')
            continue

        sol_mat = matrix(GF_q, o * l, r, lambda i, j: solution[i * r + j])
        vinegar += T12 * sol_mat
        sig_gf = [vinegar[mi * l + i1, j1] for mi in range(v) for i1 in range(l) for j1 in range(r)]
        sig_gf += solution

        # Check for symmetric signature
        ok = True
        if l == r and not ASYMMETRIC_PUBMAT:
            for idx in range(n):
                is_sym = True
                for i1 in range(l):
                    for j1 in range(i1 + 1, l):
                        if sig_gf[idx * l * r + i1 * r + j1] != sig_gf[idx * l * r + j1 * r + i1]:
                            is_sym = False
                if is_sym:
                    ok = False
        if ok:
            break

    return compress_gf(sig_gf, n * l * r) + salt


# Verify signature of message

def verify(pk, sig_bytes, msg):
    # Verify signature of message

    # Decode sig
    salt = sig_bytes[-16:]
    gfsig = expand_gf(sig_bytes[:-16], n * l * r)
    if len(gfsig) < n * l * r:
        raise Exception('Verify failed.')
    sig = matrix(GF_q, n * l, r, lambda i, j: gfsig[i * r + j])

    # Check for symmetric signature
    if l == r and not ASYMMETRIC_PUBMAT:
        for idx in range(n):
            is_sym = True
            for i1 in range(l):
                for j1 in range(i1 + 1, l):
                    if sig[idx * l + i1, j1] != sig[idx * l + j1, i1]:
                        is_sym = False
            if is_sym:
                raise Exception('Verify failed!')

    # Expand pubkey
    pk_seed = pk[:16]
    P11, P12, P21, abqdata = expand_public(pk_seed)
    P22 = expand_p22(pk[16:])
    P = [matrix.block([[P11[mi], P12[mi]], [P21[mi], P22[mi]]]) for mi in range(m1)]
    A, B, q1, q2 = gen_ABQ(abqdata)

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
    msg_hash, _ = hash_combined(msg, pk_seed, salt)

    if msg_hash != sig_hash:
        raise Exception('Verify failed')


################################################################

# Generate KATs
entropy_input = bytearray(48)
for i in range(48):
    entropy_input[i] = i
drbg = nistrng.rng(entropy_input)

print('# SNOVA', v, o, q, l, 'AES' if aes else 'SHAKE', r, m1)
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
