from .aes256_ctr_drbg import AES256_CTR_DRBG
from . import pyaes

def rng(seed):
    return AES256_CTR_DRBG(bytes(seed))

def aesctr(seed, num, iv=0):
    res = bytearray()
    # Qorks for QR-UOV as iv < 255
    block_i = iv << 64
    cipher = pyaes.AESModeOfOperationECB(bytes(seed))
    while len(res) < num:
        blockseed = bytearray()
        for j in range(16):
            blockseed.append((block_i >> (8 * (15 - j))) % 256)
        res += cipher.encrypt(bytes(blockseed))
        block_i += 1
    return bytes(res[:num])
