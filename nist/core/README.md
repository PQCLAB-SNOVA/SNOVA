# SNOVA reference implementation

Requirements: GNU make and a C11 compiler (gcc/clang).
The default optimized build (`ARCH=x86_avx2`) targets x86-64 with AVX2 + GFNI
(Ice Lake+ / Zen 4+). On x86-64 without GFNI use `ARCH=x86_avx2 GFNI=0`; on any other
platform use the portable reference backend `ARCH=ref`. All backends produce
bit-identical KAT output.

## Parameter sets: the complete 7-tuple

A SNOVA parameter set is the 7-tuple `(v, o, q, l, r, m1, m2)`, and its name is those
seven values in that order  -  e.g. `27_5_16_4_4_5_80` means v=27, o=5, q=16, l=4, r=4,
m1=5, m2=80. Every build passes all seven:

    SNOVA_V=<v> SNOVA_O=<o> SNOVA_Q=<q> SNOVA_L=<l> SNOVA_R=<r> SNOVA_M1=<m1> SNOVA_M2=<m2>

There is no short form. `SNOVA_R` and `SNOVA_M1` are not optional extras for "some"
sets: m1 is an independent parameter that differs from floor(o*r/l) in several of the
shipped sets, so leaving it out silently builds a *different* scheme that will not match
the reference KAT. `SNOVA_M2` is derived (m2 = o*l*r) and can never be overridden  - 
supplying it only feeds a compile-time `_Static_assert` cross-check that the tuple you
typed is self-consistent, which is why every example below carries it.

`snova_config.src` lists the frozen complete 7-tuple of every shipped parameter set,
in both XOF variants.

Every submitted parameter set uses q=16. The implementation is written for general q,
so the sources also contain the prime-field (odd q) branches, including their own AVX2
kernels; these compile and are byte-exact against the reference backend, but no submitted
set selects them. They are retained deliberately for future parameter choices and are not
dead code.

## Build + NIST KAT generation

    make PQCgenKAT SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80 ARCH=ref   # 27_5_16_4_4_5_80, portable reference
    make PQCgenKAT SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80            # 27_5_16_4_4_5_80, optimized (default; needs AVX2+GFNI)
    make PQCgenKAT SNOVA_V=29 SNOVA_O=3 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=8 SNOVA_M1=5 SNOVA_M2=96            # 29_3_16_4_8_5_96 (here m1=5, not floor(o*r/l)=6)

- Output: `build/<config>/PQCsignKAT_<ALGNAME>.{req,rsp}`, where `<ALGNAME>` is
  `SNOVA_<v>_<o>_<q>_<l>_<r>_<m1>_<m2>` (with `_AES` appended for the AES variant).
- `AES=1` selects the AES public-key expansion variant; `AES=0` (the default) selects
  SHAKE. `AES` is the only XOF knob  -  `PK_EXPAND_SHAKE` is derived from it and is
  ignored if given on the command line.

## Building every parameter set

This folder builds one parameter set at a time. To build them all, use the expander one
level up: `snova_config.src` there lists every set, and the Makefile beside it gives each its
own folder with the build flags already frozen in.

    cd ..
    make                     # every set listed in snova_config.src, both XOF variants
    make SNOVA_29_3_16_4_8_5_96   # just one set (the folder name is the target)

## KAT: reference digests + self-verification

The authoritative KAT reference is `ref_kat/KAT_DIGESTS.sha256`  -  the SHA-256 of each
parameter's response file (both the SHAKE and `_AES` variants). Full `PQCsignKAT_*.req/.rsp`
files (~1 MB per set) are intentionally not shipped, to keep the package small; regenerate
and verify them yourself:

    make kat-verify SNOVA_V=27 SNOVA_O=5 SNOVA_Q=16 SNOVA_L=4 SNOVA_R=4 SNOVA_M1=5 SNOVA_M2=80   # one set: regenerate .rsp + check sha256 vs ref_kat/
    cd .. && make            # every set, via the expander one level up

`make kat-verify` regenerates `PQCsignKAT_<ALGNAME>.{req,rsp}` under `build/<config>/` and
checks the .rsp SHA-256 against `ref_kat/KAT_DIGESTS.sha256`; a MATCH confirms byte-exact KAT.
All backends (x86_avx2 GFNI / no-GFNI / ref / portable_opt) and both XOFs produce identical
KAT for a given parameter.

## Note on optimizations

This implementation includes several original constant-time, table-free SIMD
optimizations for GF-arithmetic in the SNOVA signing path. A full technical
description and benchmarks will be published in a forthcoming IACR ePrint
(expected late 2026).

## License

Third-party components keep their upstream licenses:
`src/third_party/` (XKCP/Keccak), `src/nistkat/` (NIST rng/KAT harness).
