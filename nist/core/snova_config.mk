# Frozen build configuration for this implementation folder.
#
# This folder is: AVX2 only - runs on ANY AVX2 CPU (Haswell, 2013 and later)
# When GFNI=1: AVX2+GFNI, fastest - REQUIRES GFNI (Ice Lake+ / Zen 4+); SIGILLs without it
#
# The makefile in this folder -includes this file, so plain
#
#     make PQCgenKAT SNOVA_V=.. SNOVA_O=.. SNOVA_Q=.. SNOVA_L=.. SNOVA_R=.. SNOVA_M1=.. SNOVA_M2=..
#
# already builds the avx2 variant. You do not need to pass ARCH= or GFNI=,
# and the Additional_Implementations/Makefile one level up passes the same values anyway.
# Changing the lines below builds something other than what this folder is meant to hold.

ARCH=x86_avx2
GFNI=0
