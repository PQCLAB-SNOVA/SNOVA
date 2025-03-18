# KAT files for SNOVA

The file `KAT-SHAKE256` in this folder contains the official 64 byte SHAKE256 digests of the KAT response (.rsp) files. The input request file `PQCsignKAT_SNOVA.req` is the same for all parameter sets. The KAT files can also be downloaded at https://github.com/PQCLAB-SNOVA/SNOVA_KAT

To generate a digest use e.g.
```
make clean PQCgenKAT SNOVA_V=37 SNOVA_O=17 SNOVA_L=2
openssl dgst -shake256 -xoflen=64 PQCsignKAT_SNOVA_37_17_2_ESK.rsp
```

This will output
```
SHAKE-256(PQCsignKAT_SNOVA_37_17_2_ESK.rsp)= cd5f1711dbc6a780effc5e0a2b97af92880c643be1640188eded8cd52c3f3efbdbdf417f43f006ddf1d5b0dc690cd21ce58e18daac30547a1494e78f8f25287c
```
