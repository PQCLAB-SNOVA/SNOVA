# KAT files for SNOVA

The file `KAT-SHAKE256` in this folder contains 24 byte SHAKE256 digests of the KAT response (.rsp) files. The input request file `PQCsignKAT_SNOVA.req` is the same for all parameter sets.

To generate a digest use e.g.
```
make clean PQCgenKAT SNOVA_V=37 SNOVA_O=17 SNOVA_L=2
openssl dgst -shake256 -xoflen=24 PQCsignKAT_SNOVA_37_17_2_ESK.rsp
```

This will output
```
SHAKE-256(PQCsignKAT_SNOVA_37_17_2_ESK.rsp)= cd5f1711dbc6a780effc5e0a2b97af92880c643be1640188
```

The KAT files can also be downloaded at https://github.com/PQCLAB-SNOVA/SNOVA_KAT
