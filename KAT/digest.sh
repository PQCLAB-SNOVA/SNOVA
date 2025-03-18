#!/bin/bash

for i in "$@" 
do
openssl dgst -shake256 -xoflen=64 $i | awk -F 'PQCsignKAT_' '{print $2}' | awk '{ sub(/.rsp)=/, " "); print }'
done
