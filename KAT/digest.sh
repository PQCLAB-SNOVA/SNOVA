#!/bin/bash

for i in "$@" 
do
openssl dgst -shake256 -xoflen=24 $i | awk -F 'PQCsignKAT_' '{print $2}' | awk '{ sub(/)=/, " "); print }'
done
