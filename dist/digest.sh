#!/bin/bash

sortedargs=`printf '%s\n' "$@" | sort`

for i in $sortedargs
do
openssl dgst -shake256 -xoflen=24 $i | awk -F 'PQCsignKAT_' '{print $2}' | awk '{ sub(/)=/, " "); print }'
done
