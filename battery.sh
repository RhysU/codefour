#!/bin/bash
# Run a battery of test cases in parallel
for code in weno32 weno54
do
    echo $code spawning...
    for resolution in 8 16 32 64 128 256
    do
        resolution_name=`printf "%03d" $resolution`
        ./${code}.x "${code}.${resolution_name}.h5" $resolution &
    done
done
wait
