#!/bin/bash
# Run a battery of test cases in parallel
for code in weno32 weno54
do
    echo $code running...
    for resolution in 25 50 100 125 150 175 200 250 300 350
    do
        ./${code}.x "${code}.${resolution}.h5" $resolution &
    done
    wait
done