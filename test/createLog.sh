#!/bin/sh

# creates the log file against which we test the output of SMSmain.py. Use with care!

cd .. && python ./SMSmain.py | tee test/SMSmain.log
