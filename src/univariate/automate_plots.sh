#!/usr/bin/env bash

find . -iregex ".*sim_fisher.*[0-9]" -print0 | xargs -0 -P4 -I% ./plot_output.r %

rm -rf Rplots.pdf
