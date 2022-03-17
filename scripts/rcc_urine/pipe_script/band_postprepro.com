#!/bin/csh
# preprocessing of the NMR spectra
clean.com
proc.com
# filter fid into frequency bands and convert the format
bands.com
postproc.com
