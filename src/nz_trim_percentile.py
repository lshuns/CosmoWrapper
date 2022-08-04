#!/usr/bin/env python3
import argparse
import os

import numpy as np


parser = argparse.ArgumentParser(
    description="Truncate a redshift distribution")
parser.add_argument(
    "infile", help="input file with redshift distribution")
parser.add_argument(
    "outfile", help="output file path")
parser.add_argument(
    "percentile", type=float,
    help="percentile (in percent) at which the distribution is truncated")
args = parser.parse_args()

print("loading n(z) file, clipping {:.1f}-th percentile: {:}".format(
    args.percentile, args.infile))
indata = np.loadtxt(args.infile)
if os.path.exists(args.outfile):
    raise OSError("output file exists: {:}".format(args.outfile))
# compute simple PDF and CDF
riemann_sum = np.diff(indata[:, 0]) * indata[:-1, 1]
norm = riemann_sum.sum()
pz = np.append(riemann_sum, 0.0) / norm
Pz = np.cumsum(pz)
# find the index at which the percentile is reached
idx_clip = np.argmin(np.abs(Pz - args.percentile / 100.0))
# clip the data and restore original normalisation
outdata = indata.copy()
outdata[idx_clip:, 1] = 0.0
riemann_sum = np.diff(outdata[:, 0]) * outdata[:-1, 1]
outdata[:, 1] *= norm / riemann_sum.sum()  # renormalise
# write to output path
print("writing n(z) file: {:}".format(args.outfile))
np.savetxt(args.outfile, outdata)
