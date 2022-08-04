#!/usr/bin/env python3
import argparse
import os

import numpy as np


parser = argparse.ArgumentParser(
    description="Apply a shift to a redshift distribution")
parser.add_argument(
    "-i", "--input", nargs="*",
    help="input files with redshift distributions")
parser.add_argument(
    "-o", "--output",
    help="path where files are stored with same names as input files")
parser.add_argument(
    "-s", "--shifts", nargs="*", type=float,
    help="shifts to apply to the input files (provided in same order")
parser.add_argument(
    "--method", choices=("interp", "shift"), default="interp",
    help="shifting method: 'shift' modifies the redshift sample values, "
         "'interp' resamples the histogram values (default: %(default)s)")
args = parser.parse_args()
if len(args.shifts) != len(args.input):
    raise parser.error("number of shift parameterss do not match input files")
if not os.path.exists(args.output):
    raise OSError("output folder does not exist: {:}".format(args.output))

for infile, shift in zip(args.input, args.shifts):
    print("loading n(z) file, shifting by dz={:.3f}: {:}".format(
        shift, infile))
    z, p = np.loadtxt(infile).T  # load and unpack n(z)
    outfile = os.path.join(args.output, os.path.basename(infile))
    if os.path.exists(outfile):
        raise OSError("output file exists: {:}".format(outfile))
    if args.method == "shift":
        # apply shift to the redshift sampling points
        z += shift
        # uniformely distribute all negative samples between 0 and first
        # positive z value, interpolate probability accordingly
        nbad = np.count_nonzero(z <= 0.0)
        if nbad > 0:
            interp = np.arange(0.0, 1.0, 1.0/nbad)
            z[:nbad] = interp * z[nbad]
            p[:nbad] = interp * p[nbad]
    else:
        # interpolate the probability on z values shifted in the opposite
        # direction, preserves the original redshift sampling
        p = np.interp(z - shift, z, p, left=0.0, right=0.0)
    outdata = np.column_stack([z, p])
    # write to output path
    print("writing n(z) file: {:}".format(outfile))
    np.savetxt(outfile, outdata)
