#!/usr/bin/env python3
import argparse
import os

from astropy.io import fits as pyfits


ldac_template = os.path.join(os.path.dirname(os.path.abspath(__file__)), "LDAC.cat")

parser = argparse.ArgumentParser(
    description="Create a fake LDAC table file from a single FITS table "
                "extension.")
parser.add_argument(
    "fitsfile", help="path to the FITS file")
parser.add_argument(
    "-o", "--output", required=False,
    help="output LDAC file path (default: input file with .cat extension)")
parser.add_argument(
    "-n", type=int, default=1,
    help="index of the table HDU in the input file (default: %(default)s)")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.output is None:
        args.output = os.path.splitext(args.fitsfile)[0] + ".cat"
    # load the data file
    with pyfits.open(args.fitsfile) as fits:
        data = fits[args.n]
        data.name = "OBJECTS"
        # combine with template LDAC HDU
        with pyfits.open(ldac_template) as ldac:
            ldac.append(data)
            print("==> LDAC file schema:")
            ldac.info()
            print("==> Writing file to: " + args.output)
            ldac.writeto(args.output)
