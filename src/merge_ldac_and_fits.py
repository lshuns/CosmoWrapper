#  We've got two catalogues in the same order
#  One is ldac-fits, one is fits
#  We want to merge them and create one final master ldac catalogue

import ldac
import sys
import numpy as np
from astropy.io import fits

#==============================

# command line input/output files
if len(sys.argv) != 6: 
    print ("Usage: %s InputCat_Ldac InputCat_Fits OutputCat 'Blind List' 'Columns to copy'" % sys.argv[0] )
    sys.exit(1)
else:
    infile_ldac = sys.argv[1]
    infile_fits = sys.argv[2]
    outfile = sys.argv[3] 
    blindlist = sys.argv[4]
    blindlist = blindlist.split()
    print("Reading user blindlist: %s" % blindlist)
    columnlist = sys.argv[5]
    columnlist = columnlist.split()
    seen = set()
    result = []
    for item in columnlist:
        if item not in seen:
            seen.add(item)
            result.append(item)
    columnlist = result
    print("Reading user columnlist: %s" % columnlist)



print(infile_ldac)
print(infile_fits)
print(outfile)
print(blindlist)

#==============================

# Open the ldac catalogue using functions in ldac.py
ldac_cat = ldac.LDACCat(infile_ldac)
ldac_table = ldac_cat['OBJECTS']

# Open the fits catalogue using astropy
f=fits.open(infile_fits)
fitsext = 1 # assuming the info that we want is in the first extension

# Lets perform a sanity check before we go any further
SeqNr_LDAC = ldac_table['SeqNr']
SeqNr_FITS = f[fitsext].data['SeqNr']
THELI_INT_LDAC = ldac_table['THELI_INT']
THELI_INT_FITS = f[fitsext].data['THELI_INT']
#ESO_ID_LDAC = ldac_table['ID']
#ESO_ID_FITS = f[fitsext].data['ID']

#if np.array_equal(ESO_ID_LDAC, ESO_ID_FITS):
#    print ("Everything is well with the world - continue")
#else:
#    print ("ESO_IDs do not match - exit")
#    sys.exit(1)

if np.array_equal(SeqNr_LDAC, SeqNr_FITS):
    print ("Everything is well with the world - continue")
else:
    print ("SeqNrs do not match - exit")
    sys.exit(1)

if np.array_equal(THELI_INT_LDAC, THELI_INT_FITS):
    print ("Everything is well with the world - continue")
else:
    print ("THELI_INT's do not match - exit")
    sys.exit(1)

# These are the columns that we want to merge into the 
# master ldac fits table, for blinds A, B and C
print("Looping over Columns")
for col in columnlist: 
    for blind in blindlist:
        colname="%s_%s"%(col,blind)
        comment="SOM Flag for case %s and blind %s"%(col,blind)
        print ("adding column %s to ldac table" %colname)
        #add the column into the ldactable
        ldac_table[colname] = np.round(f[fitsext].data[colname]).astype("int")
        #add the comments and units
        ldac_table.set_comment(colname,comment)
        ldac_table.set_unit(colname, ' ')

# and save the new ldac catalogue automatically including FIELDS table:
print ("writing out file %s" %outfile)
ldac_cat.saveas(outfile)


