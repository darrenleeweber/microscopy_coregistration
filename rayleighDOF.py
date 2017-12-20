#!/usr/bin/env python

from optparse import OptionParser
import sys


def main():
# Parse command line options
p = OptionParser()
p.add_option("-w", "--wavelength", dest="wavelength", default=0.5,
action="store", type="float",
help="Light wavelength (default is 0.5 micrometers)", 
metavar="WAVELENGTH")
p.add_option("-r", "--refractiveIndex", dest="refractiveIndex", default=1.0,
action="store", type="float",
help="Refractive index of immersion medium (1.0 air, 1.33 water, 1.515 oil)",
metavar="REFRACTIVE_INDEX")
p.add_option("-n", "--numericalAperture", dest="numericalAperture", 
action="store", type="float",
help="Numerical aperture of the objective lens", 
metavar="NUMERICAL_APERTURE")
(options, args) = p.parse_args()
print options

if __name__ == "__main__":
    main()
    sys.exit()

