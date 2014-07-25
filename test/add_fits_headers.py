#!/usr/bin/python
"""
Add fits headers that are required by the code.
"""
import argparse
from generate_inputs import add_fits_headers



def main():
    """main method for using as a script"""
    #Command line parser
    parser = argparse.ArgumentParser(
        description='Generate fits and bin files to test APS code.')
    parser.add_argument("fits_file", help="Path fits file.")
    parser.add_argument("-nside", type=int, dest="nside",
        help="Value for fits argument NSIDE", nargs=1)
    parser.add_argument("-ngalaxy", type=int, dest="ngalaxy",
        help="Value for fits argument NGALAXY.", nargs=1)
    args = parser.parse_args()

    header_dict = {}
    header_dict['NGALAXY'] = (args.ngalaxy[0], "Total number of galaxies")
    header_dict['NSIDE'] = (args.nside[0], "Healpix subdivision")

    add_fits_headers(args.fits_file, header_dict)

if __name__ == "__main__":
    main()