#!/usr/bin/python
"""
Generates healpix maps and band power files for test input into the
angular power spectrum code.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import healpy as hp
import pyfits
import os
import argparse

def load_catalog(catalog_path):
    """
    Load catalog of galaxies as comma seperated x, y, z
    and convert to angular coordinates.
    """
    x, y, z, _, _, _ = np.loadtxt(catalog_path, delimiter=',', unpack=True)
    x = x - 62.5/2
    y = y - 62.5/2
    z = z - 62.5/2
    v = np.array([x, y, z]).T
    t, p = hp.vec2ang(v)
    return t, p

def plot_galaxies_and_density(t, p, M):
    """Plots galaxies and pixels from polar form."""
    plt.figure(1, dpi=100)
    hp.mollview(M, fig=1, nest=True, title='Density')
    plt.figure(2, dpi=100)
    hp.mollview(M*0-1, fig=2, nest=True, cbar=False, cmap=cm.binary,
        min=-1, max=1, title='Galaxies')
    hp.projscatter(t, p, marker='.', facecolor='white', s=0.1, alpha=0.9)

def plot_anafast_band_powers(ell, cl):
    """Plot the Band powers."""
    plt.figure(5).add_subplot(111)
    plt.plot(ell, cl, 'bo', label='Band Powers')
    plt.yscale('log')
    plt.ylim(10**(-4.5), 1e-1)
    plt.xlabel(r'$\ell$', fontsize=18)
    plt.ylabel(r'$C_\ell$', fontsize=18)
    plt.legend(loc=0, numpoints=1)


def random_galaxies(total_galaxies):
    """Generate galaxies distributed with numpy random."""
    temp = np.random.rand(total_galaxies)
    t = np.arccos(temp*2. - 1)
    p = np.random.rand(total_galaxies)*np.pi*2.
    return t, p


def generate_map(t, p, nside, plot=False):
    """Generate map and bands for the given galaxy locations"""
    npix = 12*nside**2
    ix = hp.ang2pix(nside, t, p, nest=True)
    M = np.zeros(npix)
    for i in xrange(len(ix)):
        M[ix[i]] += 1.
    ave = 1.*len(t)/npix
    M = M/ave-1.0 #Calculate overdensity

    cl = hp.anafast(M)
    ell = np.arange(len(cl))+1
    ell2 = (ell[::3]+ell[1::3]+ell[2::3])/3.
    cl2 = (cl[::3]+cl[1::3]+cl[2::3])/3.

    if plot:
        plot_galaxies_and_density(t, p, M)
        plot_anafast_band_powers(ell2, cl2)
        plt.show()

    ell3 = ell2[1:]
    cl3 = cl2[1:]
    cl_err = np.ones(len(ell3))*0.00002
    index = np.arange(len(ell3))
    ell_min = ell3-1
    ell_max = ell3+1

    bands = zip(index, ell3, ell_min, ell_max, cl3, cl_err, cl3)

    return M, bands

BAND_FORMAT = '%d %d %d %d %.6f %.6f %.6f'
def write_bands(filename, bands):
    """Write bands to file name"""
    np.savetxt(filename, bands, fmt=BAND_FORMAT)

def write_fits(filename, overdensity_map, header):
    """
    Write generated map to a fits file using supplied key:val in header.
    """
    hp.write_map(filename, overdensity_map, nest=True, coord='C')
    hdulist = pyfits.open(filename, mode='update')
    prihdr = hdulist[1].header
    for key, value in header.iteritems():
        prihdr.set(key, *value)
    hdulist.close()

def main():
    """main method for using as a script"""
    #Command line parser
    parser = argparse.ArgumentParser(
        description='Generate fits and bin files to test APS code.')
    parser.add_argument("data_directory", help="Path to output directory.")
    parser.add_argument("nside", help="Pixels per edge (Power of two).",
        type=int)
    parser.add_argument("-p", "--plot", dest="plot", action="store_true",
        help="display plots of data")
    parser.add_argument("-r", "--random", type=int, dest="n_random",
        help="generate random dataset with N_RANDOM galaxies", nargs=1)
    parser.add_argument("-c", "--catalog", nargs=1, dest="catalog_file",
        help="path to catalog.dat")
    args = parser.parse_args()

    #Check Arguments
    if args.data_directory[-1] == '/':
        args.data_directory = args.data_directory[:-1]

    is_power_of_two = lambda num: num > 0 and not num & (num - 1)
    assert is_power_of_two(args.nside)

    if not os.path.isdir(args.data_directory):
        os.system('mkdir  -p ' + args.data_directory)
    
    nside = args.nside
    data_dir = args.data_directory

    #Generate output files
    if args.catalog_file:
        t, p = load_catalog(args.catalog_file[0])
        M, bands = generate_map(t, p, nside, plot=args.plot)
        #Write bands
        write_bands("{}/CL_{}_lcdm.dat".format(data_dir, nside),
            bands)
        #Write overdensity map
        fits_keys = {"NGALAXY":(len(t), "Total number of galaxies")}
        write_fits("{}/{}_{}_lcdm.fits".format(data_dir, nside, len(t)),
            M, fits_keys)

    if args.n_random:
        t, p = random_galaxies(args.n_random)
        #Write bands
        M, bands = generate_map(t, p, nside, plot=args.plot)
        write_bands("{}/CL_{}_random.dat".format(data_dir, nside),
            bands)
        #Write overdensity map
        fits_keys = {"NGALAXY":(len(t), "Total number of galaxies")}
        write_fits("{}/{}_{}_random.fits".format(data_dir, nside, len(t)),
            M, fits_keys)

if __name__ == "__main__":
    main()
