import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pyfits

def is_power_of_two(num):
    return num > 0 and not num & (num - 1)

def add_fits_headers(filename, header_dict):
    """Add fits headers to a fits file from dictionary"""
    hdulist = pyfits.open(filename, mode='update')
    prihdr = hdulist[1].header
    for key, value in header_dict.iteritems():
        prihdr.set(key, value[0], value[1])
    hdulist.close()

class aps_input:
    def __init__(self, source_file, nside):
        l, cl = np.loadtxt(source_file, unpack=True)
        cl = 2.0*cl*np.pi / (l*(l + 1.0))
        self.cl = cl
        self.l = l
        self.nside = nside
        self.map = hp.synfast(cl, nside)
        self.make_mean_zero()
        self.get_bands()

    def patch_mask(self, mask):
        mask = hp.ud_grade(mask, nside_out=self.nside, 
                order_in='NESTED', order_out='NESTED')
        self.idx_masked = np.where(mask == 0.)[0]
        self.map[self.idx_masked] = -2.0
        self.make_mean_zero()

    def reduce_to_count(self, count):
        self.count_pixels()
        if self.pixels < count:
            print "Cannot increase count {} -> {}".format(self.pixels, count)
        idx = 0
        while self.pixels > count:
            if self.map[idx] >= -1.0:
                self.map[idx] = -2.0
                self.pixels += -1
            idx += 1


    def make_mean_zero(self):
        idx = np.where(self.map >= -1.0)[0]
        mean = np.mean(self.map[idx])
        self.map = self.map - mean
        # Synfast may always generate mean = 0 maps?
        # mean = np.mean(self.map[idx])
        # print mean

    def count_pixels(self):
        self.pixels = len(np.where(self.map >= -1.0)[0])

    def get_bands(self):
        delta_l = 4
        min_l = 7
        self.bands = int(np.floor( (self.nside*3.0 - min_l)/delta_l ))
        self.set_bands(self.bands, delta_l, min_l)

    def set_bands(self, bands, delta_l, min_l):
        self.bands = bands
        self.l_start = np.array([min_l+delta_l*i for i in xrange(self.bands)])
        self.l_end = self.l_start+delta_l-1
        self.l_center = self.l_start+delta_l/2
        er = np.ones(self.bands)*0.000001
        self.expected_cl = self.cl[self.l_center-1]+np.random.randn(self.bands)*er[0]
        self.initial_cl = np.ones(self.bands)*1e-5+np.random.randn(self.bands)*er[0]

    def write(self, fits_file, band_file, ngalaxies=1e6):
        self.count_pixels()
        print "saving %s" % band_file
        np.savetxt(band_file,
                zip(range(self.bands), self.l_center, self.l_start, self.l_end,
                          self.initial_cl, np.zeros(self.bands), self.initial_cl ),
                fmt = '%d %d %d %d %.6f %.6f %.6f')
        hp.write_map(fits_file, self.map, nest=True, coord='C')

        header_dict = {}
        header_dict['NGALAXY'] = (ngalaxies, "Total number of galaxies")
        header_dict['NSIDE'] = (self.nside, "Healpix subdivision")
        add_fits_headers(fits_file, header_dict)


    def moll_graph(self):
        M = np.copy(self.map)
        M[self.idx_masked] = hp.UNSEEN
        hp.mollview(M,  xsize = 1000,  nest = True)
        plt.show()
