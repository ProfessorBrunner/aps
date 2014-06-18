#!/usr/bin/python
__author__='Matias Carrasco Kind'
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import healpy as hp
import os
import argparse


##    Get Arguments
###############################################################################

parser = argparse.ArgumentParser(
    description='Generate fits and bin files to test APS code.')
parser.add_argument("catalog", help="Path to catalog.dat.")
parser.add_argument("data", help="Path to output directory.")
parser.add_argument("nside", help="Pixels per edge (Power of two).", type=int)

parser.add_argument("-p", "--plot", dest="plot", action="store_true",
        help="Display plots of data.")

args = parser.parse_args()

if args.data[-1] != '/':
    args.data = args.data + '/'

is_power_of_two = lambda num: num > 0 and not (num & (num - 1))
assert is_power_of_two(args.nside)

Nside=args.nside
npix=12*Nside**2

if not os.path.isdir(args.data):
    os.system('mkdir  -p ' + args.data)

x,y,z,j,j,j=np.loadtxt(args.catalog, delimiter=',', unpack=True)
x=x-62.5/2
y=y-62.5/2
z=z-62.5/2

v=np.array([x,y,z]).T
t,p=hp.vec2ang(v)

temp=np.random.rand(len(x))
tr=np.arccos(temp*2.-1)
pr=np.random.rand(len(x))*np.pi*2.

ix=hp.ang2pix(Nside,t,p,nest=True)
ixr=hp.ang2pix(Nside,tr,pr,nest=True)

M=np.zeros(npix)
Mr=np.zeros(npix)

for i in xrange(len(ix)):
    M[ix[i]]+=1.
    Mr[ixr[i]]+=1.

ave=1.*len(x)/npix
aver=1.*len(x)/npix

M=M/ave-1.
Mr=Mr/aver-1.

fig1=plt.figure(1,dpi=100)
hp.mollview(M,fig=1,nest=True,title='Density')
fig2=plt.figure(2,dpi=100)
hp.mollview(M*0-1,fig=2,nest=True,cbar=False,cmap=cm.binary,min=-1,max=1,title='Galaxies')
hp.projscatter(t,p,marker='.',facecolor='white',s=0.1,alpha=0.9)

fig3=plt.figure(3,dpi=100)
hp.mollview(Mr,fig=3,nest=True,title='Density')
fig4=plt.figure(4,dpi=100)
hp.mollview(Mr*0-1,fig=4,nest=True,cbar=False,cmap=cm.binary,min=-1,max=1,title='Galaxies')
hp.projscatter(tr,pr,marker='.',facecolor='white',s=0.1,alpha=0.9)


fig5=plt.figure(5)
ax=fig5.add_subplot(111)
clr=hp.anafast(Mr)
cl=hp.anafast(M)

ell=np.arange(len(cl))+1

ell2=(ell[::3]+ell[1::3]+ell[2::3])/3.
cl2=(cl[::3]+cl[1::3]+cl[2::3])/3.
clr2=(clr[::3]+clr[1::3]+clr[2::3])/3.

plt.plot(ell2,cl2,'bo',label='LCDM')
plt.plot(ell2,clr2,'ro',label='random')

plt.yscale('log')
#plt.xscale('log')

plt.ylim(10**(-4.5),1e-1)
plt.xlabel(r'$\ell$',fontsize=18)
plt.ylabel(r'$C_\ell$',fontsize=18)
plt.legend(loc=0,numpoints=1)


hp.write_map(args.data+str(Nside)+'_'+str(len(x))+'_'+'lcdm.fits',M,nest=True,coord='C')
hp.write_map(args.data+str(Nside)+'_'+str(len(x))+'_'+'random.fits',Mr,nest=True,coord='C')

ell3=ell2[1:]
cl3=cl2[1:]
clr3=clr2[1:]
cl_err=np.ones(len(ell3))*0.00002
index=np.arange(len(ell3))
ell_min=ell3-1
ell_max=ell3+1

np.savetxt(args.data+'CL_'+str(Nside)+'_lcdm.dat',zip(index,ell3,ell_min,ell_max,cl3,cl_err,cl3),fmt='%d %d %d %d %.6f %.6f %.6f')
np.savetxt(args.data+'CL_'+str(Nside)+'_random.dat',zip(index,ell3,ell_min,ell_max,cl3,cl_err,cl3),fmt='%d %d %d %d %.6f %.6f %.6f')

if args.plot:
    plt.show()
      
