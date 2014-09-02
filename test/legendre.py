#!/usr/bin/python
from math import radians, cos
from numpy.polynomial import legendre as leg
from tabulate import tabulate

def recuranceLegendre(x, n):
	"""
	Calculate Legendre from recurance relation
	x  cos of angle
	n  which Legendre polynomial
	"""
	hold1, hold2 = x, 1.0
	if n == 0:
		return hold2
	elif n == 1:
		return hold1
	else:
		for ii in xrange(2, n+1):
			new = ((2.0*ii-1.0)*x*hold1 - (ii-1.0)*hold2)/(1.0*ii)
			hold2 = hold1
			hold1 = new
		return new

def legendre(x,n):
	return leg.legval(x, [int(ii == n) for ii in xrange(n+2)])

print recuranceLegendre(1, 8)

samples = [1.0, 0.991741, 0.967171, 0.962783, 0.967171, 0.926852, 0.94062, 0.900524, 0.871662, 0.867519]
norm_table = []
coef_table = []
for ell in xrange(8):
	if not ell:
		continue
	row = []
	coef = (2.0*ell + 1) / (2.0*ell * (ell+1))
	print "band: {} coef: {}".format(ell, coef)
	for sample in samples:
		row.append(legendre(sample, ell)+.00001) # this is so the table works better
	norm_table.append([ell] + row)
	coef_table.append([ell] + [coef * elem for elem in row])
print "Normal"
print tabulate(norm_table)
print "After Coefficient"
print tabulate(coef_table)

band = [coef_table[3],coef_table[4],coef_table[5]]

band_result = [0 for x in range(len(samples))]
for x in (3,4,5):
	for i in range(len(samples)):
		band_result[i] += coef_table[x][i+1]

band_result = [0] + band_result
temp = [[0]+[0 for s in samples]]
temp.append(band_result)
print tabulate(temp)
#print temp

# f_string = "cos({:3}) = {:6}  recuranceLegendre({:3}, {:2}) = {:6}  legendre(cos({:3}), {:2}) = {:6}"
# for ii in xrange(0, 10):
# 	#print "L_{}".format(ii)
# 	for jj in xrange(0,180, 10):
# 		print f_string.format(
# 			            jj, round(cos(radians(jj)),3),
# 			            jj, ii, round(recuranceLegendre(cos(radians(jj)), ii), 3), 
# 			            jj, ii, round(legendre(cos(radians(jj)), ii), 3))

