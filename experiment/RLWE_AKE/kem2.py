import numpy as np
from numpy.polynomial import polynomial as p



n = 4 		# Z/<x^4+1>
q = 12289	# Sufficiently big such that you don't have to do %q ever.

print "\n-Params-"
print "n:",n
print "q:",q


#KEM Key Generation
print "\n-Keygen-"

f = np.asarray([1,1,9,2])

print "f:",f

g = np.asarray([1,0,1,1])

print "g:",g

h = p.polydiv(f,g)

print "h=(f/g):",h

h=h[0]

#KEM Encapsulation

print "\n-Key Encapsulation-"

r = np.asarray([1,0,1,4])

print "r:",r

e = np.asarray([1,3,-2,1])

print "e:",e

t = p.polymul(h,r)

print "[c] t=(2hr)+e:",t,"+",e,"=",h+e

t = h+e

print "[k]:",e%2

#KEM Decapsulation

print "\n-Key Decapsulation-"

print "k=( gt(%2) )/g (%2)"

gt = p.polymul(g,t)

print "gt:",gt

print "gt%2:",gt%2

print "(gt%2)/g",p.polydiv((gt%2),g)

print "((gt%2)/g)%2", p.polydiv((gt%2),g)[0]%2

print "k :",e%2
print "k':",p.polydiv((gt%2),g)[0]%2