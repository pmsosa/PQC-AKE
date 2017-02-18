import numpy as np
from numpy.polynomial import polynomial as p


n = 4 		# Z/<x^4+1>
q = 12289	# Sufficiently big such that you don't have to do %q ever.

print "\n-Params-"
print "n:",n
print "q:",q


#KEM Key Generation
print "\n-Keygen-"

l = 2.0

f = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))
while (f[0] == 0):
	f = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))
#np.asarray(np.random.normal(0,3.0,size=(n)),dtype=int)

print "f:",f

g = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))
while (g[0] == 0):
	g = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))

print "g:",g

h = p.polydiv(f,g)

print "h=(f/g):",h

h=np.asarray(h[0],dtype=int)

#KEM Encapsulation

print "\n-Key Encapsulation-"

r = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))
while (r[0] == 0):
	r = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))

print "r:",r

e = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))
while (e[0] == 0):
	e = np.asarray(np.asarray(np.random.normal(0,l,size=(n)),dtype=int))

print "e:",e

t = p.polymul(h,r)

print "[c] t=(2hr)+e:",t,"+",e,"=",p.polyadd(h,e)

t = p.polyadd(h,e)

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
print "k':",np.asarray(p.polydiv((gt%2),g)[0]%2,dtype=int)
print "k==k':",np.array_equal(e%2,np.asarray(p.polydiv((gt%2),g)[0]%2,dtype=int))