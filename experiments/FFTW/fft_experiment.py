from IPython import embed
import numpy as np
from numpy.polynomial import polynomial as p
import timeit

#############################################
def gen_poly(n,q):
    global hlpr
    l = 0 #Gamma Distribution Location (Mean "center" of dist.)
    poly = np.floor(np.random.normal(l,size=(n)))
    while (len(poly) != n):
        poly = np.floor(np.random.normal(l,size=(n)))
        poly = np.floor(p.polydiv(poly,hlpr)[1]%q)
    return poly
#############################################


#############################################
def mult(a,b,type=1):
    kex_start = timeit.default_timer()
    if type == 1:
        c = p.polymul(a,b)
    else:
        c = np.fft.irfft(np.fft.rfft(a,2*n)*np.fft.rfft(b,2*n),2*n)
    kex_end = timeit.default_timer()
    return c,kex_end-kex_start 

#############################################


n = 1024        #n and q, based on D. Stebila's talk
q = 2**32-1
hlpr = [1] + [0] * (n-1) + [1]


A = np.floor(np.random.random(size=(n))*q)%q
A = np.floor(p.polydiv(A,hlpr)[1])

s = gen_poly(n,q)


(r1,t1) = mult(A,s)
(r2,t2) = mult(A,s,2)

print r1
print r2

embed()