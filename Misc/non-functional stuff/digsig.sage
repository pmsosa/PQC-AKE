####      Digital Signature        ###
### Hash-and-Sign Message Recovery ###

from IPython import embed
import numpy as np
import hashlib
from sage.stats.distributions.discrete_gaussian_lattice import DiscreteGaussianDistributionLatticeSampler as GaussDistSampler

q = 12889
n = 4 # < n^n+1 > ; n = 2^i

R.<x> = ZZ['x'];
#R = PolynomialRing(ZZ,"x")
Rz = R.quotient(x^n+1)

_Rq = PolynomialRing(GF(q),"x")
Rq = _Rq.quotient(x^n+1)

R2 = PolynomialRing(GF(2),"x")
R2 = R2.quotient(x^n+1)




################################################################
# Helper Functions
################################################################

# Generate a polynomial
def __genPolynomial(n,q,check_inverse=False):

    poly = [0]
    while (poly[len(poly)-1] == 0): #Force polynomial to be filled
        poly = np.asarray(np.random.normal(0,2.99,size=n),dtype=int)
        
        if (check_inverse):
            try:
                Rq(list(poly))**-1
                R2(list(poly))**-1
            except:
                poly = [0] #Force while loop to continue;
    #return self.__modCoeffs(np.asarray(self.R(list(poly)).list(),dtype=int))
    return list(poly)

################################################################
# Actual KEM
################################################################

######################## Signature Key Generation ##############

f = __genPolynomial(n,q,check_inverse=False);
g = __genPolynomial(n,q,check_inverse=True);


h = Rq(f)*(Rq(g)**-1)
h = h.list()

Ks = (f,g)
Kv = h


######################## Sign ###################################

m = [1,2,3,4]

sigma = float(1.77*sqrt(q))

t = int(hashlib.sha256(str(m)).hexdigest(),16)

t = [1,2,3,4]

s1 = vector(ZZ,t) - GaussDistSampler(ZZ^n,sigma,t)()
s2 = vector(ZZ,[0]*n) - GaussDistSampler(ZZ^n,sigma)()

embed()



