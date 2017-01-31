##             ##
# KEM Prototype #
#     Sage      #
##             ##

import numpy as np
from numpy.polynomial import polynomial as p


class KEM:


    def __init__(self,n=512,q=12289):
        self.q  = q
        self.n  = n

        #Assuming n = 512 / q = 12289
        self.l = 0      #Df,De - Distribution Mean
        self.o = 2.991  #Df,De - Distribution Variance

        self.R  = PolynomialRing(GF(self.q),'x')
        self.R  = self.R.quotient(x^self.n+1)

        self.R2 = PolynomialRing(GF(2),'x')
        self.R2 = self.R2.quotient(x^self.n+1)

        self.errors = 0



    # Reduction a modulo q that maps coeffs: [-(q-1)/2 , -(q-1)/2]
    def __modCoeffs(self,f,q=None):
        if (q is None): q = self.q
        q2 = q/2
        for i in range(len(f)):
            f[i] = f[i] % q
            if f[i] > q2:
                f[i] -= q
        return f;

    def __modPoly(self,f):
        return self.__modCoeffs(np.asarray(self.R(list(f)).list(),dtype=int))

    
    # Check if a given polynomial p has inverses in Rx = { Zq[X]/<x^n+1> , Z2[X]/<x^n+1> }
    def __HasInverse(self,p,Rx):
        try:
            a = Rx(list(p))
            a = a^-1
            print "ACCEPTED! in",Rx
            return True
        except Exception, e:
            print "<No inverse found>: ",len(p),e
            self.errors += 1
            if self.errors % 50==0:
                if(raw_input("Stop?")=="y"):
                    exit; #Crashes
            return False

    # Generate a polynomial
    def __genPolynomial(self,n=None,q=None,check_inverse=True):
        if (n is None): n = self.n
        if (q is None): q = self.q

        poly = [0]
        while (poly[0] == 0 or not self.__HasInverse(poly,self.R)):
            poly = np.asarray(np.random.normal(self.l,self.o,size=(n)),dtype=int)
            if (check_inverse):
                if (not self.__HasInverse(poly,self.R2)):
                    poly[0] = 0
        #return self.__modCoeffs(np.asarray(self.R(list(poly)).list(),dtype=int))
        return self.__modPoly(poly)
    

    def Keygen(self,l):
        # l has to come into play to provide strong security.
        # Not included yet.

        print "creating f"
        f = self.__genPolynomial(check_inverse=True)
        print "creating g"
        g = self.__genPolynomial()

        print "creating h"
        h = np.asarray(p.polydiv(f,g)[0], dtype=int)
        h = self.__modPoly(h)

        #h = self.R(list(h)).list()
        #h = self.__modCoeffs(np.asarray(h,dtype=int))

        print "---Keygen:"
        print "(..)  f: ",f,len(f)
        print "(OUT) g: ",g,len(g)
        print "(OUT) h: ",p.polydiv(f,g),h,len(h)
        print "\n"

        return (g,h)

    def Encapsulate(self,h):

        r = self.__genPolynomial()
        e = self.__genPolynomial()

        t = p.polymul(h,r)
        t = self.__modPoly(t)
        t = p.polyadd(2*t,e)
        t = self.__modPoly(t)


        print "---ENC:"
        print "(IN) h:",h,len(h)
        print "(..) r: ",r,len(r)
        print "(..) e: ",e,len(e)
        print "(..) t: ",t,len(t)
        print "(..) e%2: ",e%2
        print "\n"

        return (t,e%2)

    def Decapsulate(self,g,t):

        k = p.polymul(g,t)
        print "Mult g*t=",k,len(k)
        k = self.__modPoly(t)

        print "gt/g",p.polydiv((k),g)
        k = np.array(p.polydiv((k)%2,g)[0]%2,dtype=int)
        k = self.__modPoly(t)%2
        print "---DEC:"
        print "(IN)  g: ",g,len(g)
        print "(IN)  t: ",t,len(t)
        print "(OUT) k: ",k,len(k)
        print "\n"

        return k

if __name__ == "__main__":
    k  = 0
    k2 = -1
    i  = 0

    while( not np.array_equal(k,k2) and (i == 0)):
        i = i+1

        print "(",i,") ---------------------------------------------------------------------------------------V"

        kem = KEM(); #Default params {n=1024 ; q=12289}

        (Kd,Ke) = kem.Keygen(None) #Fix usage of l
        
        (c, k) = kem.Encapsulate(Ke)
        
        (k2) = kem.Decapsulate(Kd,c)

        print "\nResults:"
        print "k  :",k,len(k)
        print "k' :",k2,len(k2)
        print "Equal: ",np.array_equal(k,k2)
        print "\n"