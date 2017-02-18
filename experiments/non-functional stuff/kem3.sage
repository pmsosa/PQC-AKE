##             ##
# KEM Prototype #
#     Sage      #
##             ##

import numpy as np
from numpy.polynomial import polynomial as p
from IPython import embed


#Modular Half Reduction and g^-1
np.random.seed(20)

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
        return np.asarray(self.R(list(f)).list(),dtype=int)

    
    # Check if a given polynomial p has inverses in Rx = { Zq[X]/<x^n+1> , Z2[X]/<x^n+1> }
    def __HasInverse(self,p,Rx):
        try:
            a = Rx(list(p))
            a = a^-1
            #print "ACCEPTED! in",Rx
            return True
        except Exception, e:
            #print "<No inverse found>: ",len(p),e
            self.errors += 1
            if self.errors % 50==0:
                if(raw_input("Stop?")=="y"):
                    exit; #Crashes
            return False


    #Create an Inverse
    def __Invert(self,p):
        a = self.R(list(p))
        a = a^-1
        a = self.__modCoeffs(self.__modPoly(np.asarray(a.list())))
        return np.asarray(a)
        #np.assarray(((R(list(a)))^-1).list(),dtype=int)

    # Generate a polynomial
    def __genPolynomial(self,n=None,q=None,check_inverse=False):
        if (n is None): n = self.n
        if (q is None): q = self.q

        poly = [0]
        while (poly[len(poly)-1] == 0):
            poly = np.asarray(np.random.normal(self.l,self.o,size=n-1),dtype=int)
            
            if (check_inverse):
                if (not self.__HasInverse(poly,self.R2) or not self.__HasInverse(poly,self.R)):
                    poly[len(poly)-1] = 0
        #return self.__modCoeffs(np.asarray(self.R(list(poly)).list(),dtype=int))
        return self.__modCoeffs(self.__modPoly(poly))
    

    def Keygen(self,l):
        # l has to come into play to provide strong security.
        # Not included yet.
        h = [0]*self.n

        while (np.array_equal(h,[0]*self.n)):
            


            f = self.__genPolynomial(check_inverse=True)
            
            g = self.__genPolynomial(check_inverse=True)

            #h = np.asarray(p.polydiv(f,g)[0], dtype=int)
            h1 = np.asarray(p.polymul(f,self.__Invert(g)),dtype=int)

            h2 = self.__modPoly(h1)
            h = self.__modCoeffs(h2)


        #h = self.R(list(h)).list()
        #h = self.__modCoeffs(np.asarray(h,dtype=int))

        print "---Keygen:\n"
        print "(..)  f: ",f,len(f),"\t",self.R(list(h))
        print "(OUT) g: ",g,len(g),"\t",self.R(list(g))
        print "(OUT) h (raw): ",h1,len(h1)
        print "(OUT) h (polymod): ",h2,len(h2)
        print "(OUT) h (coefmod): ",h,len(h)
        print "\n"

        return (g,h)

    def Encapsulate(self,h):


        r = self.__genPolynomial()
        e = self.__genPolynomial()

        t1 = p.polymul(h,r)

        t2 = self.__modPoly(t1)
        t3 = self.__modCoeffs(t2)
        #print "h*r:",t
        t4 = p.polyadd(2*t3,e)
        #print "2h*r+e:",t
        t5 = self.__modPoly(t4)
        t = self.__modCoeffs(t5)


        print "---Encapsulation:\n"
        print "(IN) h: ",h,len(h)
        print "(..) r: ",r,len(r),"\t",self.R(list(h))
        print "(..) e: ",e,len(e),"\t",self.R(list(e))

        print "(..) h*r (raw): ",t1,len(t1)
        print "(..) h*r (polymod): ",t2,len(t2)
        print "(..) h*r (coefmod): ",t3,len(t3)
        print "(..) 2h*r+e: ",t4,len(t4)
        print "(OUT) 2h*r+e (polymod): ",t,len(t)
        print ""

        
        print "(OUT) e%2: ",e%2
        print "(MISC) t%2: " ,t%2
        print "\n"

        return (t,e%2)

    def Decapsulate(self,g,t):

        
        k1 = np.array(p.polymul(g,t),dtype=int)
        k2 = self.__modPoly(k1)
        k3 = self.__modCoeffs(k2)
        k4 = k3%2

        
        k5 = np.array(p.polymul(k4,self.__Invert(g)),dtype=int)
        k6 = self.__modPoly(k5)
        k7 = self.__modCoeffs(k6)
        k = k7%2
        

        print "---Decapsulation:\n"
        print "(IN)  g: ",g,len(g)
        print "(IN)  t: ",t,len(t)

        print "(..) g*t (raw):",k1,len(k1)
        print "(..) g*t (modPoly):",k2,len(k2)
        print "(..) g*t (modCoef):",k3,len(k3)
        print "(..) g*t mod 2:",k4,len(k4)
        print "(..) (g*t%2)/g (raw):",k5,len(k5)
        print "(..) (g*t%2)/g (modPoly):",k6,len(k6)
        print "(..) (g*t%2)/g (modCoef):",k7,len(k7)
        print "(OUT) k = (g*t%2)/g mod 2:", k,len(k)
        print "\n"

        return k

if __name__ == "__main__":
    k  = 0
    k2 = -1
    i  = 0

    while( not np.array_equal(k,k2) and (i == 0)):
        i = i+1

        print "(",i,") -----------------------------------------------------------------------------------V\n"

        kem = KEM(n=4); #Default params {n=1024 ; q=12289}

        (Kd,Ke) = kem.Keygen(None) #Fix usage of l
        
        (c, k) = kem.Encapsulate(Ke)
        
        (k2) = kem.Decapsulate(Kd,c)

        print "\nResults:"
        print "k  :",k,len(k)
        print "k' :",k2,len(k2)
        print "c  :",c%2
        print "Equal: ",np.array_equal(k,k2)
        print "Equal: ",np.array_equal(k,c%2)
        print "\n"
        embed()