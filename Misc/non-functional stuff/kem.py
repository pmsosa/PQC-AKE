##############################################
## Key Encapsulation Mechanism (KEM) Module ##
##############################################


import numpy as np
from numpy.polynomial import polynomial as p
import random

#ISSUE: Your main issue is that you are not properly checking that everything stays in Zq[x]/<x^n+-1>
#ISSUE: len(hlpr) != poly

class KEM:

    def __init__(self,n=1024,q=12289):
        self.n = n
        self.q = q
        self.hlpr = [1] + [0] * (n-1) + [1] # x^1024 + 1

    def __gen_poly(self,n,q):

        #For now. If not 1024 just crash
        if (n == 1024):
            template = [0]*166 + [1]*150 + [2]*118 + [3]*80 + [4]*46 + [5]*22 + [6]*9 + [7]*3 + [8]*1 \
                               +[-1]*150 +[-2]*118 +[-3]*80 +[-4]*46 +[-5]*22 +[-6]*9 +[-7]*3 +[-8]*1 

            np.random.shuffle(template)
            #Avoid having an element of lesser degree.
            while(template[0] == 0):
                np.random.shuffle(template)

        else:
            return self.__gen_poly_dynamic(n,q)
            #temp_sign = [(-1)**random.randint(0,1) for i in range(1024)]

        return np.array(template)

    def polymod(self,p1):
        while(len(p1) >= len(self.hlpr)):
            print "modding",len(p1),len(self.hlpr)
            print p.polydiv(p1,self.hlpr)
            p1 = (p.polydiv(p1,self.hlpr)[1])
        return p1

    def coefmod(self,p1):
        for i in range(0,len(p1)):
            if (abs(p1[i]) > (q-1)/2):
                p1[i]=p1[i]%q


    #not being used
    def __gen_poly_dynamic(self,n,q):
        l = 0 #Gamma Distribution Location (Mean "center" of dist.)
        o = 3.467
        poly = np.floor(np.random.normal(l,size=(n)))
        while (len(poly) != n or poly[0] == 0):
            poly = (np.random.normal(l,o,size=(n))).astype(int)
            poly = (p.polydiv(poly,self.hlpr)[1]).astype(int)
        return poly


    #KeyGen(l) = (Ke,Kd)
    def Keygen(self,l):
        #TODO: l has to come into play properly to provide strong security
        #For now I'm letting the polynomials be small, uniformly distributed centered at 0
        
        f = np.asarray([-7,1,6,3])#self.__gen_poly(self.n,self.q)          # f <-$- Df
        g = np.asarray([-2,-1,-1,-2])#self.__gen_poly(self.n,self.q)          # g <-$- Df
        # We gotta check that f,g are invertible

        h = (p.polydiv(f,g)[0]).astype(int) # h = f/g (mod q)

        print "---Keygen:"
        print "(..)  f: ",f,len(f)
        print "(OUT) g: ",g,len(g)
        print "(OUT) h: ",p.polydiv(f,g),h,len(h)
        print "\n"

        return (g,h)                                # Return (Kd, Ke) = (g,h)

    #Enc(h) = (c,k)
    def ENC(self,h):
        r = np.asarray([2,-2,0,0])#self.__gen_poly(self.n,self.q)          # r <-$- De
        e = np.asarray([-1,-2,-1,-1])#self.__gen_poly(self.n,self.q)          # e <-$- De
        
        t = self.polymod(p.polymul(h,r))           # t = 2hr + e (mod q)
        t = p.polyadd(2*t,e)
        t = t.astype(int)    
        
        print "---ENC:"
        print "(IN) h:",h,len(h)
        print "(..) r: ",r,len(r)
        print "(..) e: ",e,len(e)
        print "(..) t: ",t,len(t)
        print "(..) e%2: ",e%2
        print "\n"

        return (t,e%2)                              # Return (c,k) = (t,e mod 2)

    #Dec(Ke,k)
    def DEC(self,g,t):


        # Return k = (gt (mod q) (mod 2) )/g (mod 2)
        k = p.polymul(g,t)
        print "Mult g*t=",k,len(k)
        k = self.polymod(k)
        print k
        print "gt/g",p.polydiv((k),g)
        k = (p.polydiv((k)%2,g)[0]%2).astype(int)

        
        print "---DEC:"
        print "(IN)  g: ",g,len(g)
        print "(IN)  t: ",t,len(t)
        print "(OUT) k: ",k,len(k)
        print "\n"

        return k




#Quick Test
if __name__ == "__main__":
    k = 0
    k2 = -1
    i = 0
    while (not np.array_equal(k,k2) and (i == 0)):
        i = i+1

        print "(",i,") ---------------------------------------------------------------------------------------V"

        kem = KEM(n=4); #Default params {n=1024 ; q=12289}

        (Kd,Ke) = kem.Keygen(None) #Fix usage of l
        
        (c, k) = kem.ENC(Ke)
        
        (k2) = kem.DEC(Kd,c)


        print "\nResults:"
        print "k  :",k,len(k)
        print "k' :",k2,len(k2)
        print "Equal: ",np.array_equal(k,k2)
        print "\n"