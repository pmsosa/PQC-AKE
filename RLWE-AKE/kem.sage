# Lightweight KEM Example #

from IPython import embed
import numpy as np
#import sys

q = 12889
n = 512 # < n^n+1 > ; n = 2^i

R.<x> = ZZ['x'];
#R = PolynomialRing(ZZ,"x")
Rz = R.quotient(x^n+1)

Rq = PolynomialRing(GF(q),"x")
Rq = Rq.quotient(x^n+1)

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


# Reduction a modulo q that maps coeffs: [-(q-1)/2 , -(q-1)/2]
def modCoeffs(f,pp):
    clist=f.list()
    p2=int(pp/2)
    for i in range(len(clist)):
        clist[i] = int(clist[i])%pp
        #print clist[i],p2,clist[i]>p2,int(clist[i])-pp,type(clist[i]),clist[i]*-1
        if clist[i]>p2:
            clist[i]= int(clist[i]) - pp
            #print clist[i]
    return R(clist)



################################################################
# Actual KEM
################################################################
def kem(n,q):

    print "\nNote: [5,6,7,8] refers to 8x^3 + 7x^2 + 6x + 5"
    ######################## Key Generation ########################
    print "\nKeyGen"

    #Randomly choosing f and g (wiht small norms)
    f = __genPolynomial(n,q,check_inverse=True) # [-3,-1, 0, 1] # 1x^3 + 0x^2 - 1x - 3    ## __genPolynomial(n,q,check_inverse=True)
    g = __genPolynomial(n,q,check_inverse=True) # [ 1,-2, 1, 3] # 1x^3 + 1x^2 - 2x + 1    ## __genPolynomial(n,q,check_inverse=True)

    print "f:",f
    print "g:",g

    # Checking both f,g are invertible in Zp/<x^n+1> and Z2/<x^n+1>
    # Otherwise the following two lines would crash
    Rq(f)**-1 ; R2(f)**-1
    Rq(g)**-1 ; R2(g)**-1
    print "g & f are invertible"

    g_inv = Rq(g)**-1
    g_inv = modCoeffs(g_inv,q).list()

    h = R(f)*R(g_inv)
    h = h%(x**n+1)
    h = modCoeffs(h,q)
    print "h: ",h
    h = h.list()


    ######################## Encapsulation #########################
    print "\nEncapsulation"

    #Randomly Generating r and e (with small norms)
    r = __genPolynomial(n,q,check_inverse=False) # [ 1,-2, 0, 1] ## __genPolynomial(n,q,check_inverse=False)
    e = __genPolynomial(n,q,check_inverse=False) #[ 2,-1,-1, 1] ##__genPolynomial(n,q,check_inverse=False)
    print "r:",r
    print "e:",e

    t = R(h)*R(r)
    t = t%(x**n+1)
    t = modCoeffs(t,q)
    t = 2*t
    t = modCoeffs(t,q)

    print "t1:,",t

    t = t+R(e)
    t = modCoeffs(t,q)

    print "t:,",t

    t = t.list()

    k = R2(e).list()


    ######################## Decapsulation ##########################

    k2 = R(g)*R(t)
    k2 = k2%(x**n+1)
    k2 = modCoeffs(k2,q)

    print "k''':",k2

    k2 = modCoeffs(k2,2)
    print "k'':",k2

    k2 = k2*R((R2(g)**-1).list())
    k2 = k2%(x**n+1)

    print "k':",k2
    k2 = modCoeffs(k2,2)

    k2 = k2.list()


    print "\nResults:----------\n"
    print "k: ", k
    print "k2: ", k2
    print "t: ", t
    print "t%2: ",R2(t).list()



    ####
    print "\n\n\n\n Testing"
    print "2fr+ge",Rq(2*R(f)*R(r)+R(g)*R(e))
    print "gt mod q",Rq(R(g)*R(t))

    print "Valid:",k == (k2 + [0]*(len(k)-len(k2)))
    k2
    print len(k),len(k2)





kem(n,q)