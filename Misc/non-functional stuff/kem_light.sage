# Lightweight KEM Example #

from IPython import embed
#import numpy as np
#import sys

q = 12889
n = 4 # < n^n+1 > ; n = 2^i

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

# Maps coefficients to [-(p-1/2),(p-1)/2]
# In the case pp=2 notice that it will just perform regular modulo 2
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

print "\nNote: [5,6,7,8] refers to 8x^3 + 7x^2 + 6x + 5"
######################## Key Generation ########################
print "\nKeyGen"

#Randomly choosing f and g (wiht small norms)
f = [-3,-1, 0, 1] # 1x^3 + 0x^2 - 1x - 3    ## __genPolynomial(n,q,check_inverse=True)
g = [ 1,-2, 1, 1] # 1x^3 + 1x^2 - 2x + 1    ## __genPolynomial(n,q,check_inverse=True)

print "f:",f
print "g:",g

# Checking both f,g are invertible in Zp/<x^n+1> and Z2/<x^n+1>
# Otherwise the following two lines would crash
Rq(f)**-1 ; R2(f)**-1
Rq(g)**-1 ; R2(g)**-1
print "g & f are invertible"



# G inverse
#g_inv = (inv_poly_mod_prime_pow__(R(g))[1]).list()
#print Rq(Rz(g)*Rz(g_inv))
g_inv = (Rq(g)**-1)
g_inv = modCoeffs(g_inv,q)
print "g^-1:", g_inv.list()

print "(g*g^-1):",Rq(Rz(g)*Rz(g_inv))

#print "Inverse_Sam:",inv_poly_mod_prime_pow__(R(g))[1] # Doesn't seem to work
#print "Inverse_Rq:",Rq(g)**-1 #What im using right now
#print "Inverse_R2:",R2(g)**-1 #Doesn't seem to work + breaks the encapsulation

h = R(f)*(R(g_inv))
h = modCoeffs(Rz(h),q)
#h = Rq(h)
h = h.list()  # (Doing this to move things around as lists)
print "h:",h


######################## Encapsulation #########################
print "\nEncapsulation"

#Randomly Generating r and e (with small norms)
r = [ 1,-2, 0, 1] ## __genPolynomial(n,q,check_inverse=False)
e = [ 2,-1,-1, 1] ##__genPolynomial(n,q,check_inverse=False)
print "r:",r
print "e:",e

# Notice, unlike f & g, we don't care if r & e are invertible.

# t = 2hr + e (mod q)
t = 2*R(h)*R(r) + R(e)  # t = 2hr + e
t = Rz(t)               # t modulo the polynomial <x^n+1>
t = modCoeffs(t,q)      # coeffs of t modulo q
#t = Rq(t)
t = t.list()           
print "t (c):",t

k = R2(e).list()        # k = e % 2
print "k:",k

######################## Decapsulation ##########################
print "\nDecapsulation"

# k = ( (gt mod q mod 2)/ g )(mod 2)
k2 = R(g)*R(t)          # k = gt
k2 = Rz(k2)             # gt modulo the polynomial <x^n+1>
k2 = modCoeffs(k2,q)    # k = (gt) mod p
#k2 = Rq(k2).list()
k2 = R2(k2).list()      # k = (gt mod p) mod 2

#print "Nominator (gt mod q mod 2): ",k2 

k2 = R(k2)*R(g_inv)     # (gt mod p mod 2 )/g
k2 = Rz(k2)             # (gt mod p mod 2 )/g modulo the polynomial <x^n+1>
k2 = R2(k2)             # ((gt mod p mod 2 )/g) mod 2
k2 = k2.list()

### RESULTS ###
print "[Goal: k == k' ; t != k]"
print "k: ",k
print "k':",k2
if (k==k2): print "Valid!"
else: print "Invalid!"
print "t: ",t,"\t%2= ",modCoeffs(R(t),2).list()
embed()

##################################################################
##################################################################
##################################################################
# THIS CODE BELOW IS NOT IN USE
##################################################################
# Code Borrowed from Sam Green's NTRU Implementation
##################################################################

def posResidue(poly, base):
    clist = poly.list()
    for i in range(len(clist)):
        clist[i] = clist[i]%base
        if clist[i] < 0:
            clist[i] += base
    return R(clist)

def inv_poly_mod2__(poly):
    k=0;b=1;c=0*x;
    f=poly;g=x^n-1
    f=modCoeffs(f, 2)
    res=False
    while True:
        while f(0)==0 and not f.is_zero():
            f=f.shift(-1)
            c=c.shift(1)
            #c=self.modCoeffs(c, 2)
            c = posResidue(c,2)
            k+=1
        if f.is_one():
            e=(-k)%n
            retval= x^e*b 
            res=True
            break
        elif f.degree()==-1 or f.is_zero():
            break
        if f.degree()<g.degree():
            f,g=g,f
            b,c=c,b
        f=f+g
        b=b+c
        #f=self.modCoeffs(f, 2)
        f=posResidue(f,2)
        #c=self.modCoeffs(c, 2)
        c=posResidue(c,2)
    if res:
        retval=retval%(x^n-1)
        #retval=self.modCoeffs(retval, 2)
        retval = posResidue(retval, 2)
        return True, retval
    else:
        return False,0

def inv_poly_mod_prime_pow__(poly):
    res,b=inv_poly_mod2__(poly)
            #print "Inside __inv_poly_mod_prime_pow__(): res={}, b={}".format(res,b)
    if res:
        qr=2
        while qr<q:
            qr=qr^2
            b=b*(2-poly*b)
            b=b%(x^n-1)
            #b=self.modCoeffs(b, self.q)
            b = posResidue(b, q)
        return True,b
    else:
        return False,0


################ RANDOM RAMBLINGS ################
#a = modCoeffs(Rq(r)*modCoeffs(Rq(f)/Rq(g),q) + Rq(e),q)
#a = modCoeffs(modCoeffs(a*Rq(g),q),2)
#a = modCoeffs(a/Rq(g),2)
#a = modCoeffs( ( (Rq(r)*Rq(f)/Rq(g) + Rq(e) )*Rq(g) ), 2)
#a = Rq(f)/Rq(g)*Rq(r)+Rq(e)