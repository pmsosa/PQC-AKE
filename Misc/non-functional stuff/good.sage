# Lightweight KEM Example #

from IPython import embed
import numpy as np
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

def modCoeffs(f,pp):
    clist=f.list()
    p2=int(pp/_sage_const_2 )
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

g_inv = Rq(g)**-1
print "g_inv:",g_inv

h = R(f)*g_inv
print "h:",h


######################## Encapsulation #########################
print "\nEncapsulation"

#Randomly Generating r and e (with small norms)
r = [ 1,-2, 0, 1] ## __genPolynomial(n,q,check_inverse=False)
e = [ 2,-1,-1, 1] ##__genPolynomial(n,q,check_inverse=False)
print "r:",r
print "e:",e

t = 2*h*R(r)+R(e)
print "t:",t
k = R2(e)

######################## Decapsulation ##########################

k2 = R(g)*t
k2 = R2(k2.list())
k2 = k2*(R2(g)**-1)
k2 = R2(k2)


print "\nResults:----------\n"
print "k: ", k.list()
print "k2: ", k2.list()
print "t: ", t.list()
print "t%2: ",R2(t.list()).list()


print "\n\nTest"
print "2fr+ef in Z  :",(2*(R(f)*R(r))%(x**4+1)+(R(e)*R(g))%(x**4+1))
print "2fr+eg in Zq :",Rq(2*R(f)*R(r)+R(e)*R(g))
print "gt mod q     :",Rq(R(g)*R(t.list()))
print "\n\n"
embed()