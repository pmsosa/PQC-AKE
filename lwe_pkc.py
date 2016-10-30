#Learning With Errors

import numpy as np


'''
		Publicly Available
		A <- Z (n*n, q)
			n*n = size or Matrix
			  q = "max" size of coefficient (mod q) 

Alice 								Bob
  s <- X^n 						r <- X^n (vector of size n, sampled randomly from error dist.)

  	|__________u = s*A + e__________=> ( e is an error; sampled from error distribution)
  	<=_________b = A*r + x___________| ( x is an error; sampled from error distribution)

  	<=____b' = u*r+x' + (bit)*q/2____|


  b - s ~ bit*(q/2)

'''




def gen_A(n,q):
	location = 0 #Gamma Distribution Location (Mean "center" of dist.)
	#a = np.floor(np.random.normal(0,size=(n,n))%q)
	a = np.floor(np.random.random(size=(n,n))*q)%q
	return np.matrix(a);

def gen_vect(n):
	s = np.floor(np.random.normal(0,size=(n,1)))
	return np.matrix(s);

def gen_num():
	s = np.floor(np.random.normal(0))
	return np.matrix(s);




#Step 0. (UNIVERSE) Define n and q
n = 1024   #NewHope Standards
q = 12289  #NewHope Standards
# k = 1
bit = 1 # Bit that bob wants to send to alice

#Step 1. (UNIVERSE) Generate A
A = gen_A(n,q) # N*N matrix whith integer matrices that are {0, ..., q-1}
#print "A: ", A

#Step 2. (ALICE)
s = gen_vect(n)	#Alice Generate Secret <s1 .... sn>
e = gen_vect(n)	#Alice Generate Error  <e1 .....en>
#print "S: ",s,"E:",e


## St       A
## (1*5)  (5*5)
ut = np.transpose(s)*A + np.transpose(e) #Alice send to Bob (Public Key)
#print "ut: ",ut



#Step 3. (BOB)
r = gen_vect(n)	  #Bob Generate Secret <r1 .... rn>
x = gen_vect(n)   #Bob Generate Error  <x1 .....xn>
xp = gen_num()    #Bob Generates x' Error

b = A*r + x #Bob sends to Alice (CipherText Preamble)
#print "b:",b
bp = ut*r + xp + bit*np.floor(q/2) #Bob sends to Alice (Payload)
#print "bp",bp



#Step 4. (ALICE)
info = bp - np.transpose(s)*b 


print "------"
print "Alice:",float(info[0][0])," ~~ Bob:",bit*np.floor(q/2)
if (np.round((bit*np.floor(q/2))/info) == bit):
	print "Agreed secret:", bit
print "------"

