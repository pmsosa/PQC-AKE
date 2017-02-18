import numpy as np
from numpy.polynomial import polynomial as p

n = 1024
q = 2**32-1
hlpr = [1] + [0] * (n-1) + [1]


def gen_poly(n,q):
    global hlpr
    l = 0 #Gamma Distribution Location (Mean "center" of dist.)
    poly = np.floor(np.random.normal(l,size=(n)))
    while (len(poly) != n):
        poly = np.floor(np.random.normal(l,size=(n)))
        poly = np.floor(p.polydiv(poly,hlpr)[1]%q)
    return poly



#Generate A
A = np.floor(np.random.random(size=(n))*q)%q

A = np.floor(p.polydiv(A,hlpr)[1])

#Alice (Secret & Error)
sA = gen_poly(n,q)
eA = gen_poly(n,q)

bA = p.polymul(A,sA)%q
bA = np.floor(p.polydiv(sA,hlpr)[1])
bA = p.polyadd(bA,eA)%q


#Bob
sB = gen_poly(n,q)
eB = gen_poly(n,q)


bB = p.polymul(A,sB)%q
bB = np.floor(p.polydiv(sB,hlpr)[1])
bB = p.polyadd(bB,eB)%q


#Shared Secret
#Alice
sharedAlice = np.floor(p.polymul(sA,bB)%q)
sharedAlice = np.floor(p.polydiv(sharedAlice,hlpr)[1])%q  #TODO FIX THIS HAS TO BE DIVED BY HELPER
sharedBob = np.floor(p.polymul(sB,bA)%q)
sharedBob = np.floor(p.polydiv(sharedBob,hlpr)[1])%q

#Error Rounding
#--Bob
u = np.asarray([0] * n)
i = 0

while (i < len(u)):
	if (len(bB) <= i): break;
	if (int(bB[i]/(q/4)) == 0): u[i] = 0
	elif (int(bB[i]/(q/2)) == 0): u[i] = 1
	elif (int(bB[i]/(3*q/4)) == 0): u[i] = 0
	elif (int(bB[i]/(q)) == 0): u[i] = 1
	else:
		print "error! (1)"
	i+=1


i = 0
while (i < len(u)):
	#Region 0 (0 --- q/4 and q/2 --- 3q/4)
	if (u[i] == 0):
		if (sharedBob[i] >= q*0.125 and sharedBob[i] < q*0.625):
			sharedBob[i] = 1
		else:
			sharedBob[i] = 0


	#Region 1 (q/4 --- q/2 and 3q/4 --- q)
	elif (u[i] == 1):
		if (sharedBob[i] >= q*0.875 and sharedBob[i] < q*0.375):
			sharedBob[i] = 0
		else:
			sharedBob[i] = 1

	else:
		print "error! (2)"

	i += 1

#--Alice
i = 0
while (i < len(u)):
	#Region 0 (0 --- q/4 and q/2 --- 3q/4)
	if (u[i] == 0):
		if (sharedAlice[i] >= q*0.125 and sharedAlice[i] < q*0.625):
			sharedAlice[i] = 1
		else:
			sharedAlice[i] = 0


	#Region 1 (q/4 --- q/2 and 3q/4 --- q)
	elif (u[i] == 1):
		if (sharedAlice[i] >= q*0.875 and sharedAlice[i] < q*0.375):
			sharedAlice[i] = 0
		else:
			sharedAlice[i] = 1

	else:
		print "error! (3)"
	i += 1





#
print "A:",len(A),"|",A
print "\n-Alice---"
print " s:",len(sA),"|",sA
print " e:",len(eA),"|",eA
print " b:",len(bA),"|",bA
print "\n-Bob---"
print " s':",len(sB),"|",sB
print " e':",len(eB),"|",eB
print " b':",len(bB),"|",bB
print " u :",len(u),"|",u
print "\n"
print "Shared Secret Alice:",len(sharedAlice),"|",sharedAlice
print "Shared Secret Bob:",len(sharedBob),"|",sharedBob


print "\n\n--Verification--"
i = 0
while (i < len(sharedBob)):
	if (sharedAlice[i] != sharedBob[i]):
		print "Error at index",i
	i+=1