###                                         ###
## Ring Learning With Errors Client (Alice)  ##
##       Pedro Miguel Sosa (pmsosa)          ##
###                                         ###

import socket
import pickle
import numpy as np
from numpy.polynomial import polynomial as p


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



TCP_IP = '127.0.0.1'
TCP_PORT = 5005
BUFFER_SIZE = 50000

n = 1024        #n and q, based on D. Stebila's talk
q = 2**32-1

hlpr = [1] + [0] * (n-1) + [1]


skt = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

#1. CONNECT 
skt.connect((TCP_IP, TCP_PORT))

#2. SEND A
A = np.floor(np.random.random(size=(n))*q)%q
A = np.floor(p.polydiv(A,hlpr)[1])

skt.send(pickle.dumps(A))

#3. SEND b= (A * s + e)
s = gen_poly(n,q)
e = gen_poly(n,q)

b = p.polymul(A,s)%q
b = np.floor(p.polydiv(s,hlpr)[1])
b = p.polyadd(b,e)%q

skt.send(pickle.dumps(b))


#4. RECIEVE b' = (s' * A + e') and U (for error correction)
bB = pickle.loads(skt.recv(BUFFER_SIZE))
u  = pickle.loads(skt.recv(BUFFER_SIZE))

#5. Calculate Shared Secret & Error Correct
sharedAlice = np.floor(p.polymul(s,bB)%q)
sharedAlice = np.floor(p.polydiv(sharedAlice,hlpr)[1])%q

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

print "Shared Secret", sharedAlice


skt.close()
