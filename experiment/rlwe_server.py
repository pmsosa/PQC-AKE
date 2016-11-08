###                                         ###
## Ring Learning With Errors Server (Bob)    ##
##       Pedro Miguel Sosa (pmsosa)          ##
###                                         ###

import socket
import pickle
import timeit
import numpy as np
from numpy.polynomial import polynomial as p
from Crypto.Hash import SHA256
from Crypto.Cipher import AES
from Crypto import Random



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

#1. Listen and Connect
print "Waiting for connection..."
skt = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
skt.bind((TCP_IP, TCP_PORT))
skt.listen(1)

conn, addr = skt.accept()
print 'Connection address:', addr


#2. Get A and b (from Alice)
A = pickle.loads(conn.recv(BUFFER_SIZE)) #In reality A would be KNOWN by both users

kex_start = timeit.default_timer()

b = pickle.loads(conn.recv(BUFFER_SIZE))

#3. Calculate and send b' = A * s' + e'

s = gen_poly(n,q)
e = gen_poly(n,q)

bB = p.polymul(A,s)%q
bB = np.floor(p.polydiv(s,hlpr)[1])
bB = p.polyadd(bB,e)%q 

conn.send(pickle.dumps(bB))


#5. Calculate and send u
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
    
conn.send(pickle.dumps(u))    

#4. Calculate Shared Secret (And do error correction)
sharedBob = np.floor(p.polymul(s,b)%q)
sharedBob = np.floor(p.polydiv(sharedBob,hlpr)[1])%q

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

    
kex_end = timeit.default_timer()
print "KEX Time (Bob): ",kex_end-kex_start    
print "Shared Secret", sharedBob

##############################>----KEX ENDS
BS = 16
pad = lambda s: s + (BS - len(s) % BS) * chr(BS - len(s) % BS) 
unpad = lambda s : s[0:-ord(s[-1])]


enc_key = SHA256.new(pickle.dumps(sharedBob)).hexdigest().decode("hex")
print enc_key.encode("hex")

iv = Random.new().read(AES.block_size)
obj = AES.new(enc_key, AES.MODE_CBC, iv)

while 1:
    ciphertext = conn.recv(BUFFER_SIZE)
    message = unpad(obj.decrypt(ciphertext))
    print message
    if (message == "quit"): break;


conn.close()