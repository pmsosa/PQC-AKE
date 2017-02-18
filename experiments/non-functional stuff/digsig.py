


class Signature:

	def __init__(self):
		return 0;

	def __gen_poly(self,n,q):
		return 0;

	def __hash(self,m):
		return 0;

	def __check_addition(s1,s2,t):
		return 0;

	#SigKeyGen --> (Ks,Kv) = ((f,g),h)
	def SigKeyGen(self):
		f = self.__gen_poly(self.n,self.q) 
		g = self.__gen_poly(self.n,self.q) 
		h = (p.polydiv(f,g)[0]).astype(int)
		return ((f,g),h)

	#Sig((f,g),m) -> o=(s1,m)
	def Sig(self,(f,g),m):
		t = self.__hash(m)
		s1 = self.__gen_poly(self.n,self.q) 
		s2 = self.__gen_poly(self.n,self.q) 
		self.__check_addition(s1,s2,t)
		return (s1,m)

	def Ver(self,h,(s1,m)):
		t = self.__hash(m)

		s2 = p.polymult(h,s1)
		s2 = p.polysub(t,s2)

		### if ||(s1,s2) || < B : 
		###  return 0
		### else:
		###  return -1

		return 0