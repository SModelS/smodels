	QVAR G(DIM,DIM), Gnorm(DIM)
	integer perm(DIM)

#define XSetup(G) XLUDecompose(n, G,DIM, Gnorm, perm)
#define IN(i) in(perm(i))





	QVAR G(DIM,DIM), Ginv(DIM,DIM)

#define XSetup(G) XLUDecompose(n, G,DIM, Ginv,DIM)
#define IN(i) in(i)




sign of permutation?

