import numpy as np
import numpy.random as r
import numpy.linalg as la
import matplotlib.pyplot as plt


def conjugate_gradient(x, A, b, t, p_old, r_old):
	Ap = np.dot(A, p_old)
	alpha = -np.dot(r_old, p_old)/np.dot(p_old, Ap)
	xnew = x + alpha*p_old
	r_new = np.dot(A,xnew) - b
	p_new = -r_new + p_old*np.dot(r_new, Ap)/np.dot(p_old, Ap)


	return xnew, p_new, r_new


def main(T=int(50)):
    M1 = np.genfromtxt("ConjugateGradient/M1.csv",delimiter=',')
    b1 = np.genfromtxt("ConjugateGradient/b1.csv",delimiter=',')
    M2 = np.genfromtxt("ConjugateGradient/M2.csv",delimiter=',')
    b2 = np.genfromtxt("ConjugateGradient/b2.csv",delimiter=',')
    xtrue = np.genfromtxt("ConjugateGradient/x.csv",delimiter=',')


    x1 = np.zeros(M1.shape[1])
    r = np.dot(M1, x1) - b1
    p = -r
    error1 = []
    for t in xrange(M1.shape[1]):
        # update A (either subgradient or frank-wolfe)
        x1, p, r = conjugate_gradient(x1, M1, b1, t, p, r)
        
        # record error and l1 norm
        if (t % 1 == 0) or (t == T - 1):
            error1.append(la.norm(np.dot(x1-xtrue, r)))

            #assert not np.isnan(error1[-1])



    x2 = np.zeros(M2.shape[1])
    r = np.dot(M2, x2) - b2
    p = -r
    error2 = []
    for t in xrange(M1.shape[1]):
        # update A (either subgradient or frank-wolfe)
        x2, p, r = conjugate_gradient(x2, M2, b2, t, p, r)
        
        # record error
        if (t % 1 == 0) or (t == T - 1):
            error2.append(la.norm(np.dot(x2-xtrue, r)))

            #assert not np.isnan(error2[-1])

    # add plots for BTLS
    plt.clf()
    plt.plot(error1, label='M1')
    plt.plot(error2, label='M2')
    plt.title('Error')
    plt.legend()
    plt.gca().set_yscale('log')
    plt.savefig('error.eps')


if __name__ == "__main__":
    main()