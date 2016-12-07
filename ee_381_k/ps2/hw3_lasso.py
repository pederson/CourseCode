import numpy as np
import numpy.random as rn
import numpy.linalg as la
import matplotlib.pyplot as plt
import sys


def frank_wolfe(x, A, b, t, gam):
    # update x (your code here)

    # step size
    eta = 2.0/(t+2)
    # calculate gradient
    gr = np.dot(A.T,(np.dot(A,x)-b))
    # next step is in the minimum gradient component
    s = np.zeros(x.size)
    s[gr.argmin()] = -gam*np.sign(gr[gr.argmin()])
    dn = s-x
    # update value
    x = x + eta*dn

    if (x!=0).sum() == 0:
        print "Zero vector!"

    #print la.norm(x)
    return x


def subgradient(x, A, b, t, lam, c=1e-3):
    # update x (your code here), set c above
    
    # step size
    eta = c/np.sqrt(t+1)
    # calculate subgradient
    sg = np.dot(A.T,(np.dot(A,x)-b)) + lam*np.sign(x)
    # calculate new value
    x = x - eta*sg

    return x

# add BTLS variants and include them in main/descent below

def descent(update, A, b, reg, T=int(1e4)):
    x = np.zeros(A.shape[1])
    error = []
    l1 = []
    for t in xrange(T):
        # update A (either subgradient or frank-wolfe)
        x = update(x, A, b, t, reg)
        
        # record error and l1 norm
        if (t % 1 == 0) or (t == T - 1):
            error.append(la.norm(np.dot(A, x) - b))
            l1.append(np.sum(abs(x)))

            assert not np.isnan(error[-1])

    return x, error, l1


def main(T=int(1e2)):
    A = np.load("A.npy")
    b = np.load("b.npy")

    # modify regularization parameters below
    x_sg, error_sg, l1_sg = descent(subgradient, A, b, reg=1e-1, T=T)
    x_fw, error_fw, l1_fw = descent(frank_wolfe, A, b, reg=1.5, T=T)
    # add BTLS experiments

    # add plots for BTLS
    plt.clf()
    plt.plot(error_sg, label='Subgradient')
    plt.plot(error_fw, label='Frank-Wolfe')
    plt.title('Error')
    plt.legend()
    plt.xlabel('iteration number')
    plt.ylabel('Error')
    plt.savefig('error.eps')

    plt.clf()
    plt.plot(l1_sg, label='Subgradient')
    plt.plot(l1_fw, label='Frank-Wolfe')
    plt.title("$\ell^1$ Norm")
    plt.legend()
    plt.xlabel('iteration number')
    plt.ylabel('Norm')
    plt.savefig('l1.eps')


if __name__ == "__main__":
    main()
