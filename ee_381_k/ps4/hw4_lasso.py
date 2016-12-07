import numpy as np
import numpy.random as rn
import numpy.linalg as la
import matplotlib.pyplot as plt
import sys


def frank_wolfe(x, A, b, t, gam, l, y):
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
    return x, 0, 0


def subgradient(x, A, b, t, lam, l, y, c=1e-3):
    # update x (your code here), set c above
    
    # step size
    eta = c/np.sqrt(t+1)
    # calculate subgradient
    sg = np.dot(A.T,(np.dot(A,x)-b)) + lam*np.sign(x)
    # calculate new value
    x = x - eta*sg

    return x, 0, 0


def proximal_gradient(x, A, b, t, lam, l, y, c=1e-3):

    # step size
    eta = c

    # take gradient step
    x = x - eta*np.dot(A.T,(np.dot(A,x)-b))

    # use the prox operator
    x = np.sign(x)*(abs(x) - lam)

    return x, 0, 0


def accelerated_gradient(x, A, b, t, lam, lprev, yprev, c=1e-3):

    # step size
    eta = c

    # lambda params
    l = (1+np.sqrt(1+4*lprev))/2
    gamma = (1-lprev)/l

    # take gradient step
    x = x - eta*np.dot(A.T,(np.dot(A,x)-b))

    # use the prox operator
    y = np.sign(x)*(abs(x) - lam)

    # take the weighted sum
    x = (1-gamma)*y + gamma*yprev

    return x, l, y


# add BTLS variants and include them in main/descent below

def descent(update, A, b, reg, T=int(1e4)):
    x = np.zeros(A.shape[1])
    error = []
    l1 = []
    l = 0;
    y = x;
    for t in xrange(T):
        # update A (either subgradient or frank-wolfe)
        x, l, y = update(x, A, b, t, reg, l, y)
        
        # record error and l1 norm
        if (t % 1 == 0) or (t == T - 1):
            error.append(la.norm(np.dot(A, x) - b))
            l1.append(np.sum(abs(x)))

            assert not np.isnan(error[-1])

    return x, error, l1


def main(T=int(1e2)):
    A = np.load("A_test.npy")
    b = np.load("b_test.npy")

    # modify regularization parameters below
    x_sg, error_sg, l1_sg = descent(subgradient, A, b, reg=1e-1, T=T)
    x_fw, error_fw, l1_fw = descent(frank_wolfe, A, b, reg=1.5, T=T)
    x_pg, error_pg, l1_pg = descent(proximal_gradient, A, b, reg=1e-2, T=T)
    x_ap, error_ap, l1_ap = descent(accelerated_gradient, A, b, reg=1e-2, T=T)
    # add BTLS experiments

    # add plots for BTLS
    plt.clf()
    plt.plot(error_sg, label='Subgradient')
    plt.plot(error_fw, label='Frank-Wolfe')
    plt.plot(error_pg, label='Proximal-Gradient')
    plt.plot(error_ap, label='Accelerated Proximal-Gradient')
    plt.title('Error')
    plt.legend()
    plt.xlabel('iteration number')
    plt.ylabel('Error')
    plt.savefig('error.eps')

    plt.clf()
    plt.plot(l1_sg, label='Subgradient')
    plt.plot(l1_fw, label='Frank-Wolfe')
    plt.plot(l1_pg, label='Proximal-Gradient')
    plt.plot(l1_ap, label='Accelerated Proximal-Gradient')
    plt.title("$\ell^1$ Norm")
    plt.legend()
    plt.xlabel('iteration number')
    plt.ylabel('Norm')
    plt.savefig('l1.eps')


if __name__ == "__main__":
    main()
