import numpy as np
import numpy.random as r
import numpy.linalg as la
import matplotlib.pyplot as plt


def evaluate_function(A, x, mi, n, p, X, y):
    resid = np.dot(X,x)-y;
    sig = 0
    for i in xrange(p):
        sig = sig + la.norm(multiply_Ai(A, x, i, mi, n))

    return 0.5*np.dot(resid,resid) + sig 

def project_2_norm_ball(vec, t):
    vnorm = np.dot(vec, vec)
    if (vnorm > t):
        proj = vec*t/vnorm
    else:
        proj = vec
    return proj

def multiply_Ai(A, vec, i, mi, n):
    return np.dot(A[:, (n*i):(n*(i+1))], vec)

def multiply_Ai_t(A, vec, i, mi, n):
    # print vec.shape
    # print A.shape
    Ait = (A[:, (n*i):(n*(i+1))]).transpose()
    # Ai = A[:, (n*i):(n*(i+1))]
    # print Ait.shape
    res = Ait.dot(vec)
    # print res.shape
    return res

def sigma_Ait_lami(A, lam, mi, n, p):
    res = np.zeros(n)
    for i in xrange(p):
        mult = multiply_Ai_t(A, lam[:,i], i, mi, n)
        # print mult.shape
        # print res.shape
        res = res + mult

    return res

def sigma_Ai_x(A, x, mi, n, p):
    res = np.zeros(mi)
    for i in xrange(p):
        mult = multiply_Ai(A, x, i, mi, n)
        # print mult.shape
        # print res.shape
        res = res + mult

    return res

def calculate_xhat(A, lam, mi, n, p, X, y):
    sig = sigma_Ait_lami(A, lam, mi, n, p)
    Xty = np.dot(X.transpose(), y)
    rhs = Xty - sig
    XtX = np.dot(X.transpose(), X)
    return np.linalg.solve(XtX, rhs)


def dual_prox_step(A, lam, mi, n, p, X, y):

    t = 10
    xh = calculate_xhat(A, lam, mi, n, p, X, y)
    for i in xrange(p):
        lup = lam[:, i] + t*multiply_Ai(A, xh, i, mi, n)
        lam[:,i] = project_2_norm_ball(lup, 1)

    return lam, xh


def calculate_gradient(A, x, mi, n, p, X, y):
    gradf = np.dot(X.transpose(), np.dot(X,x)-y)
    gradh = np.zeros(n)
    for i in xrange(p):
        val = multiply_Ai(A, x, i, mi, n)
        twonorm = la.norm(val)
        if (twonorm == 0):
            twonorm = 1
        val = multiply_Ai_t(A, val, i, mi, n)
        gradh = gradh + val/twonorm

    return gradf + gradh

def gd_step(A, x, mi, n, p, X, y):

    t = 1e-4
    df = calculate_gradient(A, x, mi, n ,p, X, y)
    xnew = x - t*df
    return xnew

def main(T=int(20)):
    X = np.genfromtxt("dual_prox_grad_example/X.csv",delimiter=',')
    A = np.genfromtxt("dual_prox_grad_example/A.csv",delimiter=',')
    y = np.genfromtxt("dual_prox_grad_example/y.csv",delimiter=',')

    print X.shape
    print A.shape
    print y.shape

    p = 500
    mi = 10
    n = 1000

##################### Gradient Descent #####################
    print "PERFORMING GRADIENT DESCENT"
    x = np.zeros(n)
    error_gd = []
    for t in xrange(T):
        # update
        print t
        x = gd_step(A, x, mi, n, p, X, y)
        
        # record error and l1 norm
        fval = evaluate_function(A, x, mi, n, p, X, y)
        error_gd.append(fval)

##############################################


##################### Dual Prox #####################
    print "PERFORMING DUAL PROX DESCENT"
    lam = np.zeros((mi, p))
    error_dp = []
    for t in xrange(T):
        # update
        print t
        lam, x = dual_prox_step(A, lam, mi, n, p, X, y)
        
        # record error and l1 norm
        fval = evaluate_function(A, x, mi, n, p, X, y)
        error_dp.append(fval)

##############################################







    plt.clf()
    plt.plot(error_dp, label='Dual Prox')
    plt.plot(error_gd, label='GD')
    # plt.plot(error_svrg, label='SVRG')
    plt.title('Dual Prox Example')
    plt.legend()
    plt.xlabel('Iteration')
    plt.ylabel('Cost Function Value')
    plt.gca().set_yscale('log')
    # plt.gca().set_xscale('log')
    plt.savefig('dual_prox.eps')


if __name__ == "__main__":
    main()