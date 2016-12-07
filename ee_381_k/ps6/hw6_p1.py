import numpy as np
import numpy.random as r
import numpy.linalg as la
import matplotlib.pyplot as plt


def gd_step(x, A, b, t, mu):

    gamma = step_size_const(t)
    gradf = calc_full_gradient(x, A, b, mu)
    xnew = x - gamma*gradf
    return xnew

def sgd_step(x, A, b, t, mu):
    rn = calc_random_int(x, A, b, mu)
    gamma = step_size_const(t)
    gradf = calc_partial_gradient(x, A, b, mu, rn)
    xnew = x - gamma*gradf
    return xnew

def svrg_step(x, A, b, t, mu, full_grad, x_old):

    rn = calc_random_int(x, A, b, mu)
    gamma = step_size_const(t)
    gradfi = calc_partial_gradient(x, A, b, mu, rn)
    gradfold = calc_partial_gradient(x_old, A, b, mu, rn)
    gradf = gradfi - gradfold + full_grad
    xnew = x - gamma*gradf
    return xnew

def evaluate_function(x, A, b, mu):
    resid = np.dot(A,x)-b
    val = np.dot(resid, resid) + mu*np.dot(x,x)
    return val

def step_size_const(t):

    return 1e-7

def step_size_decr(t):
    return 1e-7*0.9**t

def calc_full_gradient(x, A, b, mu):

    resid = np.dot(A,x)-b
    grad = 2.0*np.dot(A.transpose(), resid)
    grad = grad + 2*mu*x
    return grad


def calc_partial_gradient(x, A, b, mu, i):
    grad = 2*A.shape[0]*(np.dot(A[i,:], x)-b[i])*A[i,:].transpose() + 2*mu*x
    return grad

def calc_random_int(x, A, b, mu):
    rn = np.random.randint(0, A.shape[0]-1)
    return rn

def main(T=int(30000)):
    X_train = np.genfromtxt("digits/X_digits_train.csv",delimiter='')
    y_train = np.genfromtxt("digits/y_digits_train.csv",delimiter='')
    X_test = np.genfromtxt("digits/X_digits_test.csv",delimiter='')
    y_test = np.genfromtxt("digits/y_digits_test.csv",delimiter='')

    print X_train.shape

    mu = 1e-1
##################### GD #####################
    x_gd = np.zeros(X_train.shape[1])
    fbest_gd = evaluate_function(x_gd, X_train, y_train, mu)
    xbest_gd = x_gd
    error_gd = []
    for t in xrange(T):
        # update A (CG, SGD, SVRG)
        x_gd = gd_step(x_gd, X_train, y_train, t, mu)
        
        # record error and l1 norm
        if (t % 1 == 0) or (t == T - 1):
            fval = evaluate_function(x_gd, X_train, y_train, mu)
            if (fval < fbest_gd):
                fbest_gd = fval
                xbest_gd = x_gd
            error_gd.append(fval)

            #assert not np.isnan(error1[-1])


    print evaluate_function(x_gd, X_train, y_train, mu)
##############################################


##################### SGD ####################
    x_sgd = np.zeros(X_train.shape[1])
    fbest_sgd = evaluate_function(x_sgd, X_train, y_train, mu)
    xbest_sgd = x_sgd
    error_sgd = []
    for t in xrange(T):
        # update A (CG, SGD, SVRG)
        x_sgd = sgd_step(x_sgd, X_train, y_train, t, mu)
        
        # record error and l1 norm
        if (t % 1 == 0) or (t == T - 1):
            fval = evaluate_function(x_sgd, X_train, y_train, mu)
            if (fval < fbest_gd):
                fbest_sgd = fval
                xbest_sgd = x_sgd
            error_sgd.append(fval)
            #assert not np.isnan(error1[-1])


    print evaluate_function(x_sgd, X_train, y_train, mu)
########################################################


##################### SVRG ####################
    x_svrg = np.zeros(X_train.shape[1])
    fbest_svrg = evaluate_function(x_svrg, X_train, y_train, mu)
    xbest_svrg = x_svrg
    error_svrg = []
    m = 100
    for i in xrange(T/m):
        full_grad = calc_full_gradient(x_svrg, X_train, y_train, mu)
        x_hold = x_svrg
        for t in xrange(m):
            # update A (CG, SGD, SVRG)
            x_svrg = svrg_step(x_svrg, X_train, y_train, t, mu, full_grad, x_hold)
            
            # record error and l1 norm
            if (t % 1 == 0) or (t == T - 1):
                #error_svrg.append(evaluate_function(x_svrg, X_train, y_train, mu))
                fval = evaluate_function(x_svrg, X_train, y_train, mu)
                if (fval < fbest_gd):
                    fbest_svrg = fval
                    xbest_svrg = x_svrg
                error_svrg.append(fval)
                #assert not np.isnan(error1[-1])


    print evaluate_function(x_svrg, X_train, y_train, mu)
########################################################



    plt.clf()
    plt.plot(error_gd, label='GD')
    plt.plot(error_sgd, label='SGD')
    plt.plot(error_svrg, label='SVRG')
    plt.title('Training Set (m=100)')
    plt.legend()
    plt.xlabel('Iteration')
    plt.ylabel('Cost Function Value')
    plt.gca().set_yscale('log')
    plt.gca().set_xscale('log')
    plt.savefig('func_value.eps')


if __name__ == "__main__":
    main()