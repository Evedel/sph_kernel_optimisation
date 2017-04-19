import optbase as ob
from sympy import *
import scipy.optimize as spo
from datetime import datetime as dt
import numpy as np
import os as os
import random
import matplotlib.pyplot as plt
import matplotlib.lines as mln


q = symbols('q')
A = symbols('A0:10')
Counter = 0
HiddenCounter = 0

Makepath = "/Users/sergeibiriukov/_git/moca_study/fortran/freeze"
Runnerpath = "/Users/sergeibiriukov/_git/moca_study/fortran/freeze/run-diff-err-h.sh"

minimizer_kwargs = {"method": "BFGS"}


# -------------------------
# Exponent with integration
Basefunction = A[0]*exp(-A[1]*(q + A[2])**2) + A[3] * (1. - q/2.) + A[4]
X0 = [1., 1., -1., 1., 1.]
R = 2.
MathEngine = 0
TransformType = 0
# ------------------------------
# Piecewise with differentiation
# Basefunction = Piecewise(\
#                 ((R - q)**5 + A[2]*(A[1] - q)**5 + A[3]*(A[0] - q)**5, q < A[0]),\
#                 ((R - q)**5 + A[2]*(A[1] - q)**5, q < A[1]),\
#                 ((R - q)**5, q < R),\
#                 (0, True))
# A0 = random.uniform(0., 1.)
# X0 = [A0, 2*A0, random.uniform(0., 10.),random.uniform(0., 10.)]
# R = 2.
# MathEngine = 1
# TransformType = 1
# -------------------------------
# simple square as second derivative
# Basefunction = A[0]*(q - A[1])**2
# X0 = [1., 1.]
# R = 2.
# MathEngine = 1
# TransformType = 0
# --------------------------------
# Runge function
# Basefunction = A/(B + C*(q + D)**2) + E * (1. - q/2.) + F
# X0 = [1., 1., 25., 0., 1., 1., 1./sqrt(pi), 1./pi, 1./(pi)**(3./2.)]

MathEnStr = "Mathematica" if MathEngine == 0 else "SymPy"
TransfStr = "Integration" if TransformType == 0 else "Differentiation"
DimCase = [1, 2, 3]
# DimCase = [1]
EstimateColumn = 3
AimValue = 1.
TaskType = ['chi-graddiv']
# TaskType = ['chi-graddiv', 'chi-laplace']
LX = len(X0)

def optfunc(X):
    global Counter
    global HiddenCounter

    tic = dt.now()
    Xstr = ""
    for xi in range(LX):
        if (X[xi] != 0.):
            Xstr = Xstr + "{0:.8f}".format(N(X[xi])) + " "
        else:
            Xstr = Xstr + "{0:.8f}".format(X[xi]) + " "
    Xstr = Xstr[0:-2]
    cost = ob.BrokenPenalty
    # if ((cB < cA) and (cB < 2.) and (cA >= 0.) and (cB >= 0.)):
    w = Basefunction
    for xi in range(LX):
        w = w.subs(A[xi],X[xi])
    w = w - w.subs(q,2.)
    if TransformType == 0:
        klist, ok = ob.integratekernel(w, MathEngine)
    elif TransformType == 1:
        klist, ok = ob.differentiatekernel(w, MathEngine)
    if ok:
        if HiddenCounter == 0:
            print('Saving this info for itteration %d' %Counter)
            lamw = lambdify(q, klist[0], modules=['numpy', 'sympy'])
            lamdw = lambdify(q, klist[1], modules=['numpy', 'sympy'])
            lamd2w = lambdify(q, klist[2], modules=['numpy', 'sympy'])
            lamq = np.linspace(0.,2.,100)
            lw, ldw, ld2w = [], [], []
            for i in range(len(lamq)):
                lw.append(lamw(lamq[i]))
                ldw.append(lamdw(lamq[i]))
                ld2w.append(lamd2w(lamq[i]))
            leg = []
            plt.clf()
            plt.plot(lamq, lw, color='blue')
            leg.append(mln.Line2D([], [], color='blue', label='W'))
            plt.plot(lamq, ldw, color='red')
            leg.append(mln.Line2D([], [], color='red', label='dW'))
            plt.plot(lamq, ld2w, color='green')
            leg.append(mln.Line2D([], [], color='green', label='ddW'))
            plt.legend(handles=leg)
            plt.savefig("./plots/plt-%04d.pdf" %Counter, format='pdf', close=True, verbose=True)
        cost = float(ob.solveproblem(klist, X, Counter, 0, Makepath, TaskType, Runnerpath, DimCase, EstimateColumn, AimValue))
    else:
        print("<*><*><*>\nIntegration failed\nFunction: %s\nEngine: %s\n<*><*><*>\n" %(w,MathEnStr))
    tac = dt.now()
    print("Gen : %s; Cost: %.8f; Time: %s" %(Xstr, cost, (tac-tic)))
    HiddenCounter += 1
    if HiddenCounter == (len(X0)+2):
        Counter += 1
        HiddenCounter = 0
    return cost

def main():
    if not os.path.exists("./plots"):
        os.mkdir("./plots")
    ret = spo.basinhopping(optfunc, X0, minimizer_kwargs=minimizer_kwargs, niter=200)
    print("global minimum: x = %s, f(x0) = %s" % (ret.x, ret.fun))

main()
