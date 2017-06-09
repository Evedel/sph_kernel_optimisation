import optbase as ob
import humpskernels as hk

from sympy import *
import scipy.optimize as spo
from datetime import datetime as dt
import numpy as np
import os as os
import random
import matplotlib.pyplot as plt
import matplotlib.lines as mln
import datetime


q = symbols('q')
A = symbols('A0:10')
Counter = 0
HiddenCounter = 0
InMinimaSearch = 0

Makepath = "/Users/sergeibiriukov/_git/moca_study/fortran/freeze"
Runnerpath = "/Users/sergeibiriukov/_git/moca_study/fortran/freeze/run-diff-err-h.sh"

# minimizer_kwargs = {"method": "BFGS"}
minimizer_kwargs = {"method": "L-BFGS-B"}
# minimizer_kwargs = {"method": "Newton-CG"}

# -------------------------
# Exponent with integration
# Basefunction = A[0]*exp(-A[1]*(q + A[2])**2) + A[3] * (1. - q/2.) + A[4]
# # X0 = [1., 1., 0., 0., 0.]
# X0 = [-1.20216020, 1.209727897, -1.63596466, 1.911459369, 0.34557226]
# # X0 = [88.964659192, 3.04592817425, -0.914791070728, 7.01690442862, -0.99975783763]
# R = 2.
# MathEngine = 0
# TransformType = 2
# def validx(x):
#     return True
# ------------------------------
# Piecewise with differentiation
# R = 2.
# Basefunction = Piecewise(\
#                 ((R - q)**5 + A[2]*(A[1] - q)**5 + A[3]*(A[0] - q)**5, q < A[0]),\
#                 ((R - q)**5 + A[2]*(A[1] - q)**5, q < A[1]),\
#                 ((R - q)**5, q < R),\
#                 (0, True))
# # A0 = 0.7
# A0 = random.uniform(0., 1.)
# # X0 = [A0, 2*A0, -2., -4.]
# X0 = [A0, 2*A0, random.uniform(-10., 10.),random.uniform(-10., 10.)]
# MathEngine = 1
# TransformType = 1
# def validx(x):
#     if ((x[0] > 0) and (x[0] < X[1]) and (x[1] < 2)):
#         return True
#     else:
#         return False
# -------------------------------
# simple square as second derivative
# Basefunction = A[0]*(q - A[1])**2 + A[2]
# X0 = [random.uniform(0., 1.), 2., random.uniform(0., 1.)]
# R = 2.
# MathEngine = 0
# TransformType = 2
# def validx(x):
#     return True
# --------------------------------
# Runge function
# Basefunction = A/(B + C*(q + D)**2) + E * (1. - q/2.) + F
# X0 = [1., 1., 25., 0., 1., 1., 1./sqrt(pi), 1./pi, 1./(pi)**(3./2.)]
# --------------------------------
# Cubic spline as function
# R = 2.
# Basefunction = q**2*Piecewise(\
#                 (A[2] * (R - q)**3 - A[1] * (A[0] - q)**3, q < A[0]),\
#                 (A[2] * (R - q)**3, q < R),\
#                 (0, True))
# X0 = [1., 1., 0.25]
# MathEngine = 1
# TransformType = 0
# def validx(x):
#     if ((x[0] > 0) and (x[0] < 2)):
#         return True
#     else:
#         return False
# --------------------------------
# Humps kernels, search only for multiplier power
R = 2.
Basefunction = q
X0 = [0., 1.]
MathEngine = 1
TransformType = 3
def validx(x):
    if ((x[0] >= 0) and (x[1] > 0)):
        return True
    else:
        return False

MathEnStr = "Mathematica" if MathEngine == 0 else "SymPy"

if TransformType == 0:
    TransfStr = "Integration"
elif TransformType == 1:
    TransfStr = "Differentiation"
elif TransformType == 2:
    TransfStr = "Kernel ODE Solution"
elif TransformType == 3:
    TransfStr = "Humps Kernels"
else:
    print("Not implemented: ", TransformType)
    exit(127)

DimCase = [1, 2, 3]
# DimCase = [1]
EstimateColumn = 3
AimValue = 1.
# TaskType = ['chi-graddiv']
# TaskType = ['chi-laplace']
TaskType = ['chi-graddiv', 'chi-laplace']
LX = len(X0)

def GetKernels(X):
    w = Basefunction
    LX = len(X)
    if (len(X0) != LX):
        klist = [0, 0, 0]
        ok = False
    else:
        for xi in range(LX):
            w = w.subs(A[xi],X[xi])
        w = w - w.subs(q,2.)
        if TransformType == 0:
            klist, ok = ob.integratekernel(w, MathEngine)
        elif TransformType == 1:
            klist, ok = ob.differentiatekernel(w, MathEngine)
        elif TransformType == 2:
            klist, ok = ob.solvekernelode(w)
        elif TransformType == 3:
            klist, ok = hk.humpskerlnels(X[0], X[1])
    return klist, ok

def PrintStep(x, f, accepted):
    global InMinimaSearch
    global Counter

    X = x
    LX = len(x)
    klist, ok = GetKernels(X)
    if ok:
        if ((TransformType == 0) or
            (TransformType == 1) or
            (TransformType == 3)):
            lamw = lambdify(q, klist[0], modules=['numpy', 'sympy'])
            lamdw = lambdify(q, klist[1], modules=['numpy', 'sympy'])
            lamd2w = lambdify(q, klist[2], modules=['numpy', 'sympy'])
            lamq = np.linspace(0.01,2.,100)
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
            leg.append(mln.Line2D([], [], color='red', label='W\''))
            plt.plot(lamq, ld2w, color='green')
            leg.append(mln.Line2D([], [], color='green', label='W\'\''))
            plt.legend(handles=leg)
            if TransformType == 3:
                plt.annotate(r'$q^{%s}M4$' %("{0:.3f}".format(X[0]/X[1])), xy=(1.0, -0.1))
            plt.savefig("./plots/plt-%04d.pdf" %Counter, format='pdf', close=True, verbose=True)
        if (TransformType == 2):
            for dim in range(len(DimCase)):
                lamw = lambdify(q, klist[0][dim], modules=['numpy', 'sympy'])
                lamdw = lambdify(q, klist[1][dim], modules=['numpy', 'sympy'])
                lamd2w = lambdify(q, klist[2][dim], modules=['numpy', 'sympy'])
                lamq = np.linspace(0.01,2.,100)
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
                leg.append(mln.Line2D([], [], color='red', label='W\''))
                plt.plot(lamq, ld2w, color='green')
                leg.append(mln.Line2D([], [], color='green', label='W\'\''))
                plt.legend(handles=leg)
                plt.savefig("./plots/plt-%04d-%d.pdf" %(Counter, dim+1), format='pdf', close=True, verbose=True)
        Counter += 1
    else:
        print("Transform went wrong: ", klist, X)
    Xstr = ""
    for xi in range(LX):
        if (X[xi] != 0.):
            if (X[xi] < 0):
                Xstr = Xstr + "{0:.8f}".format(N(X[xi])) + " "
            else:
                Xstr = Xstr + "{0:.9f}".format(N(X[xi])) + " "
        else:
            Xstr = Xstr + "{0:.8f}".format(X[xi]) + " "
    Xstr = Xstr[0:-2]

    coststr = ""

    if f < 0:
        coststr += "{0:.6f}".format(f)
    else:
        coststr += "{0:.7f}".format(f)
    # print(os.getcwd())
    print("%s [ %03d ]; At minimum: %s; Cost: %s; Accepted: %5s; Steps to minima: %d" \
        %(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), Counter, Xstr, coststr, \
        accepted, InMinimaSearch)\
    )
    InMinimaSearch = 0

def optfunc(X):
    global Counter
    global HiddenCounter
    global InMinimaSearch

    InMinimaSearch += 1
    cost = ob.BrokenPenalty

    if validx(X):
        klist, ok = GetKernels(X)
        if ok:
            if (TransformType == 2):
                tmpc = 0.
                for dim in range(len(DimCase)):
                    kdl = [klist[0][dim], klist[1][dim], klist[2][dim]]
                    # print(kdl[0], '\n')
                    # print(kdl[1], '\n')
                    # print(kdl[2])
                    # print('----------')
                    tmpc += float(ob.solveproblem(kdl, X, Counter, 0, Makepath, \
                        TaskType, Runnerpath, [dim + 1], EstimateColumn, AimValue))
                cost = tmpc
            else:
                cost = float(ob.solveproblem(klist, X, Counter, 0, Makepath, \
                    TaskType, Runnerpath, DimCase, EstimateColumn, AimValue))
    return cost

def func2d(x):
    f = np.cos(14.5 * x[0] - 0.3) + (x[1] + 0.2) * x[1] + (x[0] + 0.2) * x[0]
    df = np.zeros(2)
    df[0] = -14.5 * np.sin(14.5 * x[0] - 0.3) + 2. * x[0] + 0.2
    df[1] = 2. * x[1] + 0.2
    return f

def main():
    # TEST
    # PrintStep([1.0, 1.0], func2d([1.0, 1.0]), False)
    # spo.basinhopping(func2d, [1.0, 1.0], minimizer_kwargs={"method":"L-BFGS-B"}, callback=PrintStep)
    # exit(0)
    if not os.path.exists("./plots"):
        os.mkdir("./plots")
    PrintStep(X0, optfunc(X0), False)
    ret = spo.basinhopping(optfunc, X0, minimizer_kwargs=minimizer_kwargs, callback=PrintStep)
    print("global minimum: x = %s, f(x0) = %s" % (ret.x, ret.fun))

main()
