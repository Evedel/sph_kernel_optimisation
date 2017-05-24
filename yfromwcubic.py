from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from kernoptim import optbase as ob

q = symbols('q')

def plotstaff(f,withdiff):
    lamw = lambdify(q, f, modules=['numpy', 'sympy'])
    if withdiff:
        lamdw = lambdify(q, diff(f,q), modules=['numpy', 'sympy'])
        lamddw = lambdify(q, diff(f,q,q), modules=['numpy', 'sympy'])
    lamq = np.linspace(0.5,2.,100)
    lw, ldw, lddw = [], [], []
    for i in range(len(lamq)):
        lw.append(lamw(lamq[i]))
        if withdiff:
            ldw.append(lamdw(lamq[i]))
            lddw.append(lamddw(lamq[i]))
    plt.plot(lamq, lw, color='blue', label="W")
    if withdiff:
        plt.plot(lamq, ldw, color='red', label="W\'")
        plt.plot(lamq, lddw, color='green', label="W\"")
    plt.legend(loc=1)
    plt.show()

def addtop(f):
    if (isinstance(f, Piecewise)):
        ftmp = list(f.args)
        ftmp[1] = (f.args[1][0] - f.args[1][0].subs(q, 2.), f.args[1][1])
        ftmp[2] = (f.args[2][0] - f.args[2][0].subs(q,1.) + ftmp[1][0].subs(q,1.), f.args[2][1])
        f = Piecewise(*ftmp)
    return f


R = 2.

W = Piecewise((0,q>2), (0.25*(2.0 - q)**3, q > 1), (0.25*(2.0 - q)**3 - (1.0 - q)**3, q>0), (0,True))
print("W: ")
print(W)
print('-----------')
# plot(f,(q,0,2))
DDY = W
DY = integrate(DDY,q)
DY = addtop(DY)
Y = integrate(DY,q)
Y1d = addtop(Y)

print("1D")
print("F  : ", Y1d)
print("F' : ", diff(Y1d,q))
print("F'': ", diff(Y1d,q,q))
print('-----------')

Y = Function('Y')
ftmp = list(W.args)
c1, c2 = symbols('c1 c2')
# -----------------------------------------------------------------------
# print(Derivative(Y(q),q,q) + Derivative(Y(q),q)/q - W.args[1][0])
# exit(0)
# 2 --- 1
ftmp[1] = c1*log(q) + c2 - 0.01 * q**5 + 0.09375 * q**4 - 1./3. * q**3 + 0.5 * q**2
ftmp[1] = ftmp[1].subs(c1, solve(diff(ftmp[1],q),c1)[0].subs(q,2))
ftmp[1] = ftmp[1].subs(c2, -(ftmp[1] - c2).subs(q,2.))
# print(diff(ftmp[1],q).subs(q,2.))
# print(ftmp[1].subs(q,2.))
# exit(0)
# -----------------------------------------------------------------------
# print(Derivative(Y(q),q,q) + Derivative(Y(q),q)/q - W.args[2][0])
# exit(0)
# 1 --- 0
ftmp[2] = c1 * log(q) + c2 + 0.03 * q**5 - 0.09375 + 0.25 * q**2
ftmp[2] = ftmp[2].subs(c1, solve(diff(ftmp[2],q) - diff(ftmp[1],q),c1)[0].subs(q,1))
# ftmp[2] = ftmp[2].subs(c1, solve(diff(ftmp[2],q,q) - diff(ftmp[1],q,q),c1)[0].subs(q,1))
ftmp[2] = ftmp[2].subs(c2, -(ftmp[2] - ftmp[1] - c2).subs(q,1))
# print((diff(ftmp[2],q,q) - diff(ftmp[1],q,q)).subs(q, 1.))
# print((diff(ftmp[2],q) - diff(ftmp[1],q)).subs(q, 1.))
# print((ftmp[2] - ftmp[1]).subs(q, 1.))
# exit(0)
ftmp[1] = (ftmp[1], W.args[1][1])
ftmp[2] = (ftmp[2], W.args[2][1])
Y2d = Piecewise(*ftmp)
# plotstaff(Y2d, True)
# exit(0)

print("2D")
print("F  : ", Y2d)
print("F' : ", diff(Y2d,q))
print("F'': ", diff(Y2d,q,q))
print('-----------')

ftmp = list(W.args)
c1, c2 = symbols('c1 c2')
# 2 --- 1
# print(Derivative(Y(q),q,q) + 2/q*Derivative(Y(q),q) - W.args[1][0])
# exit(0)
ftmp[1] = -c1/q + c2 - 1./12. * q**5 + 0.075 * q**4 - 0.25 * q**3 + 1./3. * q**2
ftmp[1] = ftmp[1].subs(c1, solve(diff(ftmp[1],q),c1)[0].subs(q,2))
ftmp[1] = ftmp[1].subs(c2, -(ftmp[1] - c2).subs(q,2.))
# # 1 --- 0
# print(Derivative(Y(q),q,q) + 2/q*Derivative(Y(q),q) - W.args[2][0])
# exit(0)
ftmp[2] = -c1/q + c2 + 0.025 * q**5 - 0.075 * q**4 + 1./6. * q**2
ftmp[2] = ftmp[2].subs(c1, solve(diff(ftmp[2],q) - diff(ftmp[1],q),c1)[0].subs(q,1))
ftmp[2] = ftmp[2].subs(c2, -(ftmp[2] - ftmp[1] - c2).subs(q,1))
# print((ftmp[2]-ftmp[1]).subs(q,1))
ftmp[1] = (ftmp[1], W.args[1][1])
ftmp[2] = (ftmp[2], W.args[2][1])
Y3d = Piecewise(*ftmp)
print("3D")
print("F  : ", Y3d)
print("F' : ", diff(Y3d,q))
print("F'': ", diff(Y3d,q,q))
print('-----------')

ty = Y3d.args[1][0]
tw = W.args[1][0]
print(ty)
print(tw)
print(simplify(diff(ty,q,q) + 2 * diff(ty,q)/q - tw))
ty = Y3d.args[2][0]
tw = W.args[2][0]
print(ty)
print(tw)
print(simplify(diff(ty,q,q) + 2 * diff(ty,q)/q - tw))
exit(0)

plotstaff(Y1d, True)
plotstaff(Y2d, True)
plotstaff(Y3d, True)

# wlist1d, ok = ob.differentiatekernel(Y1d, 1)
# wlist2d, ok = ob.differentiatekernel(Y2d, 1)
# wlist3d, ok = ob.differentiatekernel(Y3d, 1)
# print(ob.calculatenorms(wlist1d)[0], ob.calculatenorms(wlist2d)[1], ob.calculatenorms(wlist3d)[2])
