from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import optbase as ob

q = symbols('q')

def plotstaff(f,withdiff):
    lamw = lambdify(q, f, modules=['numpy', 'sympy'])
    if withdiff:
        lamdw = lambdify(q, diff(f,q), modules=['numpy', 'sympy'])
        lamddw = lambdify(q, diff(f,q,q), modules=['numpy', 'sympy'])
    lamq = np.linspace(0.1,2.,100)
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

f = Piecewise((0,q>2), (0.25*(2.0 - q)**3, q > 1), (0.25*(2.0 - q)**3 - (1.0 - q)**3, q>0), (0,True))
# plot(f,(q,0,2))
df = diff(f,q)
Fab1D = -2 * df /q
y = integrate(Fab1D,q)
y = addtop(y)
y = integrate(y,q)
y1d = addtop(y)

print("1D")
print("F  : ", y1d)
print("F' : ", diff(y1d,q))
print("F'': ", diff(y1d,q,q))
print('-----------')

Y = Function('Y')
ftmp = list(f.args)
c1, c2 = symbols('c1 c2')
print(Derivative(Y(q),q,q) + Derivative(Y(q),q)/q + 2*df.args[1][0]/q)
# 2 --- 1
ftmp[1] = c1*log(q) + c2 + 1/6 * q**3 - 1.5*q**2 + 6 * q
ftmp[1] = ftmp[1].subs(c1, solve(diff(ftmp[1],q),c1)[0].subs(q,2))
# exit(0)
# print(solve(ftmp[1],c2))
# exit(0)
ftmp[1] = ftmp[1].subs(c2, -(ftmp[1] - c2).subs(q,2.))
# print(N(ftmp[1].subs(q,2)))
# exit(0)
# print(N(ftmp[1].subs(q,2)))
# plotstaff(ftmp[1])
# print(q * Derivative(Y(q),q,q) + Derivative(Y(q),q) + 2*df.args[2][0])
# 1 --- 0
ftmp[2] = c1 * log(q) + c2 - 0.5*q**3 + 1.5*q**2
ftmp[2] = ftmp[2].subs(c1, solve(diff(ftmp[2],q) - diff(ftmp[1],q),c1)[0].subs(q,1))
ftmp[2] = ftmp[2].subs(c2, -(ftmp[2] - ftmp[1] - c2).subs(q,1))

# ftmp[2] = ftmp[2].subs(c1, -(diff(ftmp[2],q) - c1).subs(q,1.) + ftmp[1].subs(q,1.))
# ftmp[2] = ftmp[2].subs(c2, -(ftmp[2] - c2).subs(q,1.) + ftmp[1].subs(q,1.))
# plotstaff(ftmp[2])
ftmp[1] = (ftmp[1], f.args[1][1])
ftmp[2] = (ftmp[2], f.args[2][1])
y2d = Piecewise(*ftmp)

print("2D")
print("F  : ", y2d)
print("F' : ", diff(y2d,q))
print("F'': ", diff(y2d,q,q))
print('-----------')

# print('------------')
# print(Derivative(y2d(q),q,q) + 2/q*Derivative(y2d(q),q) + 2*df.args[1][0]/q)

ftmp = list(f.args)
c1, c2 = symbols('c1 c2')
# 2 --- 1
# print(Derivative(Y(q),q,q) + 2/q*Derivative(Y(q),q) + 2*df.args[1][0]/q)
ftmp[1] = -c1/q + c2 + 0.125 * q**3 - q**2 + 3 * q
ftmp[1] = ftmp[1].subs(c1, solve(diff(ftmp[1],q),c1)[0].subs(q,2))
ftmp[1] = ftmp[1].subs(c2, -(ftmp[1] - c2).subs(q,2.))
# # 1 --- 0
# print(Derivative(Y(q),q,q) + 2/q*Derivative(Y(q),q) + 2*df.args[2][0]/q)
ftmp[2] = -c1/q + c2 - 0.375*q**3 + q**2
ftmp[2] = ftmp[2].subs(c1, solve(diff(ftmp[2],q) - diff(ftmp[1],q),c1)[0].subs(q,1))
ftmp[2] = ftmp[2].subs(c2, -(ftmp[2] - ftmp[1] - c2).subs(q,1))
# print((ftmp[2]-ftmp[1]).subs(q,1))
ftmp[1] = (ftmp[1], f.args[1][1])
ftmp[2] = (ftmp[2], f.args[2][1])
y3d = Piecewise(*ftmp)

print("3D")
print("F  : ", y3d)
print("F' : ", diff(y3d,q))
print("F'': ", diff(y3d,q,q))
print('-----------')

plotstaff(y1d, True)
plotstaff(y2d, True)
plotstaff(y3d, True)
wlist1d, ok = ob.differentiatekernel(y1d ,1)
wlist2d, ok = ob.differentiatekernel(y2d ,1)
wlist3d, ok = ob.differentiatekernel(y3d ,1)
print(ob.calculatenorms(wlist1d))
print(ob.calculatenorms(wlist2d))
print(ob.calculatenorms(wlist3d))
