from sympy import *
import optbase as ob
import humpskernels as hk


q = symbols('q')

wker = (sin(pi/2*q)**4)/(pi/2*q)**4
wlist, ok = ob.differentiatekernel(wker ,1)
wnorms = ob.calculatenorms(wlist[0])
# print(ob.printkernel(wlist, 2.0, wnorms, 'sinc4'))
# ob.calcdphidh(wlist[0])

klist, ok = hk.humpskerlnels(2.11, 1.)
k2list, ok = hk.humpskerlnels(2., 1.)
k3list, ok = hk.humpskerlnels(3., 1.)
knorms = ob.calculatenorms(klist[0])
k2norms = ob.calculatenorms(k2list[0])
k3norms = ob.calculatenorms(k3list[0])

w2ker = 88.964659192*exp(-3.04592817425*(q -0.914791070728)**2) + 7.01690442862 * (1. - q/2.) -0.99975783763
w2list, ok = ob.integratekernel(w2ker, 0)
w2norms = ob.calculatenorms(w2list[0])

m4 = Piecewise(\
        (0., q < 0.),
        (0.25 * (2. - q)**3 - (1. - q)**3, q < 1.),\
        (0.25 * (2. - q)**3, q < 2.),\
        (0, True))
m4list, ok = ob.differentiatekernel(m4, 1)
m4norm = ob.calculatenorms(m4list[0])

# print(ob.printkernel(wlist, 2., norms, 'sinc4'))
# print(wlist)
# print(norms)
# wker = Piecewise((1, q < 1.33333333), (1.5 * (2 - q), q < 2.), (0, True))
# wlist, ok = ob.integratekernel(wker,1)
# plot(wlist[0], (q, 0, 2))
# plot(wlist[1], (q, 0, 2))

# Piecewise(
#     (0, q < 0.0), \
#     (4.5*q - 3.0, q < 1.0), \
#     (-1.5*q + 3.0, q < 2.0))
ddmm4 = Piecewise( \
    (0, q < 0.0), \
    (3.0*q - 2.0, q < 1.0), \
    (2.0 - q, q < 2.0), \
    (0, True))

ddmm4list, ok = ob.integratekernel(ddmm4, 1)
ddmm4norm = ob.calculatenorms(ddmm4list[0])

w = Piecewise(\
    (0.159524 + 0.166667*(-1.36667*q + 0.5*q**2), q < 1.0), \
    (0.152381 - 0.177778*q + 0.111111*q**4 - 0.1*q**5 + 0.0333333*q**6 - 0.00396825*q**7, q < 2.0),\
    (0,True))
wl, ok = ob.differentiatekernel(w)
# nl = ob.calculatenorms(wl[0])
ob.tabulatekernel(wl, '.')
# print(ob.printkernel(wl, 2.0, nl, 'm4--'))
# plot(wl[0], wl[1], wl[2], (q,0.1,3))
exit(0)
#
# plot(ddmm4list[0],
#     ddmm4list[1],
#     ddmm4list[2],
#     (q,0,2))

# print(m4list[2])
dim = 1

plot(wnorms[dim-1]*(wlist[2] + (dim-1)*wlist[1]/q),
    w2norms[dim-1]*(w2list[2] + (dim-1)*w2list[1]/q),
    knorms[dim-1]*(klist[2] + (dim-1)*klist[1]/q),
    k2norms[dim-1]*(k2list[2] + (dim-1)*k2list[1]/q),
    k3norms[dim-1]*(k3list[2] + (dim-1)*k3list[1]/q),
    m4norm[dim-1]*(m4list[2] + (dim-1)*m4list[1]/q),
    m4norm[dim-1]*(-2*m4list[1]/q),
    (q, 0.5, 2))
