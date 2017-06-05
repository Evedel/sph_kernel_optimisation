from sympy import *
import optbase as ob

q = symbols('q')

# wker = (4 * 2/pi**1/2 * sin((pi*q)/2)**4)/(pi*q)**4
# wlist, ok = ob.differentiatekernel(wker ,1)
# norms = ob.calculatenorms(wlist)
# print(wlist)
# print(norms)
wker = Piecewise((1, q < 1.33333333), (1.5 * (2 - q), q < 2.), (0, True))
wlist, ok = ob.integratekernel(wker,1)
plot(wlist[2], (q, 0, 2))
norms = ob.calculatenorms(wlist)
print(wlist)
print(norms)
