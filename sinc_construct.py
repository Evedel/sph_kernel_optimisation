from sympy import *
import optbase as ob

q = symbols('q')

wker = (4 * 2/pi**1/2 * sin((pi*q)/2)**4)/(pi*q)**4
wlist, ok = ob.differentiatekernel(wker ,1)
norms = ob.calculatenorms(wlist)
print(wlist)
print(norms)
