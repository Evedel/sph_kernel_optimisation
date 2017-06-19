from sympy import *
import optbase as ob
import humpskernels as hk


q = symbols('q')
dim = 3

sinc4 = (sin(pi/2*q)**4)/(pi/2*q)**4
s4list, ok = ob.differentiatekernel(sinc4 ,1)
s4norms = ob.calculatenorms(s4list[0])
s4list[2] = s4norms[dim-1]*(s4list[2] + (dim-1)*s4list[1]/q)
s4list[1] = s4norms[dim-1]*s4list[1]
s4list[0] = s4norms[dim-1]*s4list[0]
ob.tabulatekernel(s4list, 2.0, '.', 'sinc4kp.dat')
# print(ob.printkernel(s4list, 2.0, wnorms, 'sinc4'))
# ob.calcdphidh(s4list[0])

# klist, ok = hk.humpskerlnels(2.11, 1.)
k2list, ok = hk.humpskerlnels(2., 1.)
# k3list, ok = hk.humpskerlnels(3., 1.)
# knorms = ob.calculatenorms(klist[0])
k2norms = ob.calculatenorms(k2list[0])

k2list[2] = k2norms[dim-1]*(k2list[2] + (dim-1)*k2list[1]/q)
k2list[1] = k2norms[dim-1]*k2list[1]
k2list[0] = k2norms[dim-1]*k2list[0]
ob.tabulatekernel(k2list, 2.0, '.', 'q2m4kp.dat')
# k3norms = ob.calculatenorms(k3list[0])

m4list, ok = hk.humpskerlnels(0., 1.)
m4list, ok = ob.differentiatekernel(m4list[2], 1)

m4norms = ob.calculatenorms(m4list[0])
m4list[2] = m4norms[dim-1]*(m4list[2] + (dim-1)*m4list[1]/q)
m4list[1] = m4norms[dim-1]*m4list[1]
m4list[0] = m4norms[dim-1]*m4list[0]
ob.tabulatekernel(m4list, 2.0, '.', 'm4kp.dat')

zipm6 = Piecewise(\
        ((3. - q)**5 - 6. * (2. - q)**5 + 15. * (1. - q)**5, q < 1.),\
        ((3. - q)**5 - 6. * (2. - q)**5, q < 2.), \
        ((3. - q)**5, q < 3.),\
        (0, True)).subs(q, 1.5*q)

zipm6list, ok = ob.differentiatekernel(zipm6, 1)
zipm6norm = ob.calculatenorms(zipm6list[0])
zipm6list[2] = zipm6norm[dim-1]*(zipm6list[2] + (dim-1)*zipm6list[1]/q)
zipm6list[1] = zipm6norm[dim-1]*zipm6list[1]
zipm6list[0] = zipm6norm[dim-1]*zipm6list[0]
ob.tabulatekernel(zipm6list, 2.0, '.', 'zipm6.dat')

# w2ker = 88.964659192*exp(-3.04592817425*(q -0.914791070728)**2) + 7.01690442862 * (1. - q/2.) -0.99975783763
# mgaus, ok = ob.integratekernel(w2ker, 0)
# w2norms = ob.calculatenorms(mgaus[0])

# m4 = Piecewise(\
#         (0., q < 0.),
#         (0.25 * (2. - q)**3 - (1. - q)**3, q < 1.),\
#         (0.25 * (2. - q)**3, q < 2.),\
#         (0, True))
# m4list, ok = ob.differentiatekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])

# m4 = Piecewise(\
#         (q/6., q < 1.),\
#         (-1./6.*(-2. + q)**3*q*q, q < 2.),\
#         (0, True))
# m4list, ok = ob.integratekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])
# print(ob.printkernel(m4list, 2.0, m4norm, ' -0M4 '))
# ob.tabulatekernel(m4list, 2.0, '.')
#
# m4 = Piecewise(\
#         (0.288888888888*(-1. + q) + q/6, q < 1.),\
#         (-1./6.*(-2. + q)**3*q*q, q < 2.),\
#         (0, True))
# m4list, ok = ob.integratekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])
# print(ob.printkernel(m4list, 2.0, m4norm, ' 101M4 '))
# ob.tabulatekernel(m4list, 2.0, '.')
#
# print(ob.printkernel(s4list, 2., norms, 'sinc4'))
# print(s4list)
# print(norms)
# wker = Piecewise((1, q < 1.33333333), (1.5 * (2 - q), q < 2.), (0, True))
# s4list, ok = ob.integratekernel(wker,1)
# plot(s4list[0], (q, 0, 2))
# plot(s4list[1], (q, 0, 2))

# print(m4list[2])
# dim = 1
#knorms[dim-1]*(klist[2] + (dim-1)*klist[1]/q),

# plot(wnorms[dim-1]*(s4list[2] + (dim-1)*s4list[1]/q),
#     w2norms[dim-1]*(mgaus[2] + (dim-1)*mgaus[1]/q),
#     k2norms[dim-1]*(k2list[2] + (dim-1)*k2list[1]/q),
#     k3norms[dim-1]*(k3list[2] + (dim-1)*k3list[1]/q),
#     m4norm[dim-1]*(m4list[2] + (dim-1)*m4list[1]/q),
#     m6norm[dim-1]*(m6list[2] + (dim-1)*m6list[1]/q),
#     m4norm[dim-1]*(-2*m4list[1]/q),
#     (q, 0.5, 3))
