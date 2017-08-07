from sympy import *
import optbase as ob
import humpskernels as hk


q = symbols('q')
dim = 1
################################################################
# sinc4 = (sin(pi/2*q)**4)/(pi/2*q)**4
# s4list, ok = ob.differentiatekernel(sinc4, 1)
# wnorms = ob.calculatenorms(s4list[0])
# print(ob.printkernel(s4list, 2.0, wnorms, 'sinc4'))
# s4list[2] = wnorms[dim-1]*s4list[2]
# s4list[1] = wnorms[dim-1]*s4list[1]
# s4list[0] = wnorms[dim-1]*s4list[0]
# s4list.append(s4list[2] + (dim-1)*s4list[1]/q)
# s4list.append(-2*s4list[1]/q)
# ob.tabulatefunctions(s4list, 2.0, '.', 'sinc4.dat')
# #
# # # k3list, ok = hk.humpskerlnels(3., 1.)
# q2m4list, ok = hk.humpskerlnels(2., 1.)
# q2m4norm = ob.calculatenorms(q2m4list[0])
# # print(ob.printkernel(q2m4list, 2.0, q2m4norm, ' q2m4 '))
# q2m4list[2] = q2m4norm[dim-1]*q2m4list[2]
# q2m4list[1] = q2m4norm[dim-1]*q2m4list[1]
# q2m4list[0] = q2m4norm[dim-1]*q2m4list[0]
# q2m4list.append(q2m4list[2] + (dim-1)*q2m4list[1]/q)
# q2m4list.append(-2*q2m4list[1]/q)
# ob.tabulatefunctions(q2m4list, 2.0, '.', 'q2m4.dat')
# #
# #######################################################################
# m4, ok = hk.humpskerlnels(0., 1.)
# wlist, ok = ob.differentiatekernel(m4[2], 1)
# clist = ob.calculatenorms(wlist[0])
# print(clist)
# print(ob.printkernel(wlist, 2.0, clist, 'M4'))
# wlist[2] = clist[dim-1]*wlist[2]
# wlist[1] = clist[dim-1]*wlist[1]
# wlist[0] = clist[dim-1]*wlist[0]
# wlist.append(wlist[2] + (dim-1)*wlist[1]/q)
# wlist.append(-2*wlist[1]/q)
# ob.tabulatefunctions(wlist, 2.0, '.', 'm4.dat')
# #######################################################################
# m6 = Piecewise(\
#     ((3. - q)**5 - 6. * (2. - q)**5 + 15. * (1. - q)**5, q < 1.),\
#     ((3. - q)**5 - 6. * (2. - q)**5, q < 2.), \
#     ((3. - q)**5, q < 3.),\
#     (0, True))
# wlist, ok = ob.differentiatekernel(m6, 1)
# clist = ob.calculatenorms(wlist[0])
# # print(ob.printkernel(wlist, 2.0, clist, ' M6 '))
# wlist[2] = clist[dim-1]*wlist[2]
# wlist[1] = clist[dim-1]*wlist[1]
# wlist[0] = clist[dim-1]*wlist[0]
# wlist.append(wlist[2] + (dim-1)*wlist[1]/q)
# wlist.append(-2*wlist[1]/q)
# ob.tabulatefunctions(wlist, 3.0, '.', 'm6.dat')
# ########################################################################
zipm6 = Piecewise(\
        ((3. - q)**5 - 6. * (2. - q)**5 + 15. * (1. - q)**5, q < 1.),\
        ((3. - q)**5 - 6. * (2. - q)**5, q < 2.), \
        ((3. - q)**5, q < 3.),\
        (0, True)).subs(q, 1.5*q)
wlist, ok = ob.differentiatekernel(zipm6, 1)
clist = ob.calculatenorms(wlist[0])
# print(ob.printkernel(wlist, 2.0, clist, 'zipM6'))
wlist[2] = clist[dim-1]*wlist[2]
wlist[1] = clist[dim-1]*wlist[1]
wlist[0] = clist[dim-1]*wlist[0]
wlist.append(wlist[2] + (dim-1)*wlist[1]/q)
wlist.append(-2*wlist[1]/q)
ob.tabulatefunctions(wlist, 2.0, '.', 'zipm6.dat')
# #######################################################################
# gauss = exp(-q**2)
# wlist, ok = ob.differentiatekernel(gauss, 1)
# clist = ob.calculatenorms(wlist[0])
# print(ob.printkernel(wlist, 2.0, clist, 'Gauss'))
# wlist[2] = clist[dim-1]*wlist[2]
# wlist[1] = clist[dim-1]*wlist[1]
# wlist[0] = clist[dim-1]*wlist[0]
# wlist.append(wlist[2] + (dim-1)*wlist[1]/q)
# wlist.append(-2*wlist[1]/q)
# ob.tabulatefunctions(wlist, 3.0, '.', 'gauss.dat')
# ########################################################################
# zipm6optlap1d = Piecewise(\
#     ((3. - q)**5 -5.94431946621*(1.49377588576 - q)**5 + 14.9996141288*(1.02685715138 - q)**5, q < 1.02685715138),\
#     ((3. - q)**5 -5.94431946621*(1.49377588576 - q)**5, q < 1.49377588576),\
#     ((3. - q)**5, q < 3.),\
#     (0, True)).subs(q, 1.5*q)
# w0 = zipm6optlap1d
# name = ' zipM6 laplace optimised '
# wlist, ok = ob.differentiatekernel(w0, 1)
# clist = ob.calculatenorms(wlist[0])
# print(ob.printkernel(wlist, 2.0, clist, name))
# ob.tabulatekernel(wlist, 2.0, '.', 'mmq0zipm6.dat')

# mmc0fzipm6 = Piecewise(\
#         (74. - 69.*q - 4.5*q**2, q < 2./3.),\
#         (-6.*(2. - 1.5*q)**5 + (3. - 1.5*q)**5, q < 4./3.), \
#         ((3. - 1.5*q)**5, q < 2.),\
#         (0, True))
# mmc0fzipm6list, ok = ob.differentiatekernel(mmc0fzipm6, 1)
# mmc0fzipm6norm = ob.calculatenorms(mmc0fzipm6list[0])
# print(ob.printkernel(mmc0fzipm6list, 2.0, mmc0fzipm6norm, ' mmq0zipm6 '))
# ob.tabulatekernel(mmc0fzipm6list, 2.0, '.', 'mmq0zipm6.dat')

# mmc1fzipm6 = Piecewise(\
#         (135000000001./2025000000. + (-5400000004.*q - 275399999994.*q**2+178199999997.*q**3)/1800000000., q < 2./3.),\
#         (-6.*(2. - 1.5*q)**5 + (3. - 1.5*q)**5, q < 4./3.), \
#         ((3. - 1.5*q)**5, q < 2.),\
#         (0, True))
# mmc1fzipm6list, ok = ob.differentiatekernel(mmc1fzipm6, 1)
# mmc1fzipm6norm = ob.calculatenorms(mmc1fzipm6list[0])
# print(ob.printkernel(mmc1fzipm6list, 2.0, mmc1fzipm6norm, ' mmq1zipm6 '))
# ob.tabulatekernel(mmc1fzipm6list, 2.0, '.', 'mmq1zipm6.dat')

# mmc21fzipm6 = Piecewise(\
#         (61. + 65.*q - (765.*q**2)/2. + 405.*q**3 - (2295.*q**4)/16., q < 2./3.),\
#         (-6.*(2. - 1.5*q)**5 + (3. - 1.5*q)**5, q < 4./3.), \
#         ((3. - 1.5*q)**5, q < 2.),\
#         (0, True))
# mmc21fzipm6list, ok = ob.differentiatekernel(mmc21fzipm6, 1)
# mmc21fzipm6norm = ob.calculatenorms(mmc21fzipm6list[0])
# print(ob.printkernel(mmc21fzipm6list, 2.0, mmc21fzipm6norm, ' mmq21zipm6 '))
# ob.tabulatekernel(mmc21fzipm6list, 2.0, '.', 'mmq21zipm6.dat')

# mmc22fzipm6 = Piecewise(\
#         (65.5272 - 7.45802*q - 117.675*q**2 + 38.5*q**3 + 32.0833*q**4, q < 2./3.),\
#         (-6.*(2. - 1.5*q)**5 + (3. - 1.5*q)**5, q < 4./3.), \
#         ((3. - 1.5*q)**5, q < 2.),\
#         (0, True))
# mmc22fzipm6list, ok = ob.differentiatekernel(mmc22fzipm6, 1)
# mmc22fzipm6norm = ob.calculatenorms(mmc22fzipm6list[0])
# print(ob.printkernel(mmc22fzipm6list, 2.0, mmc22fzipm6norm, ' mmq22zipm6 '))
# ob.tabulatekernel(mmc22fzipm6list, 2.0, '.', 'mmq22zipm6.dat')

# gauss2r = [-0.00012340980408667956 + exp(-2.25*q**2),
#             0.001110688236780116 - (4.5*q)/exp(2.25*q**2),
#             -0.009440850012630987 - 4.5/exp(2.25*q**2) + (20.25*q**2)/exp(2.25*q**2)]
#
# gnorm2r = [0.8466567770110767,
#             0.7170821936438454,
#             0.6078976495961691]
#
# print(ob.printkernel(gauss2r, 2.0, gnorm2r, ' gauss2r '))
# ob.tabulatekernel(gauss2r, 2.0, '.', 'gauss2r.dat')

# mmq2m4ci = Piecewise( \
#     (-0.13572*(-1 + q) + 1./6.*q, q<1.),
#     (((2 - q)**3*q**2)/6., q < 2.),
#     (0, True))
# mmq2m4cilist, ok = ob.integratekernel(mmq2m4ci, 1)
# mmq2m4cinorm = ob.calculatenorms(mmq2m4cilist[0])
# print(ob.printkernel(mmq2m4cilist, 2.0, mmq2m4cinorm, ' mmq2m4ci '))
# ob.tabulatekernel(mmq2m4cilist, 2.0, '.', 'mmq2m4ci.dat')

# mmq2m4c2 = Piecewise( \
#     (1./6. - 0.0773694*(-1. + q)**2., q<1.),
#     (((2 - q)**3*q**2)/6., q < 2.),
#     (0, True))
# mmq2m4c2list, ok = ob.integratekernel(mmq2m4c2, 1)
# mmq2m4c2norm = ob.calculatenorms(mmq2m4c2list[0])
# print(ob.printkernel(mmq2m4c2list, 2.0, mmq2m4c2norm, ' mmq2m4c2 '))
# ob.tabulatekernel(mmq2m4c2list, 2.0, '.', 'mmq2m4c2.dat')

# mmq2m4c3 = Piecewise( \
#     (1./6. + 0.154751*(-1. + q)**3, q < 1.),
#     (((2 - q)**3*q**2)/6., q < 2.),
#     (0, True))
# mmq2m4c3list, ok = ob.integratekernel(mmq2m4c3, 1)
# mmq2m4c3norm = ob.calculatenorms(mmq2m4c3list[0])
# print(ob.printkernel(mmq2m4c3list, 2.0, mmq2m4c3norm, ' mmq2m4c3 '))
# ob.tabulatekernel(mmq2m4c3list, 2.0, '.', 'mmq2m4c3.dat')

# mmq2m4c4 = Piecewise( \
#     (1./6. - 0.270783*(-1. + q)**4, q < 1.),
#     (((2 - q)**3*q**2)/6., q < 2.),
#     (0, True))
# mmq2m4c4list, ok = ob.integratekernel(mmq2m4c4, 1)
# mmq2m4c4norm = ob.calculatenorms(mmq2m4c4list[0])
# print(ob.printkernel(mmq2m4c4list, 2.0, mmq2m4c4norm, ' mmq2m4c4 '))
# ob.tabulatekernel(mmq2m4c4list, 2.0, '.', 'mmq2m4c4.dat')

# w2ker = 88.964659192*exp(-3.04592817425*(q -0.914791070728)**2) + 7.01690442862 * (1. - q/2.) -0.99975783763
# mgaus, ok = ob.integratekernel(w2ker, 0)
# w2norms = ob.calculatenorms(mgaus[0])
#
# m4 = Piecewise(\
#         (0., q < 0.),
#         (0.25 * (2. - q)**3 - (1. - q)**3, q < 1.),\
#         (0.25 * (2. - q)**3, q < 2.),\
#         (0, True))
# m4list, ok = ob.differentiatekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])
#
# hkl, ok = hk.humpskerlnels(2., 1.)
# hkn = ob.calculatenorms(hkl[0])
# m4 = Piecewise(\
#         (1./6., q < 1.),\
#         (-1./6.*(-2. + q)**3*q*q, q < 2.),\
#         (0, True))
# m4list, ok = ob.integratekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])
# print(ob.printkernel(m4list, 2.0, hkn, ' 11M4 '))
# ob.tabulatekernel(m4list, 2.0, '.')
#
# m4 = Piecewise(\
#         (q/6., q < 1.),\
#         (-1./6.*(-2. + q)**3*q*q, q < 2.),\
#         (0, True))
# m4list, ok = ob.integratekernel(m4, 1)
# m4norm = ob.calculatenorms(m4list[0])
# print(ob.printkernel(m4list, 2.0, m4norm, ' -0M4 '))
# ob.tabulatekernel(m4list, 2.0, '.')

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
