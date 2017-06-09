from sympy import *
import optbase as ob

def humpskerlnels(inm1, inn1):
    # print("Aim form: q^{%d/%d} M4" %(m1,n1))

    q = symbols('q')
    n1 = symbols('n1')
    m1 = symbols('m1')

    m1=inm1
    n1=inn1


    f1 = 1/6*q**(m1/n1)*(4 - 6*q**2 + 3*q**3)
    f2 = -(1/6)*(-2 + q)**3*q**(m1/n1)

    if1 = (n1*(6*m1**2*n1*q**((m1 + n1)/n1)*(6 - 7*q**2 + 3*q**3) +             \
         m1**3*q**((m1 + n1)/n1)*(4 - 6*q**2 + 3*q**3) +                        \
         m1*n1**2*q**((m1 + n1)/n1)*(104 - 84*q**2 + 33*q**3) -                 \
         6*n1**3*(-4 + 2**(4 + m1/n1) + 8*q**(3 + m1/n1) - 3*q**(4 + m1/n1) -   \
            16*q**((m1 + n1)/n1))))/(6*(m1 + n1)*(m1 + 2*n1)*(m1 +           \
         3*n1)*(m1 + 4*n1))

    if2 = -((n1*(6*m1**2*n1*(-3 + q)*(-2 + q)**2*q**((m1 + n1)/n1) +           \
           m1**3*(-2 + q)**3*q**((m1 + n1)/n1) +                               \
           m1*n1**2*q**((m1 + n1)/n1)*(-208 + 228*q - 84*q**2 + 11*q**3) +      \
           6*n1**3*(2**(4 + m1/n1) + 24*q**(2 + m1/n1) - 8*q**(3 + m1/n1) +     \
              q**(4 + m1/n1) - 32*q**((m1 + n1)/n1))))/(6*(m1 + n1)*(m1 +     \
           2*n1)*(m1 + 3*n1)*(m1 + 4*n1)))

    iif1 = (n1**2*(6*m1**2*n1*q**(2 + m1/n1)*(8 - 8*q**2 + 3*q**3) +             \
         m1**3*q**(2 + m1/n1)*(4 - 6*q**2 + 3*q**3) -                           \
         m1*n1**2*(24 - 3*2**(6 + m1/n1) + 24*(-1 + 2**(2 + m1/n1))*q -        \
            188*q**(2 + m1/n1) + 102*q**(4 + m1/n1) - 33*q**(5 + m1/n1)) -     \
         6*n1**3*(4 - 2**(5 + m1/n1) + 20*(-1 + 2**(2 + m1/n1))*q -            \
            40*q**(2 + m1/n1) + 10*q**(4 + m1/n1) -                           \
            3*q**(5 + m1/n1))))/(6*(m1 + n1)*(m1 + 2*n1)*(m1 + 3*n1)*(m1 +   \
         4*n1)*(m1 + 5*n1))

    iif2 = -((n1**2*(6*m1**2*n1*(-4 + q)*(-2 + q)**2*q**(2 + m1/n1) +           \
           m1**3*(-2 + q)**3*q**(2 + m1/n1) +                                  \
           m1*n1**2*(-2 + q)*(3*2**(5 + m1/n1) + 188*q**(2 + m1/n1) -          \
              80*q**(3 + m1/n1) + 11*q**(4 + m1/n1)) +                        \
           6*n1**3*(-2**(5 + m1/n1) + 5*2**(4 + m1/n1)*q -                     \
              80*q**(2 + m1/n1) + 40*q**(3 + m1/n1) - 10*q**(4 + m1/n1) + q**(  \
              5 + m1/n1))))/(6*(m1 + n1)*(m1 + 2*n1)*(m1 + 3*n1)*(m1 +      \
           4*n1)*(m1 + 5*n1)))

    Y = Piecewise(\
                    (f1, q <= 1.),\
                    (f2, q <= 2.),\
                    (0, True))

    iY = Piecewise(\
                    (if1, q <= 1.),\
                    (if2, q <= 2.),\
                    (0, True))

    iiY = Piecewise(\
                    (iif1, q <= 1.),\
                    (iif2, q <= 2.),\
                    (0, True))

    # plot(Y,iY,iiY, (q,0.,2.))

    klist = [iiY, iY, Y]

    return klist, True

# print(iiY.evalf(3))
# print(iif1.subs(q,1.), '\n')
# print(iif2.subs(q,1.), '\n')
# print(iif2.subs(q,2.), '\n')
#
#
#
#
# klist, ok = ob.integratekernel(Y, 1)
# norms = ob.calculatenorms(klist[0])
# plot(klist[0], (q,0,2))
# plot(klist[1], (q,0,2))
# plot(klist[2], (q,0,2))
# print(ob.printkernel(humpskerlnel(2,2), 2., ob.calculatenorms(humpskerlnel(2,2)[0]), 'q*M4'))
