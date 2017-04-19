from sympy import *
import os
from sympy.integrals.quadrature import gauss_legendre
import re
from sympy.parsing.sympy_parser import parse_expr as symparse
import subprocess as sp
import numpy as np



Pwdpath  = os.getcwd()
q, A, B, C, D, E, F, G = symbols('q A B C D E F G')
Local_env = os.environ.copy()
Local_env["OMP_NUM_THREADS"] = "1"
BrokenPenalty = 100000

refmtsub = re.compile("((?!\*\*\d+)(?!\D\D\d+\.)\D\D\d+)|((!?\*\d+)\D\d+)|(/\d+)|((?!^\.\d+)^\d+)|((?!^-\d+\.)^-\d+)")
refmtser = re.compile("(\.\d\d\d\d\d+)")
research = [re.compile('Erf'), re.compile('Erfc'), re.compile('ArcTan'), re.compile('Log'), re.compile('d0')]
reinsert = ['erf', 'erfc', 'atan', 'log', '']

def solveproblem(klist, gen, popnum, gennum, makepath, tasktype, runnerpath, dimcase, estimatecolumn, aimvalue):
    # tic = dt.now()
    if not os.path.exists(Pwdpath + "/appsrc"):
        os.mkdir(Pwdpath + "/appsrc")
    os.system("cp -r " + makepath + "/* " + Pwdpath + "/appsrc")
    norm = calculatenorms(klist)
    n2kernf90 = printkernel(klist, 2, norm, 'genesis')
    # print(n2kernf90)
    # exit(0)
    f = open(Pwdpath + "/appsrc/src/kernel/" + "n2ext.f90", 'w')
    f.write(n2kernf90)
    f.close()

    # rmoldexe = sp.Popen(["rm", Pwdpath + "/appsrc/execute"], stdout=sp.PIPE, cwd=Pwdpath+"/appsrc").wait()
    # tac = dt.now()
    # print(tac-tic)
    make = sp.Popen(["make", "-j"], stdout=sp.PIPE, stderr=sp.PIPE, cwd=Pwdpath+"/appsrc")
    make.wait()
    if make.poll() != 0:
        # print("Make failed: " + str(make.poll()))
        # print(make.stdout.read())
        # print(make.stderr.read())
        return evaluateresult("FAIL",0,0)
    # toe = dt.now()
    # print(toe-tac)
    gnstr = "%04d" %gennum
    ppstr = "%04d" %popnum
    workpath = Pwdpath + "/results/pp-" + ppstr + "/gn-" + gnstr
    if not os.path.exists(workpath):
        os.makedirs(workpath)

    os.system("rm -rf " + workpath + "/*")
    mvoldexe = sp.Popen(["mv", Pwdpath + "/appsrc/execute", workpath], stdout=sp.PIPE).wait()
    os.system("cp " + runnerpath + " " + workpath)
    f = open(workpath + "/" + "n2ext.f90", 'w')
    f.write(n2kernf90)
    f.close()
    # Save kernel profile
    qt = np.linspace(0., 2., 100)
    f = open(workpath + "/kernprofile.dat", 'w')
    genstr = "["
    for g in gen:
        genstr += str(g) + " "
    genstr += "]"
    f.write("x3y {| x | w | dw | ddw |} %s\n" %genstr)
    for i in range(len(qt)):
        f.write("%s %s %s %s\n" %(qt[i],N(klist[0].subs(q,qt[i])),N(klist[1].subs(q,qt[i])),N(klist[2].subs(q,qt[i]))))
    # f.close()
    # print(dt.now()-toe)
    cost = 0.
    for DE in dimcase:
        for TT in tasktype:
            scrrun = sp.Popen(["bash", "run-diff-err-h.sh", str(DE), TT], stdout=sp.PIPE, stderr=sp.PIPE, cwd=workpath, env=Local_env)
            scrrun.wait()
            nameoferrfile = scrrun.communicate()[0].splitlines()[-1].decode().rstrip()
            cost += evaluateresult(workpath + "/" + nameoferrfile, estimatecolumn, aimvalue)
    return cost

def evaluateresult(path, yidx, aimy):
    quality = BrokenPenalty
    if path != "FAIL":
        with open(path, 'r') as f:
            file_string = f.readline()
            y = []
            while file_string:
                file_string = f.readline()
                elements = file_string.split()
                if len(elements) != 0:
                    y.append(float(elements[yidx]))
            f.close()
            if (len(y) != 0):
                quality = sqrt(sum((i-aimy)*(i-aimy) for i in y)/len(y))
    return quality

def fixmathpy(s):
    for iro in range(0,len(research)):
        s = re.sub(research[iro], reinsert[iro], s)
    return s

def callmath(s):
    # mathcmd = "/Volumes/Mathematica/Mathematica.app/Contents/MacOS/WolframKernel -noprompt "
    mathcmd = "wolframscript -noprompt"
    mathcmd = "echo \"InFun[q_] = " + s + "; IntRes[q_] = Integrate[InFun[q], q]; \
                IntBCRes[q_] = Simplify[IntRes[q] - IntRes[2]]; Print[IntBCRes[q]]; \
                Print[FortranForm[IntBCRes[q]]];\" | " + mathcmd
    # print(mathcmd)
    mathrun = sp.Popen(["bash", "-c", mathcmd], stdout=sp.PIPE)
    mathrun.wait()
    MForm = mathrun.stdout.readline().decode().rstrip()
    FForm = mathrun.stdout.readline().decode().rstrip()
    return MForm, FForm

def printkernel(klist,R,c,name):
    w, dw, d2w = klist
    e = symbols('e')
    nexp = N(exp(1))
    w = N(w.subs(e,nexp))
    dw = N(dw.subs(e,nexp))
    d2w = N(d2w.subs(e,nexp))

    # print(w)
    # print(dw)
    # print(i2w)
    # exit(0)

    kernfile = ""
    kernfile += "module n2ext\n"
    kernfile += "  use const\n"
    kernfile += "  implicit none\n\n"
    kernfile += "  public :: n2f, n2df, n2ddf, n2C, n2R, n2Name\n\n"
    kernfile += "  private\n\n"
    normstr = str(c[0]) + ", " + str(c[1]) + ", " + str(c[2])
    kernfile += "    real :: n2C(3,2) = reshape((/ " + normstr + ",&\n"
    kernfile += "                                  " + normstr + " /),(/3,2/))\n"
    kernfile += "    character (len=10) :: n2Name=' " + name + " '\n"
    kernfile += "    real :: n2R = " + "{0:.1f}".format(R) + "\n"
    kernfile += " contains\n\n"
    kernfile += "  subroutine n2f(r, h, f)\n"
    kernfile += "    real, intent(in)  :: r, h\n"
    kernfile += "    real, intent(out) :: f\n"
    kernfile += "    real              :: q\n\n"
    kernfile += "    q = r / h\n"
    # print(type(i2w))
    # exit(0)
    if isinstance(w, Piecewise):
        for i, (e, c) in enumerate(w.args):
            if i == 0:
                kernfile += "    if (q < 0.0) then\n"
                kernfile += "      print *, 'f: something went wrong, q =', q\n"
                kernfile += "      stop\n"
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            elif i == len(w.args)-1 and c == True:
                kernfile += "    else\n"
            else:
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            kernfile += "      f  = "+ str(e) + "\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2f\n\n"
        kernfile += "  subroutine n2df(r, h, df)\n"
        kernfile += "    real, intent(in)  :: r, h\n"
        kernfile += "    real, intent(out) :: df\n"
        kernfile += "    real              :: q\n\n"
        kernfile += "    q = r / h\n"
        for i, (e, c) in enumerate(dw.args):
            e = e / q
            if i == 0:
                kernfile += "    if (q < 0.0) then\n"
                kernfile += "      print *, 'df: something went wrong, q =', q\n"
                kernfile += "      stop\n"
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            elif i == len(dw.args)-1 and c == True:
                kernfile += "    else\n"
            else:
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            kernfile += "      df  = " + str(e) + "\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2df\n\n"
        kernfile += "  subroutine n2ddf(r, h, ddf)\n"
        kernfile += "    real, intent(in)  :: r, h\n"
        kernfile += "    real, intent(out) :: ddf\n"
        kernfile += "    real              :: q\n\n"
        kernfile += "    q = r / h\n"
        for i, (e, c) in enumerate(d2w.args):
            if i == 0:
                kernfile += "    if (q < 0.0) then\n"
                kernfile += "      print *, 'ddf: something went wrong, q =', q\n"
                kernfile += "      stop\n"
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            elif i == len(d2w.args)-1 and c == True:
                kernfile += "    else\n"
            else:
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            kernfile += "      ddf  = " + str(e) + "\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2ddf\n"
        kernfile += "end module n2ext\n"
    else:
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      print *, 'f: something went wrong, q =', q\n"
        kernfile += "      stop\n"
        kernfile += "    elseif (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      f  = " + str(w) + "\n"
        kernfile += "    else\n"
        kernfile += "      f = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2f\n\n"
        kernfile += "  subroutine n2df(r, h, df)\n"
        kernfile += "    real, intent(in)  :: r, h\n"
        kernfile += "    real, intent(out) :: df\n"
        kernfile += "    real              :: q\n\n"
        kernfile += "    q = r / h\n"
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      print *, 'df: something went wrong, q =', q\n"
        kernfile += "      stop\n"
        kernfile += "    elseif (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      df  = " + str(dw) + "\n"
        kernfile += "    else\n"
        kernfile += "      df = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2df\n\n"
        kernfile += "  subroutine n2ddf(r, h, ddf)\n"
        kernfile += "    real, intent(in)  :: r, h\n"
        kernfile += "    real, intent(out) :: ddf\n"
        kernfile += "    real              :: q\n\n"
        kernfile += "    q = r / h\n"
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      print *, 'ddf: something went wrong, q =', q\n"
        kernfile += "      stop\n"
        kernfile += "    elseif (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      ddf  = " + str(d2w) + "\n"
        kernfile += "    else\n"
        kernfile += "      ddf = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine n2ddf\n"
        kernfile += "end module n2ext\n"
    return kernfile

def integratekernel(w, enginetype):
    ok = False
    # +-------------------------+
    # | Mathematica INtegration |
    # +-------------------------+
    if enginetype == 0:
        WM = printing.mathematica.mathematica_code(w)
        IWM, IWF = callmath(WM)
        I2WM, I2WF = callmath(IWM)
        try:
            iw = symparse(fixmathpy(IWF))
            i2w = symparse(fixmathpy(I2WF))
            ok = True
        except:
            iw, i2w = "",""
            pass
    elif enginetype == 1:
    # +-------------------------+
    # |    Sympy INtegration    |
    # +-------------------------+
        iw  = simplify(integrate(w,q))
        if (isinstance(iw, Piecewise)):
            iwtmp = list(iw.args)
            iwq2 = iw.args[len(iw.args)-2][0].subs(q, 2.)
            for i, (e, c) in enumerate(iw.args):
                if c != True:
                    iwtmp[i] = (e - iwq2, c)
            iw = Piecewise(*iwtmp)
        else:
            iw = iw - iw.subs(q,2.)

        i2w = simplify(integrate(iw,q))
        if (isinstance(i2w, Piecewise)):
            i2wtmp = list(i2w.args)
            i2wq2 = i2w.args[len(i2w.args)-2][0].subs(q, 2.)
            for i, (e, c) in enumerate(i2w.args):
                if c != True:
                    i2wtmp[i] = (e - i2wq2, c)
            i2w = Piecewise(*i2wtmp)
        else:
            i2w = i2w - i2w.subs(q,2.)
        ok = True
    # print(w)
    # print(iw)
    # print(i2w)
    return [i2w, iw, w], ok

def differentiatekernel(w, enginetype):
    ok = False
    # +-----------------------------+
    # |    Sympy Differentiation    |
    # +-----------------------------+
    dw  = simplify(diff(w,q))
    d2w = simplify(diff(dw,q))
    ok = True
    return [w, dw, d2w], ok

def stupidintegration(f,x,a,b,eps):
    dx = (b-a)/2.
    iprev = N(f.subs(x,dx) * dx)
    inext = iprev * 10
    while(abs(iprev - inext) > eps):
        iprev = inext
        inext = 0.
        dx = dx/2
        n  = int((b-a)/dx)
        for i in range(0,n+1):
            inext += N(f.subs(x, a + i*dx) * dx)
        print(n, inext)
    return inext

def quadintegration(f,x,a,b,n):
    t, w = gauss_legendre(n, 8)
    dt = a - 1
    inext = 0.
    for i in range(0,n):
        inext += N(f.subs(x,t[i] - dt) * w[i])
    iprev = inext * 10
    # while(abs(iprev - inext) > eps):
    #     iprev = inext
    #     inext = 0.
    #     n *= 2
    #     t, w = gauss_legendre(n, 8)
    #     for i in range(0,n):
    #         inext += N(f.subs(x,t[i] - dt) * w[i])
    #     print(n, inext)
    return inext

def calculatenorms(kernlist):
    w, dw, d2w = kernlist
    n = 20
    tg, wg = gauss_legendre(n, 8)
    dt = -1
    inext = 0.
    for i in range(0,n):
        inext += N(w.subs(q,tg[i] - dt) * wg[i])
    c1D = 1./(2*inext)
    inext = 0.
    for i in range(0,n):
        inext += N((2*pi*q*w).subs(q,tg[i] - dt) * wg[i])
    c2D = 1./inext
    inext = 0.
    for i in range(0,n):
        inext += N((4*pi*q*q*w).subs(q,tg[i] - dt) * wg[i])
    c3D = 1./inext
    return [c1D, c2D, c3D]
