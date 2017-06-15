from sympy import *
import os
from sympy.integrals.quadrature import gauss_legendre
import re
from sympy.parsing.sympy_parser import parse_expr as symparse
import subprocess as sp
import numpy as np
import time as t


Pwdpath  = os.getcwd()
q, A, B, C, D, E, F, G = symbols('q A B C D E F G')
Local_env = os.environ.copy()
Local_env["OMP_NUM_THREADS"] = "1"
BrokenPenalty = 100000

refmtsub = re.compile("((?!\*\*\d+)(?!\D\D\d+\.)\D\D\d+)|((!?\*\d+)\D\d+)|(/\d+)|((?!^\.\d+)^\d+)|((?!^-\d+\.)^-\d+)")
refmtser = re.compile("(\.\d\d\d\d\d+)")
research = [re.compile('Erf'),
            re.compile('Erfc'),
            re.compile('ArcTan'),
            re.compile('Log'),
            re.compile('Gamma'),
            re.compile('d0'),
            re.compile('Cos'),
            re.compile('Sin')]
reinsert = ['erf',
            'erfc',
            'atan',
            'log',
            'gamma',
            '',
            'cos',
            'sin']

def waitwithtimeout(process, timeout):
    t0 = t.time()
    tn = t0
    dt = 0.5
    while ((tn < t0 + timeout) and (process.poll() == None)):
        t.sleep(dt)
        tn += dt
    if (process.poll() == None):
        process.kill()

def solveproblem(klist, gen, popnum, gennum, makepath, tasktype, runnerpath, dimcase, estimatecolumn, aimvalue):
    w, dw, d2w = klist
    if (isinstance(d2w.args[0], float)):
        print(klist)
        print("Something wrong in solveproblem")
        exit(0)
    if not os.path.exists(Pwdpath + "/appsrc"):
        os.mkdir(Pwdpath + "/appsrc")
    os.system("cp -r " + makepath + "/* " + Pwdpath + "/appsrc")
    norm = calculatenorms(w)
    n2kernf90 = printkernel(klist, 2, norm, 'genesis')

    # print(n2kernf90)
    # exit(0)
    f = open(Pwdpath + "/appsrc/src/kernel/" + "n2ext.f90", 'w')
    f.write(n2kernf90)
    f.close()

    # rmoldexe = sp.Popen(["rm", Pwdpath + "/appsrc/execute"], stdout=sp.PIPE, cwd=Pwdpath+"/appsrc").wait()
    # tac = dt.now()
    # print(tac-tic)
    make = sp.Popen(["make", "-j", "debug=no"], stdout=sp.PIPE, stderr=sp.PIPE, cwd=Pwdpath+"/appsrc")
    waitwithtimeout(make, 30)
    # make.wait()
    # print()
    if make.poll() != 0:
        # print("Make failed: " + str(make.poll()))
        # print(make.stdout.read())
        # print(make.stderr.read())
        # exit(0)
        return evaluateresult("FAIL",0,0)
    # toe = dt.now()
    # print(toe-tac)
    gnstr = "%04d" %gennum
    ppstr = "%04d" %popnum
    workpath = Pwdpath + "/results/pp-" + ppstr + "/gn-" + gnstr
    if not os.path.exists(workpath):
        os.makedirs(workpath)

    if len(dimcase) == 1:
        if (dimcase[0] == 1):
            os.system("rm -rf " + workpath + "/*")
    else:
        os.system("rm -rf " + workpath + "/*")

    mvoldexe = sp.Popen(["mv", Pwdpath + "/appsrc/execute", workpath], stdout=sp.PIPE).wait()
    os.system("cp " + runnerpath + " " + workpath)
    kernfilename = "n2ext.f90"
    if len(dimcase) == 1:
        kernfilename = "n2ext-" + str(dimcase[0]) + "D.f90"
    f = open(workpath + "/" + kernfilename, 'w')
    f.write(n2kernf90)
    f.close()
    # Save kernel profile
    qt = np.linspace(0., 2., 100)
    kernprofilefilename = "kernprofile.dat"
    if len(dimcase) == 1:
        kernprofilefilename = "kernprofile-" + str(dimcase[0]) + "D.dat"
    f = open(workpath + "/"  + kernprofilefilename, 'w')
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
    processes = []
    for DE in dimcase:
        for TT in tasktype:
            processes.append(sp.Popen(["bash", "run-diff-err-h.sh", str(DE), TT], stdout=sp.PIPE, stderr=sp.PIPE, cwd=workpath, env=Local_env))
    for prc in processes:
        nameoferrfile = prc.communicate()[0].splitlines()[-1].decode().rstrip()
        lastcost = evaluateresult(workpath + "/" + nameoferrfile, estimatecolumn, aimvalue)
        cost += lastcost
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
                # quality = sqrt(sum((yi-aimy)**2 for yi in y))/len(y)
                quality = sum(abs(yi-aimy) for yi in y)
    return quality

def fixmathpy(s):
    for iro in range(0,len(research)):
        s = re.sub(research[iro], reinsert[iro], s)
    return s

def callmathdif(s):
    mathcmd = "wolframscript -noprompt"
    mathcmd = "echo \"InFun[q_] = " + s + "; IntRes[q_] = Simplify[D[InFun[q], q]]; \
                Print[IntRes[q]]; \
                Print[FortranForm[IntRes[q]]];\" | " + mathcmd

    mathrun = sp.Popen(["bash", "-c", mathcmd], stdout=sp.PIPE)
    mathrun.wait()
    MForm = mathrun.stdout.readline().decode().rstrip()
    FForm = mathrun.stdout.readline().decode().rstrip()
    return MForm, FForm

def linetruncaion(line):
    breakePos = 50
    newLine = ''
    if (len(line) > breakePos):
        spaces = [pos for pos, char in enumerate(line) if char == ' ']
        curentBreake = breakePos
        previousBreak = 0
        while( curentBreake < spaces[-1]):
            while (spaces[0] < curentBreake):
                del spaces[0]
            curentBreake += breakePos
            newLine += line[previousBreak:spaces[0] + 1] + ' &\n'
            previousBreak = spaces[0] + 1
        newLine += line[previousBreak:]
    else:
        newLine = line
    return newLine

def printkernel(klist,R,c,name):
    w, dw, d2w = klist
    e = symbols('e')
    nexp = N(exp(1))
    w = N(w.subs(e,nexp)).evalf(3)
    dw = N(dw.subs(e,nexp)).evalf(3)
    d2w = N(d2w.subs(e,nexp)).evalf(3)

    kernfile = ""
    kernfile += "module n2ext\n"
    kernfile += "  use const\n"
    kernfile += "  implicit none\n\n"
    kernfile += "  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase\n\n"
    kernfile += "  private\n\n"
    normstr = str(c[0]) + ", " + str(c[1]) + ", " + str(c[2])
    kernfile += "    real :: n2C(3) = (/ " + normstr + " /)\n"
    kernfile += "    character (len=10) :: kernelname = ' " + name + " '\n"
    kernfile += "    real :: krad = " + "{0:.1f}".format(R) + ", wCv\n"
    kernfile += "    integer :: dim\n"
    kernfile += "  contains\n\n"
    kernfile += "  subroutine setdimbase(d)\n"
    kernfile += "    integer, intent(in) :: d\n"
    kernfile += "    dim = d\n"
    kernfile += "    wCv = n2C(dim)\n"
    kernfile += "  end subroutine\n\n"
    kernfile += "  pure subroutine kf(q, f)\n"
    kernfile += "    real, intent(in)  :: q\n"
    kernfile += "    real, intent(out) :: f\n\n"
    # print(type(i2w))
    # exit(0)
    if isinstance(w, Piecewise):
        for i, (e, c) in enumerate(w.args):
            if i == 0:
                kernfile += "    if (q < 0.0) then\n"
                kernfile += "      if (isnan(q)) then\n"
                kernfile += "        error stop 'q is nan'\n"
                kernfile += "      else\n"
                kernfile += "        error stop 'q is negative'\n"
                kernfile += "      end if\n"
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            elif i == len(w.args)-1 and c == True:
                kernfile += "    else\n"
            else:
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            kernfile += "      f  = "+ linetruncaion(str(e)) + "\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n\n"
        kernfile += "  pure subroutine kdf(q, df)\n"
        kernfile += "    real, intent(in)  :: q\n"
        kernfile += "    real, intent(out) :: df\n\n"
        for i, (e, c) in enumerate(dw.args):
            if i == 0:
                kernfile += "    if (q < 0.0) then\n"
                kernfile += "      if (isnan(q)) then\n"
                kernfile += "        error stop 'q is nan'\n"
                kernfile += "      else\n"
                kernfile += "        error stop 'q is negative'\n"
                kernfile += "      end if\n"
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            elif i == len(dw.args)-1 and c == True:
                kernfile += "    else\n"
            else:
                kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
            kernfile += "      df  = " + linetruncaion(str(e)) + "\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n\n"
        kernfile += "  pure subroutine kddf(q, ddf)\n"
        kernfile += "    real, intent(in)  :: q\n"
        kernfile += "    real, intent(out) :: ddf\n\n"
        try:
            for i, (e, c) in enumerate(d2w.args):
                if i == 0:
                    kernfile += "    if (q < 0.0) then\n"
                    kernfile += "      if (isnan(q)) then\n"
                    kernfile += "        error stop 'q is nan'\n"
                    kernfile += "      else\n"
                    kernfile += "        error stop 'q is negative'\n"
                    kernfile += "      end if\n"
                    kernfile += "    elseif (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
                elif i == len(d2w.args)-1 and c == True:
                    kernfile += "    else\n"
                else:
                    kernfile += "    else if (" + fixmathpy(fcode(c).lstrip()) + ") then\n"
                kernfile += "      ddf  = " + linetruncaion(str(e)) + "\n"
        except:
            print("Exception lines 210-226:", d2w.args)
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n"
        kernfile += "end module n2ext\n"
    else:
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      if (isnan(q)) then\n"
        kernfile += "        error stop 'q is nan'\n"
        kernfile += "      else\n"
        kernfile += "        error stop 'q is negative'\n"
        kernfile += "      end if\n"
        kernfile += "    else if (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      f  = " + linetruncaion(str(w)) + "\n"
        kernfile += "    else\n"
        kernfile += "      f = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n\n"
        kernfile += "  pure subroutine kdf(q, df)\n"
        kernfile += "    real, intent(in)  :: q\n"
        kernfile += "    real, intent(out) :: df\n\n"
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      if (isnan(q)) then\n"
        kernfile += "        error stop 'q is nan'\n"
        kernfile += "      else\n"
        kernfile += "        error stop 'q is negative'\n"
        kernfile += "      end if\n"
        kernfile += "    else if (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      df  = " + linetruncaion(str(dw)) + "\n"
        kernfile += "    else\n"
        kernfile += "      df = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n\n"
        kernfile += "  pure subroutine kddf(q, ddf)\n"
        kernfile += "    real, intent(in)  :: q\n"
        kernfile += "    real, intent(out) :: ddf\n\n"
        kernfile += "    if (q < 0.0) then\n"
        kernfile += "      if (isnan(q)) then\n"
        kernfile += "        error stop 'q is nan'\n"
        kernfile += "      else\n"
        kernfile += "        error stop 'q is negative'\n"
        kernfile += "      end if\n"
        kernfile += "    elseif (q < " + "{0:.1f}".format(R) + ") then\n"
        kernfile += "      ddf  = " + linetruncaion(str(d2w)) + "\n"
        kernfile += "    else\n"
        kernfile += "      ddf = .0\n"
        kernfile += "    end if\n"
        kernfile += "  end subroutine\n"
        kernfile += "end module\n"
    return kernfile

def tabulatekernel(wlist, path):
        qt = np.linspace(0., 2., 100)
        kernprofilefilename = "kernprofile.dat"
        f = open(path + "/"  + kernprofilefilename, 'w')
        f.write("x3y {| x | w | dw | ddw |}\n")
        for i in range(len(qt)):
            f.write("%s %s %s %s\n" \
                %(qt[i], \
                N(wlist[0].subs(q,qt[i])), \
                N(wlist[1].subs(q,qt[i])), \
                N(wlist[2].subs(q,qt[i]))))

def callmathint(s):
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

def integratekernel(w, enginetype):
    ok = False
    # +-------------------------+
    # | Mathematica Integration |
    # +-------------------------+
    if enginetype == 0:
        WM = printing.mathematica.mathematica_code(w)
        IWM, IWF = callmathint(WM)
        I2WM, I2WF = callmathint(IWM)
        try:
            iw = symparse(fixmathpy(IWF))
            i2w = symparse(fixmathpy(I2WF))
            ok = True
        except:
            iw, i2w = "",""
            pass
    elif enginetype == 1:
    # +-------------------------+
    # |    Sympy Integration    |
    # +-------------------------+
        iw  = simplify(integrate(w,q))
        if (isinstance(iw, Piecewise)):
            iwtmp = list(iw.args)
            iwtmp[len(iw.args)-2] = (iwtmp[len(iw.args)-2][0] - iw.args[len(iw.args)-2][0].subs(q, 2.), iw.args[len(iw.args)-2][1])
            iw = Piecewise(*iwtmp)
            for i in range(0, len(iw.args)-2):
                iwtmp = list(iw.args)
                iwtmp[i] = (iw.args[i][0] - (iw.args[i][0].subs(q, 1.) - iw.args[i+1][0].subs(q, 1.)), iw.args[i][1])
                iw = Piecewise(*iwtmp)
        else:
            iw = iw - iw.subs(q,2.)

        i2w = simplify(integrate(iw,q))
        if (isinstance(i2w, Piecewise)):
            iwtmp = list(i2w.args)
            iwtmp[len(i2w.args)-2] = (iwtmp[len(i2w.args)-2][0] - i2w.args[len(i2w.args)-2][0].subs(q, 2.), i2w.args[len(i2w.args)-2][1])
            i2w = Piecewise(*iwtmp)
            for i in range(0, len(i2w.args)-2):
                iwtmp = list(i2w.args)
                iwtmp[i] = (i2w.args[i][0] - (i2w.args[i][0].subs(q, 1.) - i2w.args[i+1][0].subs(q, 1.)), i2w.args[i][1])
                i2w = Piecewise(*iwtmp)
        else:
            i2w = i2w - i2w.subs(q,2.)

        ok = True
    return [i2w, iw, w], ok

def differentiatekernel(w, enginetype=1):
    ok = False
    # +-----------------------------+
    # | Mathematica Differentiation |
    # +-----------------------------+
    if enginetype == 0:
        print("Not Implemented for This Engine")
        exit(127)
        WM = "Piecewise[{"
        if (isinstance(w,Piecewise)):
            for i, (e, c) in enumerate(w.args):
                if c != True:
                    WM += "{" + printing.mathematica.mathematica_code(e) + "," + str(c) + "},"
                else:
                    WM = WM[:-1]
                    WM += "},0]"
        else:
            print("Not Implemented for Regular Functions")
            exit(127)
        DWM, DWF = callmathdif(WM)
        D2WM, D2WF = callmathdif(DWM)
        try:
            dw = symparse(fixmathpy(DWF))
            d2w = symparse(fixmathpy(D2WF))
            ok = True
        except:
            dw, d2w = "",""
            pass
        # exit(0)
    elif enginetype == 1:
    # +-----------------------------+
    # |    Sympy Differentiation    |
    # +-----------------------------+
        dw  = simplify(diff(w,q))
        d2w = simplify(diff(dw,q))
        ok = True
    return [w, dw, d2w], ok

def calcdphidh(w):
    dphi = 4*pi/q**2*integrate(w*q**2,q)
    phi = integrate(dphi - dphi.subs(q,0.),q)
    dphidh = -phi - q * diff(phi,q)
    print(w)
    print(dphi)
    print(phi)
    print(dphidh)
    exit(0)

def callmathodesol(s, dim):
    # "y[q]/.DSolve[{y''[q] + y'[q] == Exp[-q], y[2] == 0, y'[2] == 0},y[q],q][[1]]"
    mathcmd = "wolframscript -noprompt"
    mathcmd = "echo \"InFun[q_] = " + s + ";\
                IntBCRes[q_] = Simplify[y[q]/.DSolve[{y''[q] + \
                " + str(dim - 1) + " * y'[q]/q == InFun[q], y[2] == 0, y'[2] == 0},y[q],q][[1]]]; \
                Print[IntBCRes[q]]; \
                Print[FortranForm[IntBCRes[q]]];\" | " + mathcmd
    # print(mathcmd)
    mathrun = sp.Popen(["bash", "-c", mathcmd], stdout=sp.PIPE)
    mathrun.wait()
    MForm = mathrun.stdout.readline().decode().rstrip()
    FForm = mathrun.stdout.readline().decode().rstrip()
    # print(MForm)
    return MForm, FForm

def callMathDiff(s):
    #
    # Chop will erase all the terms smaller then 10^(-14)
    mathcmd = "wolframscript -noprompt"
    mathcmd = "echo \"InF[q_] = " + s + "; \
                Res[q_] = Chop[Simplify[D[InF[q],q]],10^(-14)]; \
                Print[Res[q]]; \
                Print[FortranForm[Res[q]]];\" | " + mathcmd
    mathrun = sp.Popen(["bash", "-c", mathcmd], stdout=sp.PIPE)
    mathrun.wait()
    MForm = mathrun.stdout.readline().decode().rstrip()
    FForm = mathrun.stdout.readline().decode().rstrip()
    return MForm, FForm

def solvekernelode(f):
    ok = False
    WM = printing.mathematica.mathematica_code(f)
    mfsol1d, ffsol1d = callmathodesol(WM, 1)
    mfsol2d, ffsol2d = callmathodesol(WM, 2)
    mfsol3d, ffsol3d = callmathodesol(WM, 3)

    mfdsol1d, ffdsol1d = callMathDiff(mfsol1d)
    mfdsol2d, ffdsol2d = callMathDiff(mfsol2d)
    mfdsol3d, ffdsol3d = callMathDiff(mfsol3d)

    mfd2sol1d, ffd2sol1d = callMathDiff(mfdsol1d)
    mfd2sol2d, ffd2sol2d = callMathDiff(mfdsol2d)
    mfd2sol3d, ffd2sol3d = callMathDiff(mfdsol3d)

    try:
        w1d = symparse(fixmathpy(ffsol1d))
        w2d = symparse(fixmathpy(ffsol2d))
        w3d = symparse(fixmathpy(ffsol3d))

        dw1d = symparse(fixmathpy(ffdsol1d))
        dw2d = symparse(fixmathpy(ffdsol2d))
        dw3d = symparse(fixmathpy(ffdsol3d))

        d2w1d = symparse(fixmathpy(ffd2sol1d))
        d2w2d = symparse(fixmathpy(ffd2sol2d))
        d2w3d = symparse(fixmathpy(ffd2sol3d))

        ok = True
    except:
        w1d = ""
        w2d = ""
        w3d = ""
        dw1d = ""
        dw2d = ""
        dw3d = ""
        d2w1d = ""
        d2w2d = ""
        d2w3d = ""
        print("Convertation from Mathematica to SymPy failed")
        pass
    return [[w1d, w2d, w3d], [dw1d, dw2d, dw3d], [d2w1d, d2w2d, d2w3d]], ok

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

def calculatenorms(originkern):
    w = originkern
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

# f = Piecewise((0,q>2), (0.25*(2.0 - q)**3, q > 1), (0.25*(2.0 - q)**3 - (1.0 - q)**3, q>0), (0,True))
# print(differentiatekernel(f ,1))
# wlist, ok = differentiatekernel(f ,1)
# print(calculatenorms(wlist))
