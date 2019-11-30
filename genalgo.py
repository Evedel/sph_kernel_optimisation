#!/usr/bin/env python

from sympy import *
from sympy.core.cache import *
import os
import random as rd
import matplotlib.pyplot as plt
from datetime import datetime as dt
import gc
import optbase as ob

q = symbols('q')
A = symbols('A0:10')

# Makepath = "/Users/sergeibiriukov/_git/moca_study/fortran/ND-assign"
Makepath = "/Users/sergeibiriukov/_git/moca_study/python/kernoptim/f90_src"
# Runnerpath = "/Users/sergeibiriukov/_git/moca_study/fortran/scripts/run-diff-err-h.sh"
Runnerpath = Makepath + "/run-diff-err-h.sh"
Local_env = os.environ.copy()
Local_env["OMP_NUM_THREADS"] = "1"
Popsize = 20
# Initgen1D = [1., 1., 1., 1., 1., 1./sqrt(pi)]
# Initgen2D = [6.64475184836906, 1.26725568832901, -0.131900694871420, 1.85394903650122, -0.0746141118420310, -0.0937635009463731 + 1.0/pi]
# Initgen3D = [1., 1., 1., 1., 1., 1./(pi)**(3./2.)]

# Basefunction = A[0]*exp(-A[1]*(q + A[2])**2) + A[3] * (1. - q/2.) + A[4]
# Initgen2D = [1., 1., 1., 1., 1.]
# krad = 2
# ftype = 0

# Basefunction = Piecewise(\
#                 (A[0]*(3. - q)**5 + A[1] * (2. - q)**5 + A[2] * (1. - q)**5, q < 1.),\
#                 (A[0]*(3. - q)**5 + A[1] * (2. - q)**5, q < 2.),\
#                 (A[0]*(3. - q)**5, q < 3.),\
#                 (0, True))
# krad = 3
# Initgen2D = [1., -6., 15.]
# ftype = 1

# Basefunction = Piecewise(\
#         (A[0]*(2. - q)**3 + A[1]*(1. - q)**3, q < 1.), \
#         (A[0]*(2. - q)**3, q < 2.),\
#         (0, True))
# krad = 2
# Initgen2D = [0.25, -1.]
# ftype = 1

Basefunction = Piecewise(\
                (A[0]*(3. - q)**5 + A[1] * (2. - q)**5 + A[2] * (1. - q)**5, q < 1.),\
                (A[0]*(3. - q)**5 + A[1] * (2. - q)**5, q < 2.),\
                (A[0]*(3. - q)**5, q < 3.),\
                (0, True)).subs(q, 1.5*q)
krad = 2
Initgen2D = [1., -6., 15.]
ftype = 1

DimCase = 2
LuckySpicies = int(Popsize/4)
MutationRate = 100
MutationPower = 1       # Random(0,1) * MutPow
MutationWidth = 1       # Number of gens to mutate at ones
MathEngine = 0 if (ftype==0) else 1
MathEnStr = "Mathematica" if MathEngine==0 else "SymPy"
eps = 10e-12
TaskType = "diff-laplace"
# TaskType = "diff-graddiv"
EstimateColumn = 2
AimValue = 0.

def pf(s):
    print(s, flush=True)

def fmt(e):
    import re
    s = "%s" %e
    s = re.sub(ob.refmtsub,"\g<0>.", s)
    f = sympify(s)
    g = "%s" %simplify(f)
    if re.search(ob.refmtser,g):
       return s
    else:
       return g

def initfromseed(initlist, poplen):
    pop = []
    pop.append(initlist)
    for i in range(1,poplen):
        newpop = []
        for j in range(0,len(initlist)):
            newpop.append(initlist[j] * rd.uniform(-2, 2))
        pop.append(newpop)
    return pop

def indexfromchances(pbbties):
    rndval = rd.randint(0,100)
    if (rndval < pbbties[0]):
        return 0
    for i in range(1,len(pbbties)):
        if ((rndval > pbbties[i-1]) and (rndval <= pbbties[i])):
            return i
    return 0

def isgeninlist(genlist, newgen):
    lengenlist = len(genlist)
    if lengenlist == 0:
        return False
    else:
        for i in range(lengenlist):
            listdiff = 0
            for k in range(len(genlist[i])):
                listdiff += abs(genlist[i][k] - newgen[k])
            if listdiff < eps:
                return True
        return False

def crossover(genoms, chances, luckynum):
    newgenoms = []
    gc = genoms[:]
    cc = chances[:]

    iterwall = Popsize - 1
    while ((len(newgenoms) < luckynum) and (iterwall > 0)):
        pctluck = max(cc)
        idxluck = cc.index(pctluck)
        genluck = gc[idxluck]
        if not isgeninlist(newgenoms, genluck):
            newgenoms.append(genluck)
        cc.remove(pctluck)
        gc.remove(genluck)
        iterwall -= 1

    pbbties = [chances[0]]
    for i in range(1, Popsize-1):
        pbbties.append(chances[i] + pbbties[i-1])
    pbbties.append(100)

    crossresults = ""
    trialsnum = 25
    for i in range(0, Popsize - len(newgenoms)):
        k = 0
        X, Y = 0, 0
        while ((X == Y) and (k < trialsnum)):
            X = indexfromchances(pbbties)
            Y = indexfromchances(pbbties)
            k += 1
        if (k > trialsnum):
            pf("Failed to find different parents")
        cross = rd.randint(1,len(genoms[0])-1)
        newgenoms.append(genoms[X][:cross] + genoms[Y][cross:])
        crossresults += "%d|%d|%d  " %(X, cross, Y)
    pf("  Breeding Result: [%s]" %crossresults)
    return newgenoms

def mutation(genoms, luckynum, oafit, mnfit):
    currentMR = MutationRate
    currentMP = MutationPower * mnfit
    if (abs(oafit - mnfit) < 10e-12):
        pf(" ----------------------------------------- ")
        pf("<*>>>> The catastrophe has happened! <<<<*>")
        pf(" ----------------------------------------- ")
        currentMR = 100
        # To make it proportional to current error.
        #  Bigger error  => bigger power.
        #  Smaller error => more accurate power
        #  Could force to stay in local minima?
        currentMP = 5 * mnfit
    for i in range(luckynum,len(genoms)):
        if (rd.randint(0,100) < currentMR):
            for j in range(0,MutationWidth):
                genidx = rd.randint(0,len(genoms[0])-1)
                genoms[i][genidx] += rd.uniform(-1,1) * currentMP * abs(genoms[i][genidx])
    return genoms

def main():
    popnum = 0
    if DimCase == 1:
        Initgen = Initgen1D
    elif DimCase == 2:
        Initgen = Initgen2D
    elif DimCase == 3:
        Initgen = Initgen3D
    pf("||====================----------------=================")
    pf("|| Integration Engine: %s" %MathEnStr)
    pf("|| OMP threads (Fort): %s" %Local_env["OMP_NUM_THREADS"])
    pf("||      Base Function: %s" %Basefunction)
    pf("||        Initial Gen: %s" %Initgen)
    pf("||    Population Size: %s" %Popsize)
    pf("||      Lucky Spicies: %s" %LuckySpicies)
    pf("||      Mutation Rate: %s%%" %MutationRate)
    pf("||     Mutation Power: %s" %MutationPower)
    pf("||     Mutation Width: %s" %MutationWidth)
    pf("||====================----------------=================")

    # Initgen = [-49.90003995314704, 2.267079626421994, 9.108311456646012, 1.0, 1.0]
    # popnum = 30

    population = initfromseed(Initgen,Popsize)
    overall = 1000
    cost = [overall]
    while (overall/len(cost) > 10e-6):
        overall = 0
        overallinv = 0
        cost = []
        gennum = 0
        for genom in population:

            w = Basefunction
            for xi in range(len(genom)):
                w = w.subs(A[xi],genom[xi])
            if (ftype==1):
                w = simplify(w - w.subs(q,krad))
            tic = dt.now()
            if (ftype==0):
                klist, ok = ob.integratekernel(w, MathEngine)
            if (ftype==1):
                klist, ok = ob.differentiatekernel(w, MathEngine)
            tac = dt.now()
            if ok:
                cost.append(
                    ob.solveproblem(
                        krad,
                        klist,
                        genom,
                        popnum,
                        gennum,
                        Makepath,
                        [TaskType],
                        Runnerpath,
                        [DimCase],
                        EstimateColumn,
                        AimValue))
            else:
                if (ftype==0):
                    pf("<*><*><*>\nIntegratin failed\nFunction: %s\nEngine: %s\n<*><*><*>\n" %(w,MathEnStr))
                if (ftype==1):
                    pf("<*><*><*>\nDifferentiation failed\nFunction: %s\nEngine: %s\n<*><*><*>\n" %(w,MathEnStr))
                cost.append(ob.BrokenPenalty)
            toe = dt.now()
            pf("pp-%04d/gn-%04d | Time: %s(int); %s(fit) | Error: %.8f" \
                    %(popnum, gennum, tac-tic, toe-tac, cost[-1]))
            gennum += 1
        tic = dt.now()
        for i in range(0,len(cost)):
            overall += cost[i]
            overallinv += 1/cost[i]
        chances = [(1/i)/(overallinv)*100 for i in cost]
        chancesstr = ""
        for i in chances:
            chancesstr += "%.2f  " %i
        pf("Breedings Chances: %s" %(chancesstr))
        # pf("  Avarage Fitness: %s" %(overall/len(cost)))
        # pf("     Best Fitness: %s" %(min(cost)))
        population = crossover(population, chances, LuckySpicies)
        population = mutation(population, LuckySpicies, overall/len(cost), min(cost))
        popnum += 1
        tac = dt.now()
        pf("--------------------------------[> Breeding Time %s, Best Fittness: %s" %(tac-tic, min(cost)))
        gc.collect()

pf("Started")
main()
pf("Finished")
