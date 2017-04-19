from sympy import *
from sympy.core.cache import *
import os
import random as rd
import matplotlib.pyplot as plt
from datetime import datetime as dt
import gc

q, A, B, C, D, E, F, G = symbols('q A B C D E F G')
Makepath = "/Users/sergeibiriukov/_git/moca_study/fortran/ND-assign"
Runnerpath = "/Users/sergeibiriukov/_git/moca_study/fortran/scripts/run-diff-err-h.sh"
Local_env = os.environ.copy()
Local_env["OMP_NUM_THREADS"] = "1"
Popsize = 20
Initgen1D = [1., 1., 1., 1., 1., 1./sqrt(pi)]
# Initgen2D = [1., 1., 1., 1., 1., 1./pi]
Initgen2D = [6.64475184836906, 1.26725568832901, -0.131900694871420, 1.85394903650122, -0.0746141118420310, -0.0937635009463731 + 1.0/pi]
Initgen3D = [1., 1., 1., 1., 1., 1./(pi)**(3./2.)]
Basefunction = A*exp(-B*(q + C)**2) + D * (1. - q/2.) + E
DimCase = 2
LuckySpicies = int(Popsize/4)
MutationRate = 25
MutationPower = 1       # Random(0,1) * MutPow
MutationWidth = 1       # Number of gens to mutate at ones
MathEngine = 0
MathEnStr = "Mathematica" if MathEngine==0 else "SymPy"
eps = 10e-12
TaskType = "diff-laplace"
EstimateColumn = 2
AimValue = 0.

def fmt(e):
    import re
    s = "%s" %e
    s = re.sub(refmtsub,"\g<0>.", s)
    f = sympify(s)
    g = "%s" %simplify(f)
    if re.search(refmtser,g):
       return s
    else:
       return g

def initfromseed(initlist, poplen):
    pop = []
    pop.append(initlist)
    for i in range(1,poplen):
        newpop = []
        for j in range(0,len(initlist)):
            newpop.append(initlist[j] + rd.uniform(-1, 5))
        pop.append(newpop)
    return pop

def indexfromchances(pbbties):
    rndval = rd.randint(0,100)
    if (rndval < pbbties[0]):
        return 0
    for i in range(1,len(pbbties)):
        if ((rndval > pbbties[i-1]) and (rndval <= pbbties[i])):
            return i
    # Don't return None
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
            print("Failed to find different parents")
        cross = rd.randint(1,len(genoms[0])-1)
        newgenoms.append(genoms[X][:cross] + genoms[Y][cross:])
        crossresults += "%d|%d|%d  " %(X, cross, Y)
    print("  Breeding Result: [%s]" %crossresults)
    return newgenoms

def mutation(genoms, luckynum, oafit, mnfit):
    currentMR = MutationRate
    currentMP = MutationPower * mnfit
    if (abs(oafit - mnfit) < 10e-12):
        print(" ----------------------------------------- ")
        print("<*>>>> The catastrophe has happened! <<<<*>")
        print(" ----------------------------------------- ")
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
                genoms[i][genidx] += rd.uniform(-1,1) * currentMP
    return genoms

def main():
    # gen = [1, 4.0, -0.8, -0.0773047, -0.00157556, 3./(sqrt(pi)) - 0.075, 30./(7.*pi) + 0.01, 50./(7.*pi**(3./2)) - 0.02]
    if DimCase == 1:
        Initgen = Initgen1D
    elif DimCase == 2:
        Initgen = Initgen2D
    elif DimCase == 3:
        Initgen = Initgen3D
    print("||====================----------------=================")
    print("|| Integration Engine: %s" %MathEnStr)
    print("|| OMP threads (Fort): %s" %Local_env["OMP_NUM_THREADS"])
    print("||      Base Function: %s" %Basefunction)
    print("||        Initial Gen: %s" %Initgen)
    print("||    Population Size: %s" %Popsize)
    print("||      Lucky Spicies: %s" %LuckySpicies)
    print("||      Mutation Rate: %s%%" %MutationRate)
    print("||     Mutation Power: %s" %MutationPower)
    print("||     Mutation Width: %s" %MutationWidth)
    print("||====================----------------=================")
    popnum = 0
    population = initfromseed(Initgen,Popsize)
    overall = 1000
    cost = [overall]
    while (overall/len(cost) > 10e-6):
        overall = 0
        overallinv = 0
        cost = []
        gennum = 0
        for genom in population:
            cA, cB, cC, cD, cE, ckd = genom
            w = Basefunction.subs([(A,cA), (B,cB), (C, cC), (D, cD), (E, cE)])
            tic = dt.now()
            klist, ok = integratekernel(w, MathEngine)
            tac = dt.now()
            if ok:
                cost.append(solveproblem(klist, ckd, genom, popnum, gennum, Makepath, [DimCase], [TaskType], EstimateColumn. AimValue))
            else:
                print("<*><*><*>\nIntegratin failed\nFunction: %s\nEngine: %s\n<*><*><*>\n" %(w,MathEnStr))
                cost.append(BrokenPenalty)
            toe = dt.now()
            print("pp-%04d/gn-%04d | Time: %s(int); %s(fit) | Error: %.8f" \
                    %(popnum, gennum, tac-tic, toe-tac, cost[-1]))
            gennum += 1
        # cost = [2.85852950692223, 3.27739835997030, 3.16406915975418, 2.77897308020322, 2.93216028306821]
        # cost = [2, 2, 2, 2, 2, 2, 2, 2, 2]
        tic = dt.now()
        for i in range(0,len(cost)):
            overall += cost[i]
            overallinv += 1/cost[i]
        chances = [(1/i)/(overallinv)*100 for i in cost]
        chancesstr = ""
        for i in chances:
            chancesstr += "%.2f  " %i
        print("Breedings Chances: %s" %(chancesstr))
        # print("  Avarage Fitness: %s" %(overall/len(cost)))
        # print("     Best Fitness: %s" %(min(cost)))
        population = crossover(population, chances, LuckySpicies)
        population = mutation(population, LuckySpicies, overall/len(cost), min(cost))
        popnum += 1
        tac = dt.now()
        print("--------------------------------[> Breeding Time %s, Best Fittness: %s" %(tac-tic, min(cost)))
        gc.collect()
# main()
