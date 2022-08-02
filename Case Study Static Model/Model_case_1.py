""" Model_case_1

The model for the case study.

This file contains the following components:
- Exact model formulation, without any carbon regulations (except pre-implementation CCS);
- Carbon flow and carbon intensity constraints;
- Two ways to calculate carbon-intensities post-optimization.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

""" Import packages & data """
import gurobipy as gp
from gurobipy import GRB
import time
from itertools import product

# Import the data and the parameters
from Data_case_1 import *

############################## Set up model & decision variables ######################################################

""" Initiate model """
Model = gp.Model()
# Model.params.nonConvex = 2
Model.params.timelimit = timeLimit
if not printStatus:
    Model.Params.LogToConsole = 0

""" Decision variables """

# Demand locally satisfied
DL_ig = Model.addMVar(shape=(len(I), len(G)), lb=0.0)

# Demand imported
DI_ig = Model.addMVar(shape=(len(I), len(G)), lb=0.0)

# Amount produced
P_pig = Model.addMVar(shape=(len(P), len(I), len(G)), lb=0.0)

# Total production
PT_ig = Model.addMVar(shape=(len(I), len(G)), lb=0.0)

# Flows between grids
Q_ilgg = Model.addMVar(shape=(len(I), len(L), len(G), len(G)), lb=0.0)

# Cost variables
FC = Model.addVar(lb=0.0)
FCC = Model.addVar(lb=0.0)
FOC = Model.addVar(lb=0.0)
FSC = Model.addVar(lb=0.0)
GC = Model.addVar(lb=0.0)
LC = Model.addVar(lb=0.0)
MC = Model.addVar(lb=0.0)
TCC = Model.addVar(lb=0.0)
TCE = Model.addVar(lb=0.0)
TCEFeed = Model.addVar(lb=0.0)
TCEProd = Model.addVar(lb=0.0)
TCETrans = Model.addVar(lb=0.0)
TDC = Model.addVar(lb=0.0)
TOC = Model.addVar(lb=0.0)

# Number of plants
NP_pig = Model.addMVar(shape=(len(P), len(I), len(G)), vtype=GRB.INTEGER)

# Number of transport units from grid to grid
NTU_ilgg = Model.addMVar(shape=(len(I), len(L), len(G), len(G)), vtype=GRB.INTEGER)
# Number of transport units total
NTU_il = Model.addMVar(shape=(len(I), len(L)), vtype=GRB.INTEGER)

# Binary for some transport
X_ilgg = Model.addMVar(shape=(len(I), len(L), len(G), len(G)), vtype=GRB.BINARY)

# Binary for export
Y_ig = Model.addMVar(shape=(len(I), len(G)), vtype=GRB.BINARY)

# Binary for import
Z_ig = Model.addMVar(shape=(len(I), len(G)), vtype=GRB.BINARY)

############################## Model objective & constraints ##########################################################

""" Set some variables to 0 (ig, il and gg options) """

# No local production if no feasible location (ig)
Model.addConstrs(DL_ig[(i, g)] == 0 for i in I for g in G if A_ig[(i, g)] == 0)

# No production by a plant if no feasible location (ig)
Model.addConstrs(P_pig[(p, i, g)] == 0 for p in P for i in I for g in G if A_ig[(i, g)] == 0)

# No total production in a grid plant ino feasible location (ig)
Model.addConstrs(PT_ig[(i, g)] == 0 for i in I for g in G if A_ig[(i, g)] == 0)

# No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)  # TODO  or A_gg[(g, g_prime)] == 0)

# No plants can be build on locations where plants are not feasible (ig)
Model.addConstrs(NP_pig[(p, i, g)] == 0 for p in P for i in I for g in G if A_ig[(i, g)] == 0)

# No transport units if product-mode combination does not exist (il)
Model.addConstrs(NTU_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0)

# No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(X_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)  # TODO or A_gg[(g, g_prime)] == 0)

# No export if you cannot produce the product at the location (ig)
Model.addConstrs(Y_ig[(i, g)] == 0 for i in I for l in L for g in G for g_prime in G if A_ig[(i, g)] == 0)

""" Operational part """

# (B.1) Not more demand satisfied locally than totally produced in that grid
Model.addConstrs(DL_ig[(i, g)] <= PT_ig[(i, g)] for i in I for g in G if A_ig[(i, g)] == 1)

# (B.2) NEW - local demand satisfied are the "diagonals" on the Q matrix
Model.addConstrs(gp.quicksum(Q_ilgg[(i, l, g, g)] for l in L if A_il[(i, l)] == 1) == DL_ig[(i, g)] for i in I for g in G)

# (B.3) Imported hydrogen (has NO constraint on location)
Model.addConstrs(DI_ig[(i, g)] == gp.quicksum(Q_ilgg[(i, l, g_prime, g)] for l in L for g_prime in G if g_prime != g) for i in I for g in G)

# (B.4) Total demand satisfied
Model.addConstrs(DT_g[g] == gp.quicksum(DL_ig[(i, g)] + DI_ig[(i, g)] for i in I) for g in G)

# (B.5) Balance constraints
Model.addConstrs(DL_ig[(i, g)] + DI_ig[(i, g)]
                 == PT_ig[(i, g)] + gp.quicksum(Q_ilgg[(i, l, g_prime, g)] for l in L for g_prime in G if g_prime != g)
                 - gp.quicksum(Q_ilgg[(i, l, g, g_prime)] for l in L for g_prime in G if g_prime != g) for i in I for g in G)

# (B.6)
Model.addConstrs(PT_ig[(i, g)] == gp.quicksum(P_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)

# (B.7)
Model.addConstrs(P_pig[(p, i, g)] >= PCapMin_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G if A_ig[(i, g)] == 1)
Model.addConstrs(P_pig[(p, i, g)] <= PCapMax_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G if A_ig[(i, g)] == 1)

# (B.8)
Model.addConstrs(PT_ig[(i, g)] >= gp.quicksum(PCapMin_pi[(p, i)] * NP_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)
Model.addConstrs(PT_ig[(i, g)] <= gp.quicksum(PCapMax_pi[(p, i)] * NP_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)

# (B.9)
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] >= QMin_il[(i, l)] * X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] <= QMax_il[(i, l)] * X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)

# (B.10)
Model.addConstrs(X_ilgg[(i, l, g, g_prime)] + X_ilgg[(i, l, g_prime, g)] <= 1 for i in I for l in L for g in G for g_prime in G if g != g_prime)

# (B.11)
Model.addConstrs(Y_ig[(i, g)] >= X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G if g != g_prime if A_ig[(i, g)] == 1)

# (B.12)
Model.addConstrs(Z_ig[(i, g)] >= X_ilgg[(i, l, g_prime, g)] for i in I for l in L for g_prime in G for g in G if g != g_prime)

# (B.13)
Model.addConstrs(Y_ig[(i, g)] + Z_ig[(i, g)] <= 1 for i in I for g in G)

if onlyOneImportProduct:
    # (B.14) NEW - cannot import both products (does not entirely fix issue!)
    Model.addConstrs(gp.quicksum(Z_ig[(i, g)] for i in I) <= 1 for g in G)

""" Cost part """

# (B.19)
Model.addConstr(FCC == gp.quicksum(PCC_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G))

# (B.20)
Model.addConstrs(NTU_ilgg[(i, l, g, g_prime)] >= Q_ilgg[(i, l, g, g_prime)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l])
                                                 for g in G for g_prime in G for i in I for l in L if g_prime != g if A_il[(i, l)] == 1)

# (B.21)
Model.addConstrs(NTU_ilgg[(i, l, g, g)] >= Q_ilgg[(i, l, g, g)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l])
                                           for i in I for l in L for g in G if A_il[(i, l)] == 1)

# (B.22)
Model.addConstrs(NTU_il[(i, l)] == gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] for g in G for g_prime in G) for i in I for l in L if A_il[(i, l)] == 1)

# (B.23)
Model.addConstr(TCC == gp.quicksum(TMC_il[(i, l)] * NTU_il[(i, l)] for i in I for l in L))

# (B.24)
Model.addConstr(FOC == gp.quicksum(UPC_pi[(p, i)] * P_pig[(p, i, g)] for p in P for i in I for g in G))

# (B.25)
Model.addConstr(FC == gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilgg[(i, l, g, g_prime)] * (1 / FEbetween_l[l]) * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                    + gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g)] * Q_ilgg[(i, l, g, g)] * (1 / FEwithin_l[l]) * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G if A_il[(i, l)] == 1) )

# (B.26)
Model.addConstr(LC == gp.quicksum(DW_l[l] * (Q_ilgg[(i, l, g, g_prime)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                    + gp.quicksum(DW_l[l] * (Q_ilgg[(i, l, g, g)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l]))
                                  for i in I for l in L for g in G if A_il[(i, l)] == 1) )

# (B.27)
Model.addConstr(MC == gp.quicksum(ME_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilgg[(i, l, g, g_prime)] * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1))

# (B.28)
Model.addConstr(GC == gp.quicksum(GE_l[l] * NTU_il[(i, l)] for i in I for l in L))

# (B.29) (new)
Model.addConstr(FSC == gp.quicksum(FSP_p[p] * (1/PEF_p[p]) * gp.quicksum(P_pig[(p, i, g)] for i in I for g in G) for p in P))

# (B.30)
Model.addConstr(TOC == FC + LC + MC + GC)

# (B.231)
Model.addConstr(TDC == (1 / (alpha * CCF)) * (FCC + TCC) + FOC + TOC + FSC)

# (B.32)
objective = costImportance * TDC + CO2Importance * TCE


""" CO2 part """

# (B.33) Total feed emissions
Model.addConstr(TCEFeed == gp.quicksum(CO2Feed_p[p] * (gp.quicksum(P_pig[(p, i, g)] for g in G for i in I)) for p in P))

# (B.34) Total production emissions
Model.addConstr(TCEProd == gp.quicksum(CO2Prod_pi[(p, i)] * (gp.quicksum(P_pig[(p, i, g)] for g in G)) for p in P for i in I))

# (B.35) Total transportation emissions
Model.addConstr(TCETrans == CO2Trans * gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] * 2 * L_lgg[(l, g, g_prime)] for i in I for l in L for g in G for g_prime in G))
# Model.addConstr(TCETrans == gp.quicksum(CO2Trans * (gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] * 2 * L_lgg[(l, g, g_prime)] for i in I for l in L)) for g in G for g_prime in G))

# (B.36) Total chain emissions
Model.addConstr(TCE == TCEFeed + TCEProd + TCETrans)

# (B.37) Linear constraint on total chain emissions
if TotalCO2Max != -1:
    Model.addConstr(TCE <= TotalCO2Max)

############################## Execute ##################################################################

""" Execute the model as specified in Params. """
start_time = time.time()
Model.setObjective(objective, GRB.MINIMIZE)
Model.optimize()
end_time = time.time()
total_time = end_time - start_time

########################## Carbon intensities post-optimization ##########################################

""" Options: """
# Carbon_ig : total daily CO2 emissions for product i in grid g (per day)
# Carbon_g  : total daily CO2 emissions for grid g (per day)
# CI_ig     : carbon intensity of product i in grid g (at a day)
# CI_g      : avg. carbon intensity of utilized products in grid g

# (B.38) Calculate the starting carbon intensity as parameter per plant p per product i
CI_pi = np.zeros(shape=(len(P), len(I)))  # ton CO2/ton H2
for p in P:
    for i in I:
        CI_pi[(p, i)] = CO2Feed_p[p] + CO2Prod_pi[(p, i)]

# (B.39) Calculate the carbon intensities of the all the plants producing product i in plant g
# Note that the simplifying assumption here is that what is produced totally in a grid is distributed equally over the other grids.
#  in other words, if you have one CG plant in a grid and one SMR plant in the same grid, you average the emissions of those grids for the
#  outgoing products, while in fact you could connect some cities with the CG plant and some cities with the SMR plant
CIPlants_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        if PT_ig[(i, g)].X > 0.0000001:
            CIPlants_ig[(i, g)] = gp.quicksum(P_pig[(p, i, g)]/PT_ig[(i, g)].X * CI_pi[(p, i)] for p in P).getValue()

# (B.40) carbon intensities from production units for product i at customer location g
CIProd_ig = np.zeros(shape=(len(I), len(G)))
CIProd_denominator_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        CIProd_denominator_ig [(i, g)] = gp.quicksum(Q_ilgg[(i, l_star, g_prime_star, g)] for l_star in L for g_prime_star in G).getValue()
for i in I:
    for g in G:
        if CIProd_denominator_ig[(i, g)] > 0.0000001:
            CIProd_ig[(i, g)] = gp.quicksum(Q_ilgg[(i, l, g_prime, g)]/CIProd_denominator_ig[(i, g)] * CIPlants_ig[(i, g_prime)] for l in L for g_prime in G).getValue()

CITrans_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        if (DL_ig[(i, g)].X + DI_ig[(i, g)].X) > 0.0000001:
            CITrans_ig[(i, g)] = ( (CO2Trans * gp.quicksum(NTU_ilgg[(i, l, g_prime, g)] * 2 * L_lgg[(l, g_prime, g)] for l in L for g_prime in G))/(DL_ig[(i, g)].X + DI_ig[(i, g)].X) ).getValue()

CI_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        CI_ig[(i, g)] = CIProd_ig[(i, g)] + CITrans_ig[(i, g)]

CI_g = np.zeros(shape=len(G))
for g in G:
    CI_g[g] = gp.quicksum(CI_ig[(i, g)] * (DL_ig[(i, g)] + DI_ig[(i, g)])/DT_g[g] for i in I).getValue()

Carbon_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        Carbon_ig[(i, g)] = CI_ig[(i, g)] * (DL_ig[(i, g)].X + DI_ig[(i, g)].X)

Carbon_g = np.zeros(shape=(len(G)))
for g in G:
    Carbon_g[g] = CI_g[g] * DT_g[g]

############################## Print the results #########################################################

if printVariables:
    print("===== VARIABLES =====")
    print("Num. variables   : ", Model.numVars)  # (including those that are by definition 0)
    print("Num. constraints : ", Model.numConstrs)  # (including those that make variables 0 and redundant constraints for 0 variables)
    print()

    print("=== OPERATIONAL VARIABLES ===")
    print()

    # Continuous variables
    print("D^L_ig (local demand satisfied): \n", pd.DataFrame(DL_ig.X.round(2), index=ISet, columns=GSet).to_string())
    print()

    print("D^I_ig (imported demand satisfied): \n", pd.DataFrame(DI_ig.X.round(2), index=ISet, columns=GSet).to_string())
    print()

    # Intermezzo: Integer variables
    for p in P:
        print(f"NP_[{PSet[p]}]ig (number of plant {PSet[p]}): \n {pd.DataFrame(NP_pig[p, :, :].X.round(2), index=ISet, columns=GSet).to_string()}")
        print()

    for i in I:
        for l in L:
            if A_il[(i, l)] == 1:
                print(f"NTU_ilgg (transport product {ISet[i]} using mode {LSet[l]}): \n {pd.DataFrame(NTU_ilgg[i, l, :, :].X.round(2), index=GSet, columns=GSet).to_string()}")
                print()

    print("TOTAL NTU_il (number of units): \n", pd.DataFrame(NTU_il.X.round(2), index=ISet, columns=LSet).to_string())
    print()

    # Continuous variables:
    for p in P:
        print(f"P_[{PSet[p]}]ig (production from plant {PSet[p]}): \n {pd.DataFrame(P_pig[p, :, :].X.round(2), index=ISet, columns=GSet).to_string()}")
        print()

    print("P^T_ig (total production): \n", pd.DataFrame(PT_ig.X.round(2), index=ISet, columns=GSet).to_string())
    print()

    for i in I:
        for l in L:
            if A_il[(i, l)] == 1:
                print(f"Q_[{ISet[i]}][{LSet[l]}]gg' (transport product {ISet[i]} using mode {LSet[l]}): \n {pd.DataFrame(Q_ilgg[i, l, :, :].X.round(2), index=GSet, columns=GSet).to_string()}")
                print()

    # Total feedstock needed:
    FeedNeeded_p = np.zeros(shape=len(P))
    for p in P:
        FeedNeeded_p[p] = (PEF_p[p] * gp.quicksum(P_pig[(p, i, g)] for i in I for g in G)).getValue()
    print("FeedNeeded_p \n", pd.DataFrame(FeedNeeded_p, index=PSet).to_string())
    print()

    # Binary variables
    for i in I:
        for l in L:
            if A_il[(i, l)] == 1:
                print(f"X_[{ISet[i]}][{LSet[l]}]gg' (binary for transport product {ISet[i]} using mode {LSet[l]}): \n {pd.DataFrame(X_ilgg[i, l, :, :].X, index=GSet, columns=GSet).to_string()}")
                print()

    print("Y_ig (binary export): \n", pd.DataFrame(Y_ig.X, index=ISet, columns=GSet).to_string())
    print()

    print("Z_ig (binary import): \n", pd.DataFrame(Z_ig.X, index=ISet, columns=GSet).to_string())
    print()

print()
print("===== RESULT =====")
print(f"Objective  : {round(objective.getValue(), 3)}")
print(f"Total time : {round(total_time, 3)} sec")
print()

print("=== COST VARIABLES ===")
print()
print(f"TDC  =  {round(TDC.X, 2)}  $/day")
print(f" * FCC  =  {round((FCC.X / (alpha * CCF)), 2)}  $/day  ({round(FCC.X, 2)} $ in total)")
print(f" * TCC  =  {round((TCC.X / (alpha * CCF)), 2)}  $/day  ({round(TCC.X, 2)} $ in total)")
print(f" * FOC  =  {round(FOC.X, 2)}  $/day")
print(f" * TOC  =  {round(TOC.X, 2)}  $/day")
print(f" * FSC  =  {round(FSC.X, 2)}  $/day")
print(f"    - FC   =  {round(FC.X, 2)}  $/day")
print(f"    - LC   =  {round(LC.X, 2)}  $/day")
print(f"    - MC   =  {round(MC.X)}  $/day")
print(f"    - GC   =  {round(GC.X, 2)}  $/day")
print()

print("=== CARBON VARIABLES ===")
print(f"TCE      =  {round(TCE.X, 2)} ton CO2/day")
print(f"  * TCEFeed  =  {round(TCEFeed.X, 1)}  ton CO2/day ({round(100*(TCEFeed.X/TCE.X), 1)}%)")
print(f"  * TCEProd  =  {round(TCEProd.X, 2)}  ton CO2/day ({round(100*(TCEProd.X/TCE.X), 1)}%)")
print(f"  * TCETrans =  {round(TCETrans.X, 2)} ton CO2/day ({round(100*(TCETrans.X/TCE.X), 1)}%)")
print()
print("Carbon intensities and carbon emissions per grid: ")
CI_and_Carbon = np.zeros(shape=(len(G), 6))
CI_and_Carbon[:, 0] = CI_ig[0, :]
CI_and_Carbon[:, 1] = CI_ig[1, :]
CI_and_Carbon[:, 2] = CI_g
CI_and_Carbon[:, 3] = Carbon_ig[0, :]
CI_and_Carbon[:, 4] = Carbon_ig[1, :]
CI_and_Carbon[:, 5] = Carbon_g
print(pd.DataFrame(CI_and_Carbon, index=GSet, columns=["CI_[CH2]g", "CI_[LH2]g", "CI_g", "Carbon_[CH2]g", "Carbon_[LH2]g", "Carbon_g"]).to_string())
print()
print("Validity check: ")
print("TCE (right)   =", round(TCE.X, 3))
print("sum(Carbon_g) =", round(sum(Carbon_g), 3))

if printTests:
    print("==================== TEST =============================")

    print("CI_pi: \n", pd.DataFrame(CI_pi, index=PSet, columns=ISet).to_string())
    print()

    print("CIPlants_ig: \n", pd.DataFrame(CIPlants_ig, index=ISet, columns=GSet).to_string())
    print()

    print("CIProd_denominator_ig: \n", pd.DataFrame(CIProd_denominator_ig, index=ISet, columns=GSet).to_string())
    print()

    print("CIProd_ig: \n", pd.DataFrame(CIProd_ig, index=ISet, columns=GSet).to_string())
    print()

    print("CITrans_ig: \n", pd.DataFrame(CITrans_ig, index=ISet, columns=GSet).to_string())
    print()

    print("CI_ig: \n", pd.DataFrame(CI_ig, index=ISet, columns=GSet).to_string())
    print()

    print("C_g: \n", pd.DataFrame(CI_g, index=GSet).T.to_string())
    print()

    print("Carbon_ig: \n", pd.DataFrame(Carbon_ig, index=ISet, columns=GSet).to_string())
    print()

    print("Carbon_g: \n", pd.DataFrame(Carbon_g, index=GSet).T.to_string())
    print()