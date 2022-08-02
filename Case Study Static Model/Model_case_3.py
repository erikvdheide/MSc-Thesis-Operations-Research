""" Model_case_2

NEWLY ADDED: CCS plants.

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
import math
from itertools import product

# Import the data and the parameters
from Data_case_3 import *

############################## Set up model & decision variables ######################################################

""" Initiate model """
Model = gp.Model()
if optimizeCI:
    Model.params.nonConvex = 2
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

if optimizeCI:
    # Intensity of product i produced at grid g
    CIPlants_ig = Model.addMVar(shape=(len(I), len(G)))

    # Intensity of product i satisfied at grid g from production
    CIProd_ig = Model.addMVar(shape=(len(I), len(G)))

    # Intensity of product i satisfied at grid g from transport
    CITrans_ig = Model.addMVar(shape=(len(I), len(G)))

    # Total intensity of product i at grid g (ton CO2/ton H2)
    CI_ig = Model.addMVar(shape=(len(I), len(G)))

    # Total intensity at grid g (ton CO2/ton H2)
    CI_g = Model.addMVar(shape=(len(G)))

############################## Model objective & constraints ##########################################################

""" Set some variables to 0 (ig, il and gg options) """

# No local production if no feasible location (ig)
Model.addConstrs(DL_ig[(i, g)] == 0 for i in I for g in G if A_ig[(i, g)] == 0)

# No production by a plant if no feasible location (ig)
Model.addConstrs(P_pig[(p, i, g)] == 0 for p in P for i in I for g in G if A_ig[(i, g)] == 0)

# No total production in a grid plant ino feasible location (ig)
Model.addConstrs(PT_ig[(i, g)] == 0 for i in I for g in G if A_ig[(i, g)] == 0)

# No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)

# No plants can be build on locations where plants are not feasible (ig)
Model.addConstrs(NP_pig[(p, i, g)] == 0 for p in P for i in I for g in G if A_ig[(i, g)] == 0)

# No transport units if product-mode combination does not exist (il)
Model.addConstrs(NTU_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0)

# No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(X_ilgg[(i, l, g, g_prime)] == 0 for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)

# No export if you cannot produce the product at the location (ig)
Model.addConstrs(Y_ig[(i, g)] == 0 for i in I for l in L for g in G for g_prime in G if A_ig[(i, g)] == 0)

if optimizeCI:
    # Only where you can actually build a plant you can have a carbon intensity of the plants in that grid
    Model.addConstrs(CIPlants_ig[(i, g)] == 0 for i in I for g in G if A_ig[(i, g)] == 0)

""" Operational part """

# (B.1) Not more demand satisfied locally than totally produced in that grid
Model.addConstrs(DL_ig[(i, g)] <= PT_ig[(i, g)] for i in I for g in G if A_ig[(i, g)] == 1)

# (B.2) NEW - local demand satisfied are the "diagonals" on the Q matrix
Model.addConstrs(DL_ig[(i, g)] == gp.quicksum(Q_ilgg[(i, l, g, g)] for l in L if A_il[(i, l)] == 1) for i in I for g in G)

# (B.3) Imported hydrogen (has NO constraint on location)
Model.addConstrs(DI_ig[(i, g)] == gp.quicksum(Q_ilgg[(i, l, g_prime, g)] for l in L for g_prime in G if g_prime != g) for i in I for g in G)

# (B.4) Total demand satisfied
Model.addConstrs(DT_g[g] == gp.quicksum(DL_ig[(i, g)] + DI_ig[(i, g)] for i in I) for g in G)

# (B.5) Balance constraints
Model.addConstrs(DL_ig[(i, g)] == PT_ig[(i, g)] - gp.quicksum(Q_ilgg[(i, l, g, g_prime)] for l in L for g_prime in G if g_prime != g) for i in I for g in G)
# Model.addConstrs(DL_ig[(i, g)] + DI_ig[(i, g)]
#                  == PT_ig[(i, g)] + gp.quicksum(Q_ilgg[(i, l, g_prime, g)] for l in L for g_prime in G if g_prime != g)
#                  - gp.quicksum(Q_ilgg[(i, l, g, g_prime)] for l in L for g_prime in G if g_prime != g) for i in I for g in G)

# (B.6) Total production
Model.addConstrs(PT_ig[(i, g)] == gp.quicksum(P_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)

# (B.7) Bounding production
Model.addConstrs(P_pig[(p, i, g)] >= PCapMin_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G if A_ig[(i, g)] == 1)
Model.addConstrs(P_pig[(p, i, g)] <= PCapMax_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G if A_ig[(i, g)] == 1)

# (B.8) Bounding total production
Model.addConstrs(PT_ig[(i, g)] >= gp.quicksum(PCapMin_pi[(p, i)] * NP_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)
Model.addConstrs(PT_ig[(i, g)] <= gp.quicksum(PCapMax_pi[(p, i)] * NP_pig[(p, i, g)] for p in P) for i in I for g in G if A_ig[(i, g)] == 1)

# (B.9) Bounding transportation
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] >= QMin_il[(i, l)] * X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)
Model.addConstrs(Q_ilgg[(i, l, g, g_prime)] <= QMax_il[(i, l)] * X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)

# (B.10) Import or export, not both
Model.addConstrs(X_ilgg[(i, l, g, g_prime)] + X_ilgg[(i, l, g_prime, g)] <= 1 for i in I for l in L for g in G for g_prime in G if g != g_prime)

# (B.11) Import or export, not both
Model.addConstrs(Y_ig[(i, g)] >= X_ilgg[(i, l, g, g_prime)] for i in I for l in L for g in G for g_prime in G if g != g_prime if A_ig[(i, g)] == 1)

# (B.12) Import or export, not both
Model.addConstrs(Z_ig[(i, g)] >= X_ilgg[(i, l, g_prime, g)] for i in I for l in L for g_prime in G for g in G if g != g_prime)

# (B.13) Import or export, not both
Model.addConstrs(Y_ig[(i, g)] + Z_ig[(i, g)] <= 1 for i in I for g in G)

# (B.14) NEW - cannot import both products (does not entirely fix issue!)
if onlyOneImportProduct:
    Model.addConstrs(gp.quicksum(Z_ig[(i, g)] for i in I) <= 1 for g in G)

""" Cost part """

# (B.16) Facility Capital Cost (FCC)
Model.addConstr(FCC == gp.quicksum(PCC_pi[(p, i)] * NP_pig[(p, i, g)] for p in P for i in I for g in G))

# (B.17) Number Transport Units between grids (NTU_ilgg')
Model.addConstrs(NTU_ilgg[(i, l, g, g_prime)] >= Q_ilgg[(i, l, g, g_prime)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l])
                 for g in G for g_prime in G for i in I for l in L if g_prime != g if A_il[(i, l)] == 1)

# (B.18) Number Transport Units within grids (NTU_ilgg)
Model.addConstrs(NTU_ilgg[(i, l, g, g)] >= Q_ilgg[(i, l, g, g)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l])
                 for i in I for l in L for g in G if A_il[(i, l)] == 1)

# (B.19) Number Transport Units per product and mode (NTU_il)
Model.addConstrs(NTU_il[(i, l)] == gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] for g in G for g_prime in G) for i in I for l in L if A_il[(i, l)] == 1)

# (B.20) Transportation Capital Cost (TCC)
Model.addConstr(TCC == gp.quicksum(TMC_il[(i, l)] * NTU_il[(i, l)] for i in I for l in L))

# (B.21) Facility Operating Cost (FOC)
if includeCCS:
    Model.addConstr(FOC == gp.quicksum((UPC_pi[(p, i)] + A_CCS_p[p] * CCSCost * CO2Prod_pi[(p, i)]) * P_pig[(p, i, g)] for p in P for i in I for g in G))
else:
    Model.addConstr(FOC == gp.quicksum(UPC_pi[(p, i)] * P_pig[(p, i, g)] for p in P for i in I for g in G))

# (B.22) Fuel Cost (FC)
Model.addConstr(FC == gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilgg[(i, l, g, g_prime)] * (1 / FEbetween_l[l]) * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                + gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g)] * Q_ilgg[(i, l, g, g)] * (1 / FEwithin_l[l]) * (1 / TCap_il[(i, l)]))
                              for i in I for l in L for g in G if A_il[(i, l)] == 1))

# (B.23) Labour Cost (LC)
Model.addConstr(LC == gp.quicksum(DW_l[l] * (Q_ilgg[(i, l, g, g_prime)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                + gp.quicksum(DW_l[l] * (Q_ilgg[(i, l, g, g)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l]))
                              for i in I for l in L for g in G if A_il[(i, l)] == 1))

# (B.24) Maintenance Cost (MC)
Model.addConstr(MC == gp.quicksum(ME_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilgg[(i, l, g, g_prime)] * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1))

# (B.25) General Cost (GC)
Model.addConstr(GC == gp.quicksum(GE_l[l] * NTU_il[(i, l)] for i in I for l in L))

# (B.26) - NEW - FeedStock Cost (FSC)
Model.addConstr(FSC == gp.quicksum(FSP_p[p] * PCR_p[p] * gp.quicksum(P_pig[(p, i, g)] for i in I for g in G) for p in P))

# (B.27) Transportation Operating Cost (TOC)
Model.addConstr(TOC == FC + LC + MC + GC)

# (B.28) Total Daily Cost (TDC)
Model.addConstr(TDC == (1 / (alpha * CCF)) * (FCC + TCC) + FOC + TOC + FSC)

# (B.29) Minimize TDC or TCE, or both
objective = costImportance * TDC + CO2Importance * TCE

""" CO2 part """

# (B.30) Total feed emissions
Model.addConstr(TCEFeed == gp.quicksum(CO2Feed_p[p] * (gp.quicksum(P_pig[(p, i, g)] for g in G for i in I)) for p in P))

# (B.31) Total production emissions
if includeCCS:
    Model.addConstr(TCEProd == gp.quicksum(((1 - A_CCS_p[p]) * CO2Prod_pi[(p, i)] + A_CCS_p[p] * (1 - CCSEF) * CO2Prod_pi[(p, i)]) * P_pig[(p, i, g)] for p in P for i in I for g in G))
else:
    Model.addConstr(TCEProd == gp.quicksum(CO2Prod_pi[(p, i)] * (gp.quicksum(P_pig[(p, i, g)] for g in G)) for p in P for i in I))

# (B.32) Total transportation emissions
Model.addConstr(TCETrans == CO2Trans * gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] * 2 * L_lgg[(l, g, g_prime)] for i in I for l in L for g in G for g_prime in G))
# Model.addConstr(TCETrans == gp.quicksum(CO2Trans * (gp.quicksum(NTU_ilgg[(i, l, g, g_prime)] * 2 * L_lgg[(l, g, g_prime)] for i in I for l in L)) for g in G for g_prime in G))

# (B.33) Total chain emissions
Model.addConstr(TCE == TCEFeed + TCEProd + TCETrans)

# (B.34) Linear constraint on total chain emissions
if TotalCO2Max != -1:
    Model.addConstr(TCE <= TotalCO2Max)

##################### Restrict Carbon intensities DURING-optimization ###########################################

# (B.46)-(B.53) Include CI during-optimization
if optimizeCI:

    # (B.43) Calculate the starting carbon intensity as parameter per plant p per product i
    pass  # already in data

    # (B.44) CI^{Plants}_{ig}
    Model.addConstrs(CIPlants_ig[(i, g)] * PT_ig[(i, g)] == gp.quicksum(P_pig[(p, i, g)] * CIStart_pi[(p, i)] for p in P) for i in I for g in G)

    # (B.45) CI^{Prod}_{ig}
    Model.addConstrs(CIProd_ig[(i, g)] * gp.quicksum(Q_ilgg[(i, l_star, g_prime_star, g)] for l_star in L for g_prime_star in G)
                     == gp.quicksum(Q_ilgg[(i, l, g_prime, g)] * CIPlants_ig[(i, g_prime)] for l in L for g_prime in G) for i in I for g in G)

    # (B.46) CI^{Trans}_{ig}
    Model.addConstrs(CITrans_ig[(i, g)] * (DL_ig[(i, g)] + DI_ig[(i, g)])
                     == CO2Trans * gp.quicksum(NTU_ilgg[(i, l, g_prime, g)] * 2 * L_lgg[(l, g_prime, g)] for l in L for g_prime in G) for i in I for g in G)

    # (B.47) CI_{ig}
    Model.addConstrs(CI_ig[(i, g)] == CIProd_ig[(i, g)] + CITrans_ig[(i, g)] for i in I for g in G)

    # (B.48) CI_{g}
    Model.addConstrs(CI_g[g] == gp.quicksum(CI_ig[(i, g)] * (DL_ig[(i, g)] + DI_ig[(i, g)]) * (1 / DT_g[g]) for i in I) for g in G)

    # (B.49) CI_{ig} restrictions
    for i in I:
        for g in G:
            if maxCI[(i, g)] != -1:
                Model.addConstr(CI_ig[(i, g)] <= maxCI[(i, g)])

    # (B.50) CI_{g} restrictions
    for g in G:
        if maxCI[(len(I), g)] != -1:
            Model.addConstr(CI_g[g] <= maxCI[(len(I), g)])

############################## Execute ########################################################################

""" Execute the model as specified in Params. """
start_time = time.time()
Model.setObjective(objective, GRB.MINIMIZE)
Model.optimize()
end_time = time.time()
total_time = end_time - start_time

####################### Calculate Carbon intensities POST-optimization ##########################################

""" Options: """
# Carbon_ig_post : total daily CO2 emissions for product i in grid g (per day)
# Carbon_g_post  : total daily CO2 emissions for grid g (per day)
# CI_ig_post     : carbon intensity of product i in grid g (at a day)
# CI_g_post      : avg. carbon intensity of utilized products in grid g

# (B.35) Calculate the starting carbon intensity as parameter per plant p per product i
pass  # already in data

# (B.36) Calculate the carbon intensities of the all the plants producing product i in plant g
# Note that the simplifying assumption here is that what is produced totally in a grid is distributed equally over the other grids.
# ..in other words, if you have one CG plant in a grid and one SMR plant in the same grid, you average the emissions of those grids for the
# ..outgoing products, while in fact you could connect some cities with the CG plant and some cities with the SMR plant
# ..BUT this only holds for the same product i! Depends whether you take g, (i, g) or (p, i, g)
CIPlants_ig_post = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        if PT_ig[(i, g)].X > 0.0000001:
            CIPlants_ig_post[(i, g)] = gp.quicksum(P_pig[(p, i, g)] / PT_ig[(i, g)].X * CIStart_pi[(p, i)] for p in P).getValue()

# (B.37) Carbon intensities from production units for product i at customer location g
CIProd_ig_post = np.zeros(shape=(len(I), len(G)))
CIProd_denominator_ig = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        CIProd_denominator_ig[(i, g)] = gp.quicksum(Q_ilgg[(i, l_star, g_prime_star, g)] for l_star in L for g_prime_star in G).getValue()
for i in I:
    for g in G:
        if CIProd_denominator_ig[(i, g)] > 0.0000001:
            CIProd_ig_post[(i, g)] = gp.quicksum(Q_ilgg[(i, l, g_prime, g)] / CIProd_denominator_ig[(i, g)] * CIPlants_ig_post[(i, g_prime)] for l in L for g_prime in G).getValue()

# (B.38) Carbon intensity from transportation
CITrans_ig_post = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        if (DL_ig[(i, g)].X + DI_ig[(i, g)].X) > 0.0000001:
            CITrans_ig_post[(i, g)] = ((CO2Trans * gp.quicksum(NTU_ilgg[(i, l, g_prime, g)] * 2 * L_lgg[(l, g_prime, g)] for l in L for g_prime in G)) / (DL_ig[(i, g)].X + DI_ig[(i, g)].X)).getValue()

# (B.39) Total carbon intensity is part from production and part from transportation
CI_ig_post = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        CI_ig_post[(i, g)] = CIProd_ig_post[(i, g)] + CITrans_ig_post[(i, g)]

# (B.40) Carbon intensity at a location is the resource-weighted average of both products
CI_g_post = np.zeros(shape=len(G))
for g in G:
    CI_g_post[g] = gp.quicksum(CI_ig_post[(i, g)] * (DL_ig[(i, g)] + DI_ig[(i, g)]) / DT_g[g] for i in I).getValue()

# (B.41) Total carbon emissions for which grid g is responsible, per product
Carbon_ig_post = np.zeros(shape=(len(I), len(G)))
for i in I:
    for g in G:
        Carbon_ig_post[(i, g)] = CI_ig_post[(i, g)] * (DL_ig[(i, g)].X + DI_ig[(i, g)].X)

# (B.42) Total carbon emissions for which grid g is responsible, all products combined
Carbon_g_post = np.zeros(shape=(len(G)))
for g in G:
    Carbon_g_post[g] = CI_g_post[g] * DT_g[g]


############################## Print the results #########################################################

def deleteZeroMatrix(mat):
    # Checks if a data matrix only has 0 elements
    for elem in np.nditer(mat):
        if elem > 0.0001:
            return False
    return True


def deleteZeroElements(mat):
    mat2 = np.copy(mat).astype('object')
    # Remove the elements which are 0.0
    if len(mat.shape) == 1:
        pass
        numElem = mat2.shape[0]
        for i in range(numElem):
            if -0.0001 < mat2[i] < 0.0001:
                mat2[i] = '-'
    else:
        numRows = mat2.shape[0]
        numCols = mat2.shape[1]
        for i in range(numRows):
            for j in range(numCols):
                if -0.0001 < mat2[(i, j)] < 0.0001:
                    mat2[(i, j)] = '-'
    return mat2


if printVariables:
    print()
    print("===== VARIABLES =====")
    print("Num. variables   : ", Model.numVars)  # (including those that are by definition 0)
    print("Num. constraints : ", Model.numConstrs)  # (including those that make variables 0 and redundant constraints for 0 variables)
    print()

    print("=== OPERATIONAL VARIABLES ===")
    print()

    # Continuous variables
    matrix = DL_ig.X.round(2)
    if not deleteZeroMatrix(matrix):  # if it is a zero matrix
        test = pd.DataFrame(deleteZeroElements(DL_ig.X.round(2)), index=ISet, columns=GSet).to_string()
        print("D^L_ig (local demand satisfied): \n", test)
        print()

    matrix = DI_ig.X.round(2)
    if not deleteZeroMatrix(matrix):
        print("D^I_ig (imported demand satisfied): \n", pd.DataFrame(deleteZeroElements(DI_ig.X.round(2)), index=ISet, columns=GSet).to_string())
        print()

    for p in P:
        matrix = P_pig[p, :, :].X.round(2)
        if not deleteZeroMatrix(matrix):
            print(f"P_[{PSet[p]}]ig (production from plant p={PSet[p]}): \n {pd.DataFrame(deleteZeroElements(P_pig[p, :, :].X.round(2)), index=ISet, columns=GSet).to_string()}")
            print()

    matrix = PT_ig.X.round(2)
    if not deleteZeroMatrix(matrix):
        print("P^T_ig (total production): \n", pd.DataFrame(deleteZeroElements(PT_ig.X.round(2)), index=ISet, columns=GSet).to_string())
        print()

    for i in I:
        for l in L:
            if A_il[(i, l)] == 1:
                matrix = Q_ilgg[i, l, :, :].X.round(2)
                if not deleteZeroMatrix(matrix):
                    print(f"Q_[{ISet[i]}][{LSet[l]}]gg' (transport product {ISet[i]} using mode {LSet[l]}): "
                          f"\n {pd.DataFrame(deleteZeroElements(Q_ilgg[i, l, :, :].X.round(2)), index=GSet, columns=GSet).to_string()}")
                    print()

    # Total feedstock needed (ambiguous):
    FeedNeeded_p_large = np.zeros(shape=len(P))
    for p in P:
        FeedNeeded_p_large[p] = (PCR_p[p] * gp.quicksum(P_pig[(p, i, g)] for i in I for g in G)).getValue()
    print("FeedNeeded_p \n", pd.DataFrame(deleteZeroElements(FeedNeeded_p_large), index=PSet, columns=["Amount of feed"]).to_string())
    print()

    # Integer variables
    for p in P:
        matrix = NP_pig[p, :, :].X.round(2)
        if not deleteZeroMatrix(matrix):
            print(f"NP_[{PSet[p]}]ig (number of plant {PSet[p]}): \n {pd.DataFrame(deleteZeroElements(NP_pig[p, :, :].X.round(2)), index=ISet, columns=GSet).to_string()}")
            print()

    for i in I:
        for l in L:
            if A_il[(i, l)] == 1:
                matrix = NTU_ilgg[i, l, :, :].X.round(2)
                if not deleteZeroMatrix(matrix):
                    print(f"NTU_[{ISet[i]}][{LSet[l]}]gg (#units transport product i={ISet[i]} using mode l={LSet[l]}): \n"
                          f" {pd.DataFrame(deleteZeroElements(NTU_ilgg[i, l, :, :].X.round(2)), index=GSet, columns=GSet).to_string()}")
                    print()

    print("TOTAL NTU_il (number of units): \n", pd.DataFrame(deleteZeroElements(NTU_il.X.round(2)), index=ISet, columns=LSet).to_string())
    print()

    # Binary variables
    if printBinaryVariables:
        for i in I:
            for l in L:
                if A_il[(i, l)] == 1:
                    matrix = X_ilgg[i, l, :, :].X
                    if not deleteZeroMatrix(matrix):
                        print(f"X_[{ISet[i]}][{LSet[l]}]gg' (binary for transport product {ISet[i]} using mode {LSet[l]}): "
                              f"\n {pd.DataFrame(deleteZeroElements(X_ilgg[i, l, :, :].X), index=GSet, columns=GSet).to_string()}")
                        print()

        print("Y_ig (binary export): \n", pd.DataFrame(deleteZeroElements(Y_ig.X), index=ISet, columns=GSet).to_string())
        print()

        print("Z_ig (binary import): \n", pd.DataFrame(deleteZeroElements(Z_ig.X), index=ISet, columns=GSet).to_string())
        print()

print()
print("===== RESULT =====")
print(f"Objective  : {round(objective.getValue(), 3)} (costImportance={costImportance}, CO2Importance={CO2Importance})")
print(f"Total time : {round(total_time, 3)} sec")
print()

print("=== COST VARIABLES ===")
print(f"TDC (Total Daily Cost)            =  {round(TDC.X, 2)}  $/day")
print(f" * FCC (Facility Capital Cost)    =  {round((FCC.X / (alpha * CCF)), 2)}  $/day  ({round(FCC.X, 2)} $ in total)")
print(f" * TCC (Total Capital Cost)       =  {round((TCC.X / (alpha * CCF)), 2)}  $/day  ({round(TCC.X, 2)} $ in total)")
print(f" * FOC (Facility Operating Cost)  =  {round(FOC.X, 2)}  $/day")
print(f" * FSC (FeedStock Cost)           =  {round(FSC.X, 2)}  $/day")
print(f" * TOC (Transport Operating Cost) =  {round(TOC.X, 2)}  $/day")
print(f"    - FC (Fueling Cost)           =  {round(FC.X, 2)}  $/day")
print(f"    - LC (Labour Cost)            =  {round(LC.X, 2)}  $/day")
print(f"    - MC (Maintenance Cost)       =  {round(MC.X)}  $/day")
print(f"    - GC (General Cost)           =  {round(GC.X, 2)}  $/day")
print()

print("=== CARBON VARIABLES ===")
print(f"TCE (Total Carbon Emission)  =  {round(TCE.X, 2)} ton CO2/day")
print(f"  * TCEFeed (Feedstock Em.)  =  {round(TCEFeed.X, 1)}  ton CO2/day ({round(100 * (TCEFeed.X / TCE.X), 1)}%)")
print(f"  * TCEProd (Production Em.) =  {round(TCEProd.X, 2)}  ton CO2/day ({round(100 * (TCEProd.X / TCE.X), 1)}%)")
print(f"  * TCETrans (Transport Em.) =  {round(TCETrans.X, 2)} ton CO2/day ({round(100 * (TCETrans.X / TCE.X), 1)}%)")
print()
print("POST-optimization Carbon intensities and total carbon emissions per grid: ")
CI_and_Carbon = np.zeros(shape=(len(G), 6))
CI_and_Carbon[:, 0] = CI_ig_post[0, :]
CI_and_Carbon[:, 1] = CI_ig_post[1, :]
CI_and_Carbon[:, 2] = CI_g_post
CI_and_Carbon[:, 3] = Carbon_ig_post[0, :]
CI_and_Carbon[:, 4] = Carbon_ig_post[1, :]
CI_and_Carbon[:, 5] = Carbon_g_post
print(pd.DataFrame(CI_and_Carbon, index=GSet, columns=["CI_[CH2]g", "CI_[LH2]g", "CI_g", "Carbon_[CH2]g", "Carbon_[LH2]g", "Carbon_g"]).to_string())
print()

if optimizeCI:

    def correctCI(CI_ig_input, DI_ig_input, DL_ig_input):
        """ Function to set CI values to 0 is there is no satisfied demand
        :param CI_ig_input: Carbon intensity of product i and grid g
        :param DI_ig_input: Satisfied demand imported of product i in grid g
        :param DL_ig_input: Satisfied demand locally satisfied of product in in grid g."""
        for i in I:
            for g in G:
                if (DI_ig_input[(i, g)] + DL_ig_input[(i, g)] < 0.0001) and (CI_ig_input[(i, g)] > 0):
                    CI_ig_input[(i, g)] = 0
        return CI_ig_input


    CI_ig_corr = correctCI(CI_ig.X, DI_ig.X, DL_ig.X)

    print("DURING-optimization Carbon intensities and total carbon emissions per grid: ")
    CI_and_Carbon = np.zeros(shape=(len(G), 6))
    CI_and_Carbon[:, 0] = CI_ig_corr[0, :]
    CI_and_Carbon[:, 1] = CI_ig_corr[1, :]
    CI_and_Carbon[:, 2] = CI_g.X
    CI_and_Carbon[:, 3] = CI_ig.X[0, :] * (DL_ig.X[0, :] + DI_ig.X[0, :])
    CI_and_Carbon[:, 4] = CI_ig.X[1, :] * (DL_ig.X[1, :] + DI_ig.X[1, :])
    CI_and_Carbon[:, 5] = CI_g.X * DT_g
    print(pd.DataFrame(CI_and_Carbon, index=GSet, columns=["CI_[CH2]g", "CI_[LH2]g", "CI_g", "Carbon_[CH2]g", "Carbon_[LH2]g", "Carbon_g"]).to_string())
    print()

print("Validity check: ")
print("TCE                        =", round(TCE.X, 3))
print("POST:   sum(Carbon_g_post) =", round(sum(Carbon_g_post), 3))
if optimizeCI:
    print("DURING: sum(Carbon_g)      =", round(sum(CI_g.X * DT_g), 3))
print()

if printTests:
    print("==================== TEST =============================")

    print("CIStart_pi: \n", pd.DataFrame(CIStart_pi, index=PSet, columns=ISet).to_string())
    print()

    print("CIPlants_ig_post: \n", pd.DataFrame(CIPlants_ig_post, index=ISet, columns=GSet).to_string())
    print()

    print("CIProd_denominator_ig: \n", pd.DataFrame(CIProd_denominator_ig, index=ISet, columns=GSet).to_string())
    print()

    print("CIProd_ig_post: \n", pd.DataFrame(CIProd_ig_post, index=ISet, columns=GSet).to_string())
    print()

    print("CITrans_ig_post: \n", pd.DataFrame(CITrans_ig_post, index=ISet, columns=GSet).to_string())
    print()

    print("CI_ig_post: \n", pd.DataFrame(CI_ig_post, index=ISet, columns=GSet).to_string())
    print()

    print("C_g: \n", pd.DataFrame(CI_g_post, index=GSet).T.to_string())
    print()

    print("Carbon_ig_post: \n", pd.DataFrame(Carbon_ig_post, index=ISet, columns=GSet).to_string())
    print()

    print("Carbon_g_post: \n", pd.DataFrame(Carbon_g_post, index=GSet).T.to_string())
    print()
