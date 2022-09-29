""" Model_case_dynamic

This is the dynamic model for the HSCN case study.

This file contains the following components:
- Variable declaration + variables to 0;
- Exact model - operational part;
- Exact model - cost part;
- Exact model - linear CO2 part;
- Exact model - nonlinear CI part;
- Calculate carbon intensity post-optimization.

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
from Data_HSCN import *

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
DL_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=0.0)

# Demand imported
DI_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=0.0)

# Amount produced
P_pigt = Model.addMVar(shape=(len(P), len(I), len(G), len(T)), lb=0.0)

# Total production
PT_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=0.0)

# Flows between grids
Q_ilggt = Model.addMVar(shape=(len(I), len(L), len(G), len(G), len(T)), lb=0.0)

# Cost variables
FC_t = Model.addMVar(shape=(len(T)), lb=0.0)
FCC_t = Model.addMVar(shape=(len(T)), lb=0.0)
FOC_t = Model.addMVar(shape=(len(T)), lb=0.0)
FSC_t = Model.addMVar(shape=(len(T)), lb=0.0)
GC_t = Model.addMVar(shape=(len(T)), lb=0.0)
LC_t = Model.addMVar(shape=(len(T)),lb=0.0)
MC_t = Model.addMVar(shape=(len(T)), lb=0.0)
TCC_t = Model.addMVar(shape=(len(T)), lb=0.0)
TCE_t = Model.addMVar(shape=(len(T)), lb=-10e6)
TCEFeed_t = Model.addMVar(shape=(len(T)), lb=0.0)
TCEProd_t = Model.addMVar(shape=(len(T)), lb=-10e6)
TCETrans_t = Model.addMVar(shape=(len(T)), lb=0.0)
TDC_t = Model.addMVar(shape=(len(T)), lb=0.0)
TOC_t = Model.addMVar(shape=(len(T)), lb=0.0)

TDC = Model.addVar(lb=0.0)
TCE = Model.addVar(lb=-10e6)

# Number of trips daily (not rounded)
NOT_ilggt = Model.addMVar(shape=(len(I), len(L), len(G), len(G), len(T)))

# Number of plants
NP_pigt = Model.addMVar(shape=(len(P), len(I), len(G), len(T)), vtype=GRB.INTEGER)

# Investments of plants in a time period
IP_pigt = Model.addMVar(shape=(len(P), len(I), len(G), len(T)), vtype=GRB.INTEGER)

# Upgrades of CCS
UP_pigt = Model.addMVar(shape=(len(P_prime), len(I), len(G), len(T)), vtype=GRB.INTEGER)

if overEstimateNTU:
    # Number of transport units from grid to grid
    NTU_iltggt = Model.addMVar(shape=(len(I), len(L), len(G), len(G), len(T)), vtype=GRB.INTEGER)

# Number of transport units total
NTU_ilt = Model.addMVar(shape=(len(I), len(L), len(T)), vtype=GRB.INTEGER)

# Investments of transport units in a time period
ITU_ilt = Model.addMVar(shape=(len(I), len(L), len(T)), vtype=GRB.INTEGER)

# Binary for some transport
X_ilggt = Model.addMVar(shape=(len(I), len(L), len(G), len(G), len(T)), vtype=GRB.BINARY)

# Binary for export
Y_igt = Model.addMVar(shape=(len(I), len(G), len(T)), vtype=GRB.BINARY)

# Binary for import
Z_igt = Model.addMVar(shape=(len(I), len(G), len(T)), vtype=GRB.BINARY)

if optimizeCI:
    # Intensity of product i produced at grid g
    CIPlants_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=-10e6)

    # Intensity of product i satisfied at grid g from production
    CIProd_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=-10e6)

    # Intensity of product i satisfied at grid g from transport
    CITrans_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=-10e6)

    # Total intensity of product i at grid g (ton CO2/ton H2)
    CI_igt = Model.addMVar(shape=(len(I), len(G), len(T)), lb=-10e6)

    # Total intensity at grid g (ton CO2/ton H2)
    CI_gt = Model.addMVar(shape=(len(G), len(T)), lb=-10e6)

############################## Model objective & constraints ##########################################################

""" Set some variables to 0 (ig, il and gg options) """

# (B.69)-1 No local production if no feasible location (ig)
Model.addConstrs(DL_igt[(i, g, t)] == 0 for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-2 No production by a plant if no feasible location (ig)
Model.addConstrs(P_pigt[(p, i, g, t)] == 0 for p in P for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-3 No total production in a grid plant ino feasible location (ig)
Model.addConstrs(PT_igt[(i, g, t)] == 0 for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-4 No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(Q_ilggt[(i, l, g, g_prime, t)] == 0 for i in I for l in L for g in G for g_prime in G for t in T if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)

# (B.69)-5 No plants can be build on locations where plants are not feasible (ig)
Model.addConstrs(NP_pigt[(p, i, g, t)] == 0 for p in P for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-6 No plants can be build on locations where plants are not feasible (ig)
Model.addConstrs(IP_pigt[(p, i, g, t)] == 0 for p in P for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-6 No plants can be build on locations where plants are not feasible (ig)
Model.addConstrs(UP_pigt[(p_prime, i, g, t)] == 0 for p_prime in P_prime for i in I for g in G for t in T if A_ig[(i, g)] == 0)

# (B.69)-7 No transport if: product cannot go with the mode of transport (il), product cannot be produced in the grid (ig)
Model.addConstrs(X_ilggt[(i, l, g, g_prime, t)] == 0 for i in I for l in L for g in G for g_prime in G for t in T if A_il[(i, l)] == 0 or A_ig[(i, g)] == 0)

# (B.69)-8 No export if you cannot produce the product at the location (ig)
Model.addConstrs(Y_igt[(i, g, t)] == 0 for i in I for l in L for g in G for g_prime in G for t in T if A_ig[(i, g)] == 0)

# (il)-1 No transport units if product-mode combination does not exist (il)
Model.addConstrs(NTU_ilt[(i, l, t)] == 0 for i in I for l in L for t in T if A_il[(i, l)] == 0)

# (il)-2 No transport units if product-mode combination does not exist (il)
Model.addConstrs(ITU_ilt[(i, l, t)] == 0 for i in I for l in L for t in T if A_il[(i, l)] == 0)

# (il)-3 No trips if production-mode combination does not exist (il)
Model.addConstrs(NOT_ilggt[(i, l, g, g_prime, t)] == 0 for i in I for l in L for g in G for g_prime in G for t in T if A_il[(i, l)] == 0)

if optimizeCI:
    # (B.112) Only where you can actually build a plant you can have a carbon intensity of the plants in that grid
    Model.addConstrs(CIPlants_igt[(i, g, t)] == 0 for i in I for g in G for t in T if A_ig[(i, g)] == 0)

    # (B.113) CI is 0 if there is also no demand
    for g in G:
        for t in T:
            if DT_g[g] == 0:
                Model.addConstr(CI_gt[(g, t)] == 0)
                for i in I:
                    Model.addConstr(CI_igt[(i, g, t)] == 0)


""" Operational part """

# (B.52) Not more demand satisfied locally than totally produced in that grid
Model.addConstrs(DL_igt[(i, g, t)] <= PT_igt[(i, g, t)] for i in I for g in G for t in T if A_ig[(i, g)] == 1)

# (B.53) local demand satisfied are the "diagonals" on the Q matrix (NEW)
Model.addConstrs(DL_igt[(i, g, t)] == gp.quicksum(Q_ilggt[(i, l, g, g, t)] for l in L if A_il[(i, l)] == 1) for i in I for g in G for t in T)

# (B.54) Imported hydrogen (has NO constraint on location)
Model.addConstrs(DI_igt[(i, g, t)] == gp.quicksum(Q_ilggt[(i, l, g_prime, g, t)] for l in L for g_prime in G if g_prime != g) for i in I for g in G for t in T)

# (B.55) Total demand satisfied
Model.addConstrs(gp.quicksum(DL_igt[(i, g, t)] + DI_igt[(i, g, t)] for i in I) == DT_gt[(g, t)] for g in G for t in T)

# (B.56) Balance constraints
Model.addConstrs(DL_igt[(i, g, t)] == PT_igt[(i, g, t)] - gp.quicksum(Q_ilggt[(i, l, g, g_prime, t)] for l in L for g_prime in G if g_prime != g) for i in I for g in G for t in T)

# (B.57) Total production
Model.addConstrs(PT_igt[(i, g, t)] == gp.quicksum(P_pigt[(p, i, g, t)] for p in P) for i in I for g in G for t in T if A_ig[(i, g)] == 1)

# (B.58) Bounding production
Model.addConstrs(P_pigt[(p, i, g, t)] >= PCapMin_pi[(p, i)] * NP_pigt[(p, i, g, t)] for p in P for i in I for g in G for t in T if A_ig[(i, g)] == 1)
Model.addConstrs(P_pigt[(p, i, g, t)] <= PCapMax_pi[(p, i)] * NP_pigt[(p, i, g, t)] for p in P for i in I for g in G for t in T if A_ig[(i, g)] == 1)

# (B.59) Bounding total production
Model.addConstrs(PT_igt[(i, g, t)] >= gp.quicksum(PCapMin_pi[(p, i)] * NP_pigt[(p, i, g, t)] for p in P) for i in I for g in G for t in T if A_ig[(i, g)] == 1)
Model.addConstrs(PT_igt[(i, g, t)] <= gp.quicksum(PCapMax_pi[(p, i)] * NP_pigt[(p, i, g, t)] for p in P) for i in I for g in G for t in T if A_ig[(i, g)] == 1)

# (B.60) Bounding transportation
Model.addConstrs(Q_ilggt[(i, l, g, g_prime, t)] >= QMin_il[(i, l)] * X_ilggt[(i, l, g, g_prime, t)] for i in I for l in L for g in G for g_prime in G for t in T
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)
Model.addConstrs(Q_ilggt[(i, l, g, g_prime, t)] <= QMax_il[(i, l)] * X_ilggt[(i, l, g, g_prime, t)] for i in I for l in L for g in G for g_prime in G for t in T
                 if A_il[(i, l)] == 1 and A_ig[(i, g)] == 1)

# (B.61) Import or export, not both
Model.addConstrs(X_ilggt[(i, l, g, g_prime, t)] + X_ilggt[(i, l, g_prime, g, t)] <= 1 for i in I for l in L for g in G for g_prime in G for t in T if g != g_prime)

# (B.62) Import or export, not both
Model.addConstrs(Y_igt[(i, g, t)] >= X_ilggt[(i, l, g, g_prime, t)] for i in I for l in L for g in G for g_prime in G for t in T if g != g_prime if A_ig[(i, g)] == 1)

# (B.63) Import or export, not both
Model.addConstrs(Z_igt[(i, g, t)] >= X_ilggt[(i, l, g_prime, g, t)] for i in I for l in L for g_prime in G for g in G for t in T if g != g_prime)

# (B.64) Import or export, not both
Model.addConstrs(Y_igt[(i, g, t)] + Z_igt[(i, g, t)] <= 1 for i in I for g in G for t in T)

# (B.65) cannot import both products (though production grids can still have both) (NEW)
if onlyOneImportProduct:
    Model.addConstrs(gp.quicksum(Z_igt[(i, g, t)] for i in I) <= 1 for g in G for t in T)

# (B.66)-NEW: satisfy the number of plants for the first period
Model.addConstrs(NP_pigt[(p, i, g, 0)] == NP0_pig[(p, i, g)] + IP_pigt[(p, i, g, 0)] for p in P for i in I for g in G)

# (B.67)-NEW: satisfy the number of plants in the other periods
# Note that as the not-CCS plants have the same indices as the p_prime set, we can just insert p here
Model.addConstrs(NP_pigt[(p, i, g, t)] == NP_pigt[(p, i, g, t-1)] + IP_pigt[(p, i, g, t)] - UP_pigt[(p, i, g, t)]
                 for p in P for i in I for g in G for t in T if t != 0 if A_CCS_p[p] == 0)

# (B.68)-NEW: satisfy the number of plants in the other periods
# Note that now we have to fill in for p_prime p - length of P_prime
Model.addConstrs(NP_pigt[(p, i, g, t)] == NP_pigt[(p, i, g, t-1)] + IP_pigt[(p, i, g, t)] + UP_pigt[(p-len(P_prime), i, g, t)]
                 for p in P for i in I for g in G for t in T if t != 0 if A_CCS_p[p] == 1)

""" Cost part """

# (B.70) Facility Capital Cost (FCC_t1) in the first period
Model.addConstr(FCC_t[0] == gp.quicksum(PCC_pi[(p, i)] * NP0_pig[(p, i, g)] for p in P for i in I for g in G)
                + gp.quicksum(PCC_pi[(p, i)] * IP_pigt[(p, i, g, 0)] for p in P for i in I for g in G))

# (B.71) Facility Capital Cost (FCC_t) in the first period
Model.addConstrs(FCC_t[t] == gp.quicksum(PCC_pi[(p, i)] * IP_pigt[(p, i, g, t)] for p in P for i in I for g in G) for t in T if t != 0)

# (B.72) Facility Operating Cost (FOC_t)
if includeCCS:
    Model.addConstrs(FOC_t[t] == gp.quicksum((UPC_pi[(p, i)] + A_CCS_p[p] * CCSCost * CO2Prod_pi[(p, i)]) * P_pigt[(p, i, g, t)] for p in P for i in I for g in G) for t in T)
else:
    Model.addConstrs(FOC_t[t] == gp.quicksum(UPC_pi[(p, i)] * P_pigt[(p, i, g, t)] for p in P for i in I for g in G) for t in T)

# (B.73) FeedStock Cost (FSC_t)
Model.addConstrs(FSC_t[t] == gp.quicksum(FSP_p[p] * PCR_p[p] * gp.quicksum(P_pigt[(p, i, g, t)] for i in I for g in G) for p in P) for t in T)

# (B.74) Define the number of trips (NOT_ilggt)
Model.addConstrs(NOT_ilggt[(i, l, g, g_prime, t)] == Q_ilggt[(i, l, g, g_prime, t)] * (1 / TCap_il[(i, l)]) for i in I for l in L for g in G for g_prime in G for t in T if A_il[(i, l)] == 1)

# (B.75) Number Transport Units (NTU_ilt1) in the first period
Model.addConstrs(NTU_ilt[(i, l, 0)] == NTU0_il[(i, l)] + ITU_ilt[(i, l, 0)] for i in I for l in L)

# (B.76) Number Transport Units (NTU_ilt) in the other periods
Model.addConstrs(NTU_ilt[(i, l, t)] == NTU_ilt[(i, l, t-1)] + ITU_ilt[(i, l, t)] for i in I for l in L for t in T if t != 0)

# (B.77) Invested Transport Units (ITU_ilt1) for first time period
Model.addConstrs(ITU_ilt[(i, l, 0)] >= gp.quicksum(Q_ilggt[(i, l, g, g_prime, 0)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l])
                                       for g in G for g_prime in G if g_prime != g)
                                     + gp.quicksum(Q_ilggt[(i, l, g, g, 0)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l])
                                       for g in G)
                                     - NTU0_il[(i, l)] for i in I for l in L if A_il[(i, l)] == 1)

# (B.78) Invested Transport Units (ITU_ilt) for other time periods
Model.addConstrs(ITU_ilt[(i, l, t)] >= gp.quicksum(Q_ilggt[(i, l, g, g_prime, t)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l])
                                       for g in G for g_prime in G if g_prime != g)
                                     + gp.quicksum(Q_ilggt[(i, l, g, g, t)] * (1 / TMA_l[l]) * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l])
                                       for g in G)
                                     - NTU_ilt[(i, l, t-1)] for i in I for l in L for t in T if A_il[(i, l)] == 1 if t != 0)

# (B.79) Transportation Capital Cost (TCC_t) in first period
Model.addConstr(TCC_t[0] == gp.quicksum(TMC_il[(i, l)] * (NTU0_il[(i, l)] + ITU_ilt[(i, l, 0)]) for i in I for l in L))

# (B.80) Transportation Capital Cost (TCC_t) other periods
Model.addConstrs(TCC_t[t] == gp.quicksum(TMC_il[(i, l)] * ITU_ilt[(i, l, t)] for i in I for l in L) for t in T if t != 0)

# (B.81) Fuel Cost (FC_t)
Model.addConstrs(FC_t[t] == gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilggt[(i, l, g, g_prime, t)] * (1 / FEbetween_l[l]) * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                + gp.quicksum(FP_l[l] * (2 * L_lgg[(l, g, g)] * Q_ilggt[(i, l, g, g, t)] * (1 / FEwithin_l[l]) * (1 / TCap_il[(i, l)]))
                              for i in I for l in L for g in G if A_il[(i, l)] == 1) for t in T)

# (B.82) Labour Cost (LC_t)
Model.addConstrs(LC_t[t] == gp.quicksum(DW_l[l] * (Q_ilggt[(i, l, g, g_prime, t)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g_prime)] * (1 / SPbetween_l[l]) + LUT_l[l]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1 if g_prime != g)
                + gp.quicksum(DW_l[l] * (Q_ilggt[(i, l, g, g, t)] * (1 / TCap_il[(i, l)]) * (2 * L_lgg[(l, g, g)] * (1 / SPwithin_l[l]) + LUT_l[l]))
                              for i in I for l in L for g in G if A_il[(i, l)] == 1) for t in T)

# (B.83) Maintenance Cost (MC_t)
Model.addConstrs(MC_t[t] == gp.quicksum(ME_l[l] * (2 * L_lgg[(l, g, g_prime)] * Q_ilggt[(i, l, g, g_prime, t)] * (1 / TCap_il[(i, l)]))
                                  for i in I for l in L for g in G for g_prime in G if A_il[(i, l)] == 1) for t in T)

# (B.84) General Cost (GC_t)
Model.addConstrs(GC_t[t] == gp.quicksum(GE_l[l] * NTU_ilt[(i, l, t)] for i in I for l in L) for t in T)

# (B.85) Transportation Operating Cost (TOC_t)
Model.addConstrs(TOC_t[t] == FC_t[t] + LC_t[t] + MC_t[t] + GC_t[t] for t in T)

# (B.86) Total Daily Cost (TDC_t)
Model.addConstrs(TDC_t[t] == (1 / (alpha * CCF_t[t])) * (FCC_t[t] + TCC_t[t]) + FOC_t[t] + TOC_t[t] + FSC_t[t] for t in T)

# (B.87) Total Daily Cost (TDC) averaged over the periods
Model.addConstr(TDC == gp.quicksum((CCF_t[t] / sum(CCF_t)) * TDC_t[t] for t in T))
## Alternative minimizing total costs over all the days
# Model.addConstr(TDC == gp.quicksum(TDC_t[t] * (alpha*CCF_t[t]) for t in T))

""" CO2 part """

# (B.89) Total feed emissions per time period
if newEmissionData:
    Model.addConstrs(TCEFeed_t[t] == gp.quicksum(CO2Feed_pt[(p, t)] * (gp.quicksum(P_pigt[(p, i, g, t)] for g in G for i in I)) for p in P) for t in T)
else:
    Model.addConstrs(TCEFeed_t[t] == gp.quicksum(CO2Feed_p[p] * (gp.quicksum(P_pigt[(p, i, g, t)] for g in G for i in I)) for p in P) for t in T)

# (B.90) Total production emissions per time period
if newEmissionData:
    Model.addConstrs(TCEProd_t[t] == gp.quicksum(CO2Prod_pi[(p, i)] * (gp.quicksum(P_pigt[(p, i, g, t)] for g in G)) for p in P for i in I) for t in T)
else:
    if includeCCS:
        Model.addConstrs(TCEProd_t[t] == gp.quicksum(((1 - A_CCS_p[p]) * CO2Prod_pi[(p, i)] + A_CCS_p[p] * (1 - CCSEF) * CO2Prod_pi[(p, i)]) * P_pigt[(p, i, g, t)] for p in P for i in I for g in G)
                         for t in T)
    else:
        Model.addConstrs(TCEProd_t[t] == gp.quicksum(CO2Prod_pi[(p, i)] * (gp.quicksum(P_pigt[(p, i, g, t)] for g in G)) for p in P for i in I) for t in T)

# (B.91) Total transportation emissions per time period
Model.addConstrs(TCETrans_t[t] == CO2Trans * gp.quicksum(NOT_ilggt[(i, l, g, g_prime, t)] * 2 * L_lgg[(l, g, g_prime)] for i in I for l in L for g in G for g_prime in G) for t in T)

# (B.92) Total chain emissions per time period
Model.addConstrs(TCE_t[t] == TCEFeed_t[t] + TCEProd_t[t] + TCETrans_t[t] for t in T)

# (B.93) Total chain emissions
Model.addConstr(TCE == gp.quicksum((CCF_t[t] / sum(CCF_t)) * TCE_t[t] for t in T))

# (B.94) Linear constraint on total chain emissions per period
for t in T:
    if TotalCO2Max_t[t] != -1:
        Model.addConstr(TCE_t[t] <= TotalCO2Max_t[t])

# (B.95) Linear constraint on total chain emissions per period
if TotalCO2Max != -1:
    Model.addConstr(TCE <= TotalCO2Max)

##################### Restrict Carbon intensities DURING-optimization ###########################################

# (B.104)-(B.113) Include CI during-optimization
if optimizeCI:

    # (B.104) Calculate the starting carbon intensity as parameter per plant p per product i
    pass  # already in data

    # (B.105) CI^{Plants}_{igt}
    Model.addConstrs(CIPlants_igt[(i, g, t)] * PT_igt[(i, g, t)] == gp.quicksum(P_pigt[(p, i, g, t)] * CIStart_pit[(p, i, t)] for p in P) for i in I for g in G for t in T)

    # (B.106) CI^{Prod}_{igt}
    Model.addConstrs(CIProd_igt[(i, g, t)] * gp.quicksum(Q_ilggt[(i, l_star, g_prime_star, g, t)] for l_star in L for g_prime_star in G)
                     == gp.quicksum(Q_ilggt[(i, l, g_prime, g, t)] * CIPlants_igt[(i, g_prime, t)] for l in L for g_prime in G) for i in I for g in G for t in T)

    # (B.107) CI^{Trans}_{igt}
    Model.addConstrs(CITrans_igt[(i, g, t)] * (DL_igt[(i, g, t)] + DI_igt[(i, g, t)])
                     == CO2Trans * gp.quicksum(NOT_ilggt[(i, l, g_prime, g, t)] * 2 * L_lgg[(l, g_prime, g)] for l in L for g_prime in G) for i in I for g in G for t in T)

    # (B.108) CI_{igt}
    Model.addConstrs(CI_igt[(i, g, t)] == CIProd_igt[(i, g, t)] + CITrans_igt[(i, g, t)] for i in I for g in G for t in T)

    # (B.109) CI_{gt}
    Model.addConstrs(CI_gt[(g, t)] == gp.quicksum(CI_igt[(i, g, t)] * (DL_igt[(i, g, t)] + DI_igt[(i, g, t)]) * (1 / DT_gt[(g, t)]) for i in I) for g in G for t in T if DT_gt[(g, t)] > 0)

    if not useStartSolution:
        # (B.110) CI_{ig} restrictions
        for i in I:
            for g in G:
                for t in T:
                    if maxCIDynamic[(i, g, t)] != -1:
                        Model.addConstr(CI_igt[(i, g, t)] <= maxCIDynamic[(i, g, t)])

        # (B.111) CI_{g} restrictions
        for g in G:
            for t in T:
                if maxCIDynamic[(len(I), g, t)] != -1:
                    Model.addConstr(CI_gt[(g, t)] <= maxCIDynamic[(len(I), g, t)])


# (B.88) Minimize TDC or TCE, or both
if not useStartSolution:
    # No starting solution: just minimize what you want
    objective = costImportance * TDC + CO2Importance * TCE
else:
    # Use starting solution: minimize the chain emissions first
    objective = TCE + 0

############################## Execute ########################################################################

if useStartSolution:
    """ First run: Execute minimizing the chain emissions."""
    start_time = time.time()
    Model.setObjective(objective, GRB.MINIMIZE)
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print()
    print("=== FIRST RUN: ===")
    print("Objective emission minimization:", objective.getValue(), "(total time ", round(total_time, 1), " s)")
    print("==================")
    print()

    """ Second run: Execute minimizing the costs again + restriction(s) on CI."""
    objective = TDC + 0
    # (B.110) CI_{ig} restrictions
    for i in I:
        for g in G:
            for t in T:
                if maxCIDynamic[(i, g, t)] != -1:
                    Model.addConstr(CI_igt[(i, g, t)] <= maxCIDynamic[(i, g, t)])
    # (B.111) CI_{g} restrictions
    for g in G:
        for t in T:
            if maxCIDynamic[(len(I), g, t)] != -1:
                Model.addConstr(CI_gt[(g, t)] <= maxCIDynamic[(len(I), g, t)])

start_time = time.time()
Model.setObjective(objective, GRB.MINIMIZE)
Model.optimize()
end_time = time.time()
total_time = end_time - start_time

####################### Calculate Carbon intensities POST-optimization ##########################################

""" Options: """
# Carbon_ig_post : total daily CO2 emissions for product i in grid g at time period t (per day)
# Carbon_g_post  : total daily CO2 emissions for grid g at time period t (per day)
# CI_igt_post     : carbon intensity of product i in grid g at time period t (at a day)
# CI_gt_post      : avg. carbon intensity of utilized products in grid g at time period t

# (B.96) Calculate the starting carbon intensity as parameter per plant p per product i
pass  # already in data

# (B.97) Calculate the carbon intensities of the all the plants producing product i in plant g at period t
# Note that the simplifying assumption here is that what is produced totally in a grid is distributed equally over the other grids.
# ..in other words, if you have one CG plant in a grid and one SMR plant in the same grid, you average the emissions of those grids for the
# ..outgoing products, while in fact you could connect some cities with the CG plant and some cities with the SMR plant
# ..BUT this only holds for the same product i! Depends whether you take g, (i, g) or (p, i, g) + time t
CIPlants_igt_post = np.zeros(shape=(len(I), len(G), len(T)))
for i in I:
    for g in G:
        for t in T:
            if PT_igt[(i, g, t)].X > 0.0000001:
                CIPlants_igt_post[(i, g, t)] = gp.quicksum(P_pigt[(p, i, g, t)] / PT_igt[(i, g, t)].X * CIStart_pit[(p, i, t)] for p in P).getValue()

# (B.98) Carbon intensities from production units for product i at customer location g
CIProd_igt_post = np.zeros(shape=(len(I), len(G), len(T)))
CIProd_denominator_igt = np.zeros(shape=(len(I), len(G), len(T)))
for i in I:
    for g in G:
        for t in T:
            CIProd_denominator_igt[(i, g, t)] = gp.quicksum(Q_ilggt[(i, l_star, g_prime_star, g, t)] for l_star in L for g_prime_star in G).getValue()
for i in I:
    for g in G:
        for t in T:
            if CIProd_denominator_igt[(i, g, t)] > 0.0000001:
                CIProd_igt_post[(i, g, t)] = gp.quicksum(Q_ilggt[(i, l, g_prime, g, t)] / CIProd_denominator_igt[(i, g, t)] * CIPlants_igt_post[(i, g_prime, t)] for l in L for g_prime in G).getValue()

# (B.99) Carbon intensity from transportation TODO
CITrans_igt_post = np.zeros(shape=(len(I), len(G), len(T)))
for i in I:
    for g in G:
        for t in T:
            if (DL_igt[(i, g, t)].X + DI_igt[(i, g, t)].X) > 0.0000001:
                CITrans_igt_post[(i, g, t)] = ((CO2Trans * gp.quicksum(NOT_ilggt[(i, l, g_prime, g, t)] * 2 * L_lgg[(l, g_prime, g)] for l in L for g_prime in G)) / (DL_igt[(i, g, t)].X + DI_igt[(i, g, t)].X)).getValue()

# (B.100) Total carbon intensity is part from production and part from transportation
CI_igt_post = np.zeros(shape=(len(I), len(G), len(T)))
for i in I:
    for g in G:
        for t in T:
            CI_igt_post[(i, g, t)] = CIProd_igt_post[(i, g, t)] + CITrans_igt_post[(i, g, t)]

# (B.101) Carbon intensity at a location is the resource-weighted average of both products
CI_gt_post = np.zeros(shape=(len(G), len(T)))
for g in G:
    for t in T:
        if DT_gt[(g, t)] > 0:
            CI_gt_post[(g, t)] = gp.quicksum(CI_igt_post[(i, g, t)] * (DL_igt[(i, g, t)] + DI_igt[(i, g, t)]) / DT_gt[(g, t)] for i in I).getValue()

# (B.102) Total carbon emissions for which grid g is responsible, per product
Carbon_igt_post = np.zeros(shape=(len(I), len(G), len(T)))
for i in I:
    for g in G:
        for t in T:
            Carbon_igt_post[(i, g, t)] = CI_igt_post[(i, g, t)] * (DL_igt[(i, g, t)].X + DI_igt[(i, g, t)].X)

# (B.103) Total carbon emissions for which grid g is responsible, all products combined
Carbon_gt_post = np.zeros(shape=(len(G), len(T)))
for g in G:
    for t in T:
        Carbon_gt_post[(g, t)] = CI_gt_post[(g, t)] * DT_gt[(g, t)]

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

    for t in T:
        print(f"======================================")
        print(f"======= TIME PERIOD: {TSet[t]} =======")
        print(f"======================================")
        print("=== OPERATIONAL VARIABLES ===")
        print()

        # Continuous variables
        matrix = DL_igt[:, :, t].X.round(2)
        if not deleteZeroMatrix(matrix):  # if it is a zero matrix
            print("D^L_ig (local demand satisfied): \n", pd.DataFrame(deleteZeroElements(DL_igt[:, :, t].X.round(2)), index=ISet, columns=GSet).to_string())
            print()

        matrix = DI_igt[:, :, t].X.round(2)
        if not deleteZeroMatrix(matrix):
            print("D^I_ig (imported demand satisfied): \n", pd.DataFrame(deleteZeroElements(DI_igt[:, :, t].X.round(2)), index=ISet, columns=GSet).to_string())
            print()

        for p in P:
            matrix = P_pigt[p, :, :, t].X.round(2)
            if not deleteZeroMatrix(matrix):
                print(f"P_[{PSet[p]}]ig (production from plant p={PSet[p]}): \n {pd.DataFrame(deleteZeroElements(P_pigt[p, :, :, t].X.round(2)), index=ISet, columns=GSet).to_string()}")
                print()

        matrix = PT_igt[:, :, t].X.round(2)
        if not deleteZeroMatrix(matrix):
            print("P^T_ig (total production): \n", pd.DataFrame(deleteZeroElements(PT_igt[:, :, t].X.round(2)), index=ISet, columns=GSet).to_string())
            print()

        for i in I:
            for l in L:
                if A_il[(i, l)] == 1:
                    matrix = Q_ilggt[i, l, :, :, t].X.round(2)
                    if not deleteZeroMatrix(matrix):
                        print(f"Q_[{ISet[i]}][{LSet[l]}]gg' (transport product {ISet[i]} using mode {LSet[l]}): "
                              f"\n {pd.DataFrame(deleteZeroElements(Q_ilggt[i, l, :, :, t].X.round(2)), index=GSet, columns=GSet).to_string()}")
                        print()

        # Total feedstock needed (ambiguous):
        FeedNeeded_p_large = np.zeros(shape=len(P))
        for p in P:
            FeedNeeded_p_large[p] = (PCR_p[p] * gp.quicksum(P_pigt[(p, i, g, t)] for i in I for g in G)).getValue()
        print("FeedNeeded_p \n", pd.DataFrame(deleteZeroElements(FeedNeeded_p_large), index=PSet, columns=["Amount of feed"]).to_string())
        print()

        # Integer variables
        for p in P:
            matrix = NP_pigt[p, :, :, t].X.round(2)
            if not deleteZeroMatrix(matrix):
                print(f"NP_[{PSet[p]}]ig (number of plant {PSet[p]}): \n {pd.DataFrame(deleteZeroElements(NP_pigt[p, :, :, t].X.round(2)), index=ISet, columns=GSet).to_string()}")
                print()

        for p in P:
            matrix = IP_pigt[p, :, :, t].X.round(2)
            if not deleteZeroMatrix(matrix):
                print(f"IP_[{PSet[p]}]ig (invested of plant {PSet[p]}): \n {pd.DataFrame(deleteZeroElements(IP_pigt[p, :, :, t].X.round(2)), index=ISet, columns=GSet).to_string()}")
                print()

        for i in I:
            for l in L:
                if A_il[(i, l)] == 1:
                    matrix = NOT_ilggt[i, l, :, :, t].X.round(2)
                    if not deleteZeroMatrix(matrix):
                        print(f"NOT_[{ISet[i]}][{LSet[l]}]gg (#units transport product i={ISet[i]} using mode l={LSet[l]}): \n"
                              f" {pd.DataFrame(deleteZeroElements(NOT_ilggt[i, l, :, :, t].X.round(2)), index=GSet, columns=GSet).to_string()}")
                        print()

        print("TOTAL NTU_ilt (number of units): \n", pd.DataFrame(deleteZeroElements(NTU_ilt[:, :, t].X.round(2)), index=ISet, columns=LSet).to_string())
        print()

        print("INVESTED ITU_ilt (number of units): \n", pd.DataFrame(deleteZeroElements(ITU_ilt[:, :, t].X.round(2)), index=ISet, columns=LSet).to_string())
        print()

        # Binary variables
        if printBinaryVariables:
            for i in I:
                for l in L:
                    if A_il[(i, l)] == 1:
                        matrix = X_ilggt[i, l, :, :, t].X
                        if not deleteZeroMatrix(matrix):
                            print(f"X_[{ISet[i]}][{LSet[l]}]gg' (binary for transport product {ISet[i]} using mode {LSet[l]}): "
                                  f"\n {pd.DataFrame(deleteZeroElements(X_ilggt[i, l, :, :, t].X), index=GSet, columns=GSet).to_string()}")
                            print()

            print("Y_igt (binary export): \n", pd.DataFrame(deleteZeroElements(Y_igt[:, :, t].X), index=ISet, columns=GSet).to_string())
            print()

            print("Z_igt (binary import): \n", pd.DataFrame(deleteZeroElements(Z_igt[:, :, t].X), index=ISet, columns=GSet).to_string())
            print()

def correctCI(CI_igt_input, DI_igt_input, DL_igt_input):
    """ Function to set CI values to 0 is there is no satisfied demand
    :param CI_igt_input: Carbon intensity of product i and grid g
    :param DI_igt_input: Satisfied demand imported of product i in grid g
    :param DL_igt_input: Satisfied demand locally satisfied of product in in grid g."""
    for i in I:
        for g in G:
            for t in T:
                if (DI_igt_input[(i, g, t)] + DL_igt_input[(i, g, t)] < 0.0001) and (CI_igt_input[(i, g, t)] > 0):
                    CI_igt_input[(i, g, t)] = 0
    return CI_igt_input


print("POST-optimization Carbon intensities and total carbon emissions per grid: ")
for t in T:
    print(f"======================================")
    print(f"======= TIME PERIOD: {TSet[t]} =======")
    print(f"======================================")
    CI_and_Carbon = np.zeros(shape=(len(G), 6))
    CI_and_Carbon[:, 0] = CI_igt_post[0, :, t]
    CI_and_Carbon[:, 1] = CI_igt_post[1, :, t]
    CI_and_Carbon[:, 2] = CI_gt_post[:, t]
    CI_and_Carbon[:, 3] = Carbon_igt_post[0, :, t]
    CI_and_Carbon[:, 4] = Carbon_igt_post[1, :, t]
    CI_and_Carbon[:, 5] = Carbon_gt_post[:, t]
    print(pd.DataFrame(CI_and_Carbon, index=GSet, columns=["CI_[CH2]g", "CI_[LH2]g", "CI_gt", "Carbon_[CH2]g", "Carbon_[LH2]g", "Carbon_g"]).to_string())
    print()

    if optimizeCI:

        CI_igt_corr = correctCI(CI_igt.X, DI_igt.X, DL_igt.X)

        print("DURING-optimization Carbon intensities and total carbon emissions per grid: ")
        CI_and_Carbon = np.zeros(shape=(len(G), 6))
        CI_and_Carbon[:, 0] = CI_igt_corr[0, :, t]
        CI_and_Carbon[:, 1] = CI_igt_corr[1, :, t]
        CI_and_Carbon[:, 2] = CI_gt.X[:, t]
        CI_and_Carbon[:, 3] = CI_igt.X[0, :, t] * (DL_igt.X[0, :, t] + DI_igt.X[0, :, t])
        CI_and_Carbon[:, 4] = CI_igt.X[1, :, t] * (DL_igt.X[1, :, t] + DI_igt.X[1, :, t])
        CI_and_Carbon[:, 5] = CI_gt.X[:, t] * DT_gt[:, t]
        print(pd.DataFrame(CI_and_Carbon, index=GSet, columns=["CI_[CH2]g", "CI_[LH2]g", "CI_gt", "Carbon_[CH2]g", "Carbon_[LH2]g", "Carbon_g"]).to_string())
        print()

    print("Validity check: ")
    print("TCE_t                        =", round(TCE_t[t].X, 3))
    print("POST:   sum(Carbon_g_post)_t =", round(sum(Carbon_gt_post[:, t]), 3))
    if optimizeCI:
        print("DURING: sum(Carbon_g)_t  =", round(sum(CI_gt[:, t].X * DT_gt[:, t]), 3))
    print()

print()
print("===== RESULT =====")
print(f"Objective  : {round(objective.getValue(), 3)} (costImportance={costImportance}, CO2Importance={CO2Importance})")
print(f"Total time : {round(total_time, 3)} sec")
print()

for t in T:
    print(f"======================================")
    print(f"======= TIME PERIOD: {TSet[t]} ==============")
    print(f"======================================")
    FCC_t_daily = (FCC_t[t].X / (alpha * CCF_t[t]))
    TCC_t_daily = (TCC_t[t].X / (alpha * CCF_t[t]))
    print("=== COST VARIABLES ===")
    print(f"TDC_t (Total Daily Cost)            =  {round(TDC_t[t].X, 2)}  $/day")
    print(f" * FSC_t (FeedStock Cost)           =  {round(FSC_t[t].X, 2)}  $/day ({round(100 * (FSC_t[t].X / TDC_t[t].X), 1)}%)")
    print(f" * FCC_t (Facility Capital Cost)    =  {round(FCC_t_daily, 2)}  $/day ({round(100 * (FCC_t_daily / TDC_t[t].X), 1)}%)  ({round(FCC_t[t].X, 2)} $ in total)")
    print(f" * FOC_t (Facility Operating Cost)  =  {round(FOC_t[t].X, 2)}  $/day ({round(100 * (FOC_t[t].X / TDC_t[t].X), 1)}%)")
    print(f" * TCC_t (Transport Capital Cost)   =  {round(TCC_t_daily, 2)}  $/day ({round(100 * (TCC_t_daily / TDC_t[t].X), 1)}%)  ({round(TCC_t[t].X, 2)} $ in total)")
    print(f" * TOC_t (Transport Operating Cost) =  {round(TOC_t[t].X, 2)}  $/day ({round(100 * (TOC_t[t].X / TDC_t[t].X), 1)}%)")
    print(f"    - FC_t (Fueling Cost)           =  {round(FC_t[t].X, 2)}  $/day ({round(100 * (FC_t[t].X / TOC_t[t].X), 1)}% of TOC_t), ({round(100 * (FC_t[t].X / TDC_t[t].X), 1)}% of total)")
    print(f"    - LC_t (Labour Cost)            =  {round(LC_t[t].X, 2)}  $/day ({round(100 * (LC_t[t].X / TOC_t[t].X), 1)}% of TOC_t), ({round(100 * (LC_t[t].X / TDC_t[t].X), 1)}% of total)")
    print(f"    - MC_t (Maintenance Cost)       =  {round(MC_t[t].X, 2)}  $/day ({round(100 * (MC_t[t].X / TOC_t[t].X), 1)}% of TOC_t), ({round(100 * (MC_t[t].X / TDC_t[t].X), 1)}% of total)")
    print(f"    - GC_t (General Cost)           =  {round(GC_t[t].X, 2)}  $/day ({round(100 * (GC_t[t].X / TOC_t[t].X), 1)}% of TOC_t), ({round(100 * (GC_t[t].X / TDC_t[t].X), 1)}% of total)")
    print()

    print("=== CARBON VARIABLES ===")
    print(f"TCE_t (Total Carbon Emission)  =  {round(TCE_t[t].X, 2)} ton CO2/day")
    print(f"  * TCEFeed_t (Feedstock Em.)  =  {round(TCEFeed_t[t].X, 1)}  ton CO2/day ({round(100 * (TCEFeed_t[t].X / TCE_t[t].X), 1)}%)")
    print(f"  * TCEProd_t (Production Em.) =  {round(TCEProd_t[t].X, 2)}  ton CO2/day ({round(100 * (TCEProd_t[t].X / TCE_t[t].X), 1)}%)")
    print(f"  * TCETrans_t (Transport Em.) =  {round(TCETrans_t[t].X, 2)} ton CO2/day ({round(100 * (TCETrans_t[t].X / TCE_t[t].X), 1)}%)")
    print()
print("Averaged over all time periods:")
print(f"TDC =  {round(TDC.X, 2)}  $/day")
print(f"TCE =  {round(TCE.X, 2)}  $/day")

# if printTests:
#     print("==================== TEST =============================")
#
#     print("CIStart_pi: \n", pd.DataFrame(CIStart_pi, index=PSet, columns=ISet).to_string())
#     print()
#
#     print("CIPlants_igt_post: \n", pd.DataFrame(CIPlants_igt_post, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("CIProd_denominator_ig: \n", pd.DataFrame(CIProd_denominator_ig, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("CIProd_igt_post: \n", pd.DataFrame(CIProd_igt_post, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("CITrans_igt_post: \n", pd.DataFrame(CITrans_igt_post, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("CI_igt_post: \n", pd.DataFrame(CI_igt_post, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("C_g: \n", pd.DataFrame(CI_gt_post, index=GSet).T.to_string())
#     print()
#
#     print("Carbon_ig_post: \n", pd.DataFrame(Carbon_ig_post, index=ISet, columns=GSet).to_string())
#     print()
#
#     print("Carbon_g_post: \n", pd.DataFrame(Carbon_g_post, index=GSet).T.to_string())
#     print()
