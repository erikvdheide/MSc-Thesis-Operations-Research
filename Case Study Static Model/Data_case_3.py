""" Data_2

NEWLY ADDED: CCS plants.

The data for the case study model.

The data contains the following sets:
# G = Grids/locations;
# I = Products physical form;
# L = Modes of transport;
# P = Plant types/technologies;

Thesis OR&QL - EUR 2022
@author: Erik van der Heide
"""

# Imports
import numpy as np
from math import sin, cos, sqrt, atan2, radians
from Params_case import *
import pandas as pd

""" Sets """

# The PSet basically has the attribute plant type (p), size (j) and CCS (c), but, for the dynamic case we can pool those indices together in one index p.
if multipleSizes and includeCCS:
    PSet = ["SMR-Small", "SMR-Medium", "SMR-Large", "CG-Small", "CG-Medium", "CG-Large", "BG-Small", "BG-Medium", "BG-Large", "WE-Small", "WE-Medium", "WE-Large",
            "SMR-Small-CCS", "SMR-Medium-CCS", "SMR-Large-CCS", "CG-Small-CCS", "CG-Medium-CCS", "CG-Large-CCS", "BG-Small-CCS", "BG-Medium-CCS", "BG-Large-CCS",
            "WE-Small-CCS", "WE-Medium-CCS", "WE-Large-CCS"]  # 24 elements
elif multipleSizes:
    PSet = ["SMR-Small", "SMR-Medium", "SMR-Large", "CG-Small", "CG-Medium", "CG-Large", "BG-Small", "BG-Medium", "BG-Large", "WE-Small", "WE-Medium", "WE-Large"]  # 12 elements
elif includeCCS:
    PSet = ["SMR", "CG", "BG", "WE", "SMR-CCS", "CG-CCS", "BG-CCS", "WE-CCS"]  # 8 elements
else:
    PSet = ["SMR", "CG", "BG", "WE"]  # 4 elements
P = list(range(len(PSet)))

ISet = ["CH2", "LH2"]
I = list(range(len(ISet)))

LSet = ["tube", "tanker"]
L = list(range(len(LSet)))

GSet = ["G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12", "G13", "G14", "G15", "G16", "G17", "G18", "G19", "G20", "G21", "G22", "G23", "G24", "G25"]
G = list(range(len(GSet)))

""" Incident matrices """

# Product-mode matrix
A_il = np.identity(n=len(I))

# Transport-between-grids matrix
A_gg = np.ones(shape=(len(G), len(G)))
for g in range(0, len(G)):
    A_gg[(g, g)] = 0

# Locations to build new plants (technically not incident but data matrix)
A_ig = np.zeros(shape=(len(I), len(G)))
A_ig[(0, 0)] = A_ig[(0, 4)] = A_ig[(0, 9)] = A_ig[(0, 17)] = A_ig[(0, 22)] = 1  # can build new CH2 locations at G01, G05, G10, G18, G23
A_ig[(1, 0)] = 1  # can build new LH2 locations at G01

# If it is a CCS plant
if includeCCS:
    halfway = int(len(P)/2)
    A_CCS_p = np.zeros(shape=len(P))
    for r in range(halfway, len(P)):
        A_CCS_p[r] = 1

""" Parameters """

if includeCCS:
    CCSEF = 0.9
    CCSCost = CCSCostInput

CO2Prod_pi = np.zeros(shape=(len(P), len(I)))  # ton CO2/ton H2
if not multipleSizes:
    # part without CCS, take medium values (DO NOT CHANGE!)
    CO2Prod_pi[(0, 0)] = 10.8
    CO2Prod_pi[(0, 1)] = 14.6
    CO2Prod_pi[(1, 0)] = 23.8
    CO2Prod_pi[(1, 1)] = 26.6
    CO2Prod_pi[(2, 0)] = 25.4
    CO2Prod_pi[(2, 1)] = 27.5
    CO2Prod_pi[(3, 0)] = 41.1
    CO2Prod_pi[(3, 1)] = 24.5
else:
    # SMR
    CO2Prod_pi[(0, 0)] = 11.4
    CO2Prod_pi[(0, 1)] = 15.1
    CO2Prod_pi[(1, 0)] = 10.8
    CO2Prod_pi[(1, 1)] = 14.6
    CO2Prod_pi[(2, 0)] = 10.3
    CO2Prod_pi[(2, 1)] = 14.0
    # CG
    CO2Prod_pi[(3, 0)] = 24.8
    CO2Prod_pi[(3, 1)] = 27.6
    CO2Prod_pi[(4, 0)] = 23.8
    CO2Prod_pi[(4, 1)] = 26.6
    CO2Prod_pi[(5, 0)] = 22.9
    CO2Prod_pi[(5, 1)] = 25.7
    # BG
    CO2Prod_pi[(6, 0)] = 26.9
    CO2Prod_pi[(6, 1)] = 29.0
    CO2Prod_pi[(7, 0)] = 25.4
    CO2Prod_pi[(7, 1)] = 27.5
    CO2Prod_pi[(8, 0)] = 24.0
    CO2Prod_pi[(8, 1)] = 26.2
    # WE
    CO2Prod_pi[(9, 0)] = 43.9
    CO2Prod_pi[(9, 1)] = 26.2
    CO2Prod_pi[(10, 0)] = 41.1
    CO2Prod_pi[(10, 1)] = 24.5
    CO2Prod_pi[(11, 0)] = 38.6
    CO2Prod_pi[(11, 1)] = 23.0
if includeCCS:
    # with CCS (capture corrections are done in the model):
    halfway = int(len(P)/2)
    for r in range(0, halfway):
        for i in I:
            CO2Prod_pi[(r+halfway, i)] = CO2Prod_pi[(r, i)]

if multipleSizes and includeCCS:
    # Multiple sizes + CCS (24 elem)
    CO2Feed_p = [0.58, 0.58, 0.58, 1.31, 1.31, 1.31, 0.21, 0.21, 0.21, 0, 0, 0, 0.58, 0.58, 0.58, 1.31, 1.31, 1.31, 0.21, 0.21, 0.21, 0, 0, 0]  # ton CO2/ton H2
elif multipleSizes:
    # Multiple sizes, no CCS (12 elem)
    CO2Feed_p = [0.58, 0.58, 0.58, 1.31, 1.31, 1.31, 0.21, 0.21, 0.21, 0, 0, 0]  # ton CO2/ton H2
elif includeCCS:
    # 1 size, Include CCS (8 elem)
    CO2Feed_p = [0.58, 1.31, 0.21, 0, 0.58, 1.31, 0.21, 0]  # ton CO2/ton H2
else:
    # 1 size, no CCS
    CO2Feed_p = [0.58, 1.31, 0.21, 0]  # ton CO2/ton H2

CO2Trans = 0.00075  # ton CO2/km driven (750 g/km = 0.75 kg/km = 0.00075 ton/km)

if solveTimePeriod == 1:
    CCF = 6
elif solveTimePeriod == 2 or solveTimePeriod == 3 or solveTimePeriod == 4:
    CCF = 10

# Demand in period T1, T2, T3 or T4:
if solveTimePeriod == 1:
    DT_g = [18.81921533, 3.077355966, 3.077355966, 3.077355966, 8.285189139, 5.681272552, 2.367196897, 4.497674104, 4.497674104, 0, 0, 0, 0, 0, 0, 0, 3.077355966, 0, 0, 0, 0, 0, 0, 0, 0]
elif solveTimePeriod == 2:
    DT_g = [37.62306307, 6.152198992, 6.152198992, 6.152198992, 16.56361267, 11.35790583, 4.732460763, 8.99167545, 8.99167545, 5.678952916, 5.08739532, 5.08739532, 4.732460763,
            4.732460763, 4.732460763, 4.732460763, 6.152198992, 3.076099496, 8.163494817, 3.785968611, 3.667657092, 3.549345572, 4.140903168, 4.495837725, 3.194411015]
elif solveTimePeriod == 3:
    DT_g = [134.5458014, 22.00120022, 22.00120022, 22.00120022, 59.2340006, 40.61760041, 16.92400017, 32.15560033, 32.15560033, 20.30880021, 18.19330019, 18.19330019, 16.92400017,
            16.92400017, 16.92400017, 16.92400017, 22.00120022, 11.00060011, 29.1939003, 13.53920014, 13.11610013, 12.69300013, 14.80850015, 16.07780016, 11.42370012]
elif solveTimePeriod == 4:
    DT_g = [334.5793782, 54.7110933, 54.7110933, 54.7110933, 147.2990973, 101.0050953, 42.08545638, 79.96236713, 79.96236713, 75.75382149, 67.86279842, 67.86279842, 63.12818457,
            63.12818457, 63.12818457, 63.12818457, 54.7110933, 41.03331997, 108.8961184, 50.50254766, 48.92434305, 47.34613843, 55.2371615, 59.97177535, 42.61152459]
total_demand = sum(DT_g)

DW_l = [35, 35]

FEwithin_l = [2.30, 2.30]
FEbetween_l = [2.55, 2.55]

if multipleSizes and includeCCS:
    # Multiple sizes + CCS (24 elem)
    FSP_p = [120, 120, 120, 30, 30, 30, 50, 50, 50, 50, 50, 50, 120, 120, 120, 30, 30, 30, 50, 50, 50, 50, 50, 50]
elif multipleSizes:
    # Multiple sizes, no CCS (12 elem)
    FSP_p = [120, 120, 120, 30, 30, 30, 50, 50, 50, 50, 50, 50, 50]
elif includeCCS:
    # 1 size, Include CCS (8 elem)
    FSP_p = [120, 30, 50, 50, 120, 30, 50, 50]
else:
    # 1 size, no CCS
    FSP_p = [120, 30, 50, 50]

FP_l = [1.3, 1.3]

GE_l = [8.22, 8.22]

lat_values = [51.924419, 52.06052,	52.161091, 52.070499, 52.370216, 52.387386, 52.456955, 52.090736, 52.156113, 51.985104,	52.211159, 51.812565, 51.441643, 51.697815, 51.585255, 51.571915,
              51.81155, 51.498795,	50.851368, 52.221539, 52.516773, 52.992752,	53.219383, 53.201233, 52.350784]
long_values = [4.477733, 4.49326, 4.49015, 4.3007, 4.895168, 4.646219, 4.606014, 5.12142, 5.387827, 5.89873, 5.969923, 5.837226, 5.469722, 5.303675, 5.056375, 4.768323, 4.66636, 3.610998,
               5.690972, 6.893662, 6.083022, 6.564228, 6.566502, 5.799913, 5.264702]


# Function to calculate the distance matrix based on latitude and longitude
def distanceMatrix(lat_val, long_val):
    R = 6371.0  # approximate radius of earth in km
    distanceMat = np.zeros((len(lat_val), len(lat_val)), dtype=float)

    for i in range(0, len(lat_val)):
        for j in range(0, len(long_val)):
            lat1 = radians(lat_val[i])
            lon1 = radians(long_val[i])
            lat2 = radians(lat_val[j])
            lon2 = radians(long_val[j])
            dlon = lon2 - lon1
            dlat = lat2 - lat1
            a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
            c = 2 * atan2(sqrt(a), sqrt(1 - a))
            distanceMat[(i, j)] = R * c

    return distanceMat


L_lgg = np.zeros(shape=(len(L), len(G), len(G)))
dist_matrix = distanceMatrix(lat_values, long_values)
for g in G:
    dist_matrix[(g, g)] = distanceWithinGrid
L_lgg[0, :, :] = dist_matrix
L_lgg[1, :, :] = dist_matrix

LUT_l = [4, 2]

ME_l = [0.0976, 0.0976]

PCapMin_pi = np.zeros(shape=(len(P), len(I)))
PCapMax_pi = np.zeros(shape=(len(P), len(I)))
if not multipleSizes:
    # For only 1 size, just take between 0 and 480 tpd
    for p in P:
        for i in I:
            PCapMin_pi[(p, i)] = 10
            PCapMax_pi[(p, i)] = 480
else:
    # For multiple sizes, do it manually
    smallIndices = [0, 3, 6, 9]
    mediumIndices = [1, 4, 7, 10]
    largeIndices = [2, 5, 8, 11]
    if includeCCS:
        smallIndices.append([12, 15, 18, 21])
        mediumIndices.append([13, 16, 19, 22])
        largeIndices.append([14, 17, 20, 23])
    for i in I:
        for p in smallIndices:
            PCapMin_pi[(p, i)] = 20
            PCapMax_pi[(p, i)] = 99
        for p in mediumIndices:
            PCapMin_pi[(p, i)] = 100
            PCapMax_pi[(p, i)] = 499
        for p in largeIndices:
            PCapMin_pi[(p, i)] = 500
            PCapMax_pi[(p, i)] = 1000


# TODO UPC_pi (met goede opzet)
PCC_pi = np.zeros(shape=(len(P), len(I)))
if not multipleSizes:
    # with only 1 size, take approximations from snapshot model
    PCC_pi[(0, 0)] = 379e6
    PCC_pi[(0, 1)] = 535e6
    PCC_pi[(1, 0)] = 771e6
    PCC_pi[(1, 1)] = 958e6
    PCC_pi[(2, 0)] = 907e6
    PCC_pi[(2, 1)] = 1412e6
    PCC_pi[(3, 0)] = 1000e6
    PCC_pi[(3, 1)] = 1500e6
else:
    # SMR units:
    PCC_pi[(0, 0)] = (134 * (avgCapSmall/150)**0.75) * 10e6
    PCC_pi[(0, 1)] = (230 * (avgCapSmall/150)**0.73) * 10e6
    PCC_pi[(1, 0)] = (134 * (avgCapMedium/150)**0.75) * 10e6
    PCC_pi[(1, 1)] = (230 * (avgCapMedium/150)**0.73) * 10e6
    PCC_pi[(2, 0)] = (134 * (avgCapLarge/150)**0.75) * 10e6
    PCC_pi[(2, 1)] = (230 * (avgCapLarge/150)**0.73) * 10e6
    # CG units:
    PCC_pi[(3, 0)] = (352 * (avgCapSmall/150)**0.77) * 10e6
    PCC_pi[(3, 1)] = (460 * (avgCapSmall/150)**0.76) * 10e6
    PCC_pi[(4, 0)] = (352 * (avgCapMedium/150)**0.77) * 10e6
    PCC_pi[(4, 1)] = (460 * (avgCapMedium/150)**0.76) * 10e6
    PCC_pi[(5, 0)] = (352 * (avgCapLarge/150)**0.77) * 10e6
    PCC_pi[(5, 1)] = (460 * (avgCapLarge/150)**0.76) * 10e6
    # BG units:
    PCC_pi[(6, 0)] = (360 * (avgCapSmall/150)**0.75) * 10e6
    PCC_pi[(6, 1)] = (450 * (avgCapSmall/150)**0.72) * 10e6
    PCC_pi[(7, 0)] = (360 * (avgCapMedium/150)**0.75) * 10e6
    PCC_pi[(7, 1)] = (450 * (avgCapMedium/150)**0.72) * 10e6
    PCC_pi[(8, 0)] = (360 * (avgCapLarge/150)**0.75) * 10e6
    PCC_pi[(8, 1)] = (450 * (avgCapLarge/150)**0.72) * 10e6
    # WE units:
    PCC_pi[(9, 0)] = (598 * (avgCapSmall/150)**0.85) * 10e6
    PCC_pi[(9, 1)] = (685 * (avgCapSmall/150)**0.84) * 10e6
    PCC_pi[(10, 0)] = (598 * (avgCapMedium/150)**0.85) * 10e6
    PCC_pi[(10, 1)] = (685 * (avgCapMedium/150)**0.84) * 10e6
    PCC_pi[(11, 0)] = (598 * (avgCapLarge/150)**0.85) * 10e6
    PCC_pi[(11, 1)] = (685 * (avgCapLarge/150)**0.84) * 10e6
if includeCCS:
    # with CCS (capture corrections are done in the model):
    halfway = int(len(P)/2)
    for r in range(0, halfway):
        for i in I:
            PCC_pi[(r+halfway, i)] = PCC_pi[(r, i)]

if multipleSizes and includeCCS:
    # Multiple sizes + CCS (24 elem)
    PCR_p = [4.02, 3.34, 3.16, 5.95, 5.64, 5.33, 25.5, 18.4, 11.3, 52.5, 52.5, 52.5, 4.02, 3.34, 3.16, 5.95, 5.64, 5.33, 25.5, 18.4, 11.3, 52.5, 52.5, 52.5]
elif multipleSizes:
    # Multiple sizes, no CCS (12 elem)
    PCR_p = [4.02, 3.34, 3.16, 5.95, 5.64, 5.33, 25.5, 18.4, 11.3, 52.5, 52.5, 52.5]
elif includeCCS:
    # 1 size, Include CCS (8 elem)
    PCR_p = [3.34, 5.64, 18.4, 52.5, 3.34, 5.64, 18.4, 52.5]
else:
    # 1 size, no CCS
    PCR_p = [3.34, 5.64, 18.4, 52.5]

UPC_pi = np.zeros(shape=(len(P), len(I)))
if not multipleSizes:
    # with only 1 size, take approximations from snapshot model
    UPC_pi[(0, 0)] = 0.94e3
    UPC_pi[(0, 1)] = 1.53e3
    UPC_pi[(1, 0)] = 1.06e3
    UPC_pi[(1, 1)] = 1.71e3
    UPC_pi[(2, 0)] = 1.71e3
    UPC_pi[(2, 1)] = 3.08e3
    UPC_pi[(3, 0)] = 4.00e3
    UPC_pi[(3, 1)] = 5.00e3
else:
    # SMR units:
    UPC_pi[(0, 0)] = 3.36e3
    UPC_pi[(0, 1)] = 3.36e3
    UPC_pi[(1, 0)] = 1.74e3
    UPC_pi[(1, 1)] = 1.74e3
    UPC_pi[(2, 0)] = 1.43e3
    UPC_pi[(2, 1)] = 1.43e3
    # CG units:
    UPC_pi[(3, 0)] = 2.50e3
    UPC_pi[(3, 1)] = 2.50e3
    UPC_pi[(4, 0)] = 2.03e3
    UPC_pi[(4, 1)] = 2.03e3
    UPC_pi[(5, 0)] = 1.56e3
    UPC_pi[(5, 1)] = 1.56e3
    # BG units:
    UPC_pi[(6, 0)] = 4.89e3
    UPC_pi[(6, 1)] = 4.89e3
    UPC_pi[(7, 0)] = 3.52e3
    UPC_pi[(7, 1)] = 3.52e3
    UPC_pi[(8, 0)] = 2.15e3
    UPC_pi[(8, 1)] = 2.15e3
    # WE units:
    UPC_pi[(9, 0)] = 6.82e3
    UPC_pi[(9, 1)] = 6.82e3
    UPC_pi[(10, 0)] = 5.03e3
    UPC_pi[(10, 1)] = 5.03e3
    UPC_pi[(11, 0)] = 3.24e3
    UPC_pi[(11, 1)] = 3.24e3
if includeCCS:
    # with CCS (capture corrections are done in the model):
    halfway = int(len(P)/2)
    for r in range(0, halfway):
        for i in I:
            UPC_pi[(r+halfway, i)] = UPC_pi[(r, i)]

QMin_il = np.zeros(shape=(len(I), len(L)))
QMax_il = np.zeros(shape=(len(I), len(L)))
for i in I:
    for l in L:
        QMax_il[(i, l)] = total_demand

SPwithin_l = [25, 25]
SPbetween_l = [50, 50]

TCap_il = np.zeros(shape=(len(I), len(L)))
TCap_il[(0, 0)] = 200e-3
TCap_il[(1, 1)] = 4000e-3

TMA_l = [24, 24]

TMC_il = np.zeros(shape=(len(I), len(L)))
TMC_il[(0, 0)] = 300000
TMC_il[(1, 1)] = 800000

alpha = 365

""" Data resulting from other data """

CIStart_pi = np.zeros(shape=(len(P), len(I)))  # ton CO2/ton H2
if includeCCS:
    for p in P:
        for i in I:
            CIStart_pi[(p, i)] = CO2Feed_p[p] + (1 - A_CCS_p[p]) * CO2Prod_pi[(p, i)] + A_CCS_p[p] * (1 - CCSEF) * CO2Prod_pi[(p, i)]
else:
    for p in P:
        for i in I:
            CIStart_pi[(p, i)] = CO2Feed_p[p] + CO2Prod_pi[(p, i)]

""" Print parameters used """

if printData:

    print("===== DATA PRINT =====")
    print()
    print("===Sets: ===PC")
    print("P:", PSet)
    print("I:", ISet)
    print("L:", LSet)
    print("G:", GSet)
    print()
    print("===Incident matrices: ===")
    print("A_il: \n", pd.DataFrame(A_il, index=ISet, columns=LSet).to_string())
    print()
    print("A_gg: \n", pd.DataFrame(A_gg, index=GSet, columns=GSet).to_string())
    print()
    print("A_ig (transposed): \n", pd.DataFrame(A_ig, index=ISet, columns=GSet).T.to_string())
    print()
    if includeCCS:
        print("A_CCS_p: \n", pd.DataFrame(A_CCS_p, index=PSet).to_string())
        print()

    print("===Demand and distance data: ===")
    print("D^T_g:", DT_g)
    print()
    print("Total demand: ", total_demand)
    print()
    print("L_[tube trailer]gg':\n", pd.DataFrame(L_lgg[0, :, :], index=GSet, columns=GSet))
    print()
    print("L_[tanker truck]gg':\n", pd.DataFrame(L_lgg[1, :, :], index=GSet, columns=GSet))
    print()

    print("===Transport modes data: ===")
    print("DW_l        : ", DW_l)
    print("FEwithin_l  : ", FEwithin_l)
    print("FEbetween_l : ", FEbetween_l)
    print("FP_l        : ", FP_l)
    print("GE_l        : ", GE_l)
    print("LUT_l       : ", LUT_l)
    print("ME_l        : ", ME_l)
    print("SPwithin_l  : ", SPwithin_l)
    print("SPbetween_l : ", SPbetween_l)
    print("TMA_l       : ", TMA_l)
    print()

    print("===Plant-product data: ===")
    print("PCapMin_pi: \n", pd.DataFrame(PCapMin_pi, index=PSet, columns=ISet).to_string())
    print()
    print("PCapMax_pi: \n", pd.DataFrame(PCapMax_pi, index=PSet, columns=ISet).to_string())
    print()
    print("PCC_pi: \n", pd.DataFrame(PCC_pi, index=PSet, columns=ISet).to_string())
    print()
    print("UPC_pi: \n", pd.DataFrame(UPC_pi, index=PSet, columns=ISet).to_string())
    print()

    print("===Product-mode data: ===")
    print("QMin_il: \n", pd.DataFrame(QMin_il, index=ISet, columns=LSet).to_string())
    print()
    print("QMax_il: \n", pd.DataFrame(QMax_il, index=ISet, columns=LSet).to_string())
    print()
    print("TCap_il: \n", pd.DataFrame(TCap_il, index=ISet, columns=LSet).to_string())
    print()
    print("TMC_il: \n", pd.DataFrame(TMC_il, index=ISet, columns=LSet).to_string())
    print()

    print("===Plant data: ===")
    print("FSP_p: \n", pd.DataFrame(FSP_p, index=PSet, columns=["Feed cost"]).to_string())
    print()
    print("PCR_p: \n", pd.DataFrame(PCR_p, index=PSet, columns=["Conversion rate"]).to_string())
    print()

    print("===CO2 data: ===")
    print("CO2Prod_pi: \n", pd.DataFrame(CO2Prod_pi, index=PSet, columns=ISet).to_string())
    print()
    print("CO2Feed_p: \n", pd.DataFrame(CO2Feed_p, index=PSet, columns=["CO2 feedstock"]).to_string())
    print()
    print("CO2Trans: ", CO2Trans)
    print()
    print("CIStart_pi: \n", pd.DataFrame(CIStart_pi, index=PSet, columns=ISet).to_string())
    print()

    print("Other data: ")
    print("CCF: ", CCF)
    print("alpha: ", alpha)
    print()

    if includeCCS:
        print("CCS data: ")
        print("CCSEF: ", CCSEF)
        print("CCSCost: ", CCSCost)
        print()
