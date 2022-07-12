""" Data_MRMT

The topology and data of the experimental setup of the Multi-Resource, Multi-Transport (MRMT) model as presented in the thesis.

The data contains elements from the following sets:
# R = Resources (products/materials);
# S = Source nodes;
# P = Production nodes;
# D = Depot nodes;
# T = Target nodes;
# M = Modes of transport;
# P' = Subset of production nodes where CCS is defined on;
# L = levels of invest for CCS.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

Note: a data entry with -1 means +infinity.
"""

# Imports
import gurobipy as gp
import numpy as np
import pandas as pd

""" Nodes """

# Resources
resources = ["R1", "R2", "R3", "R4", "R5", "R6"]

# Sources
sources = ["S1", "S2", "S3", "S4"]

sources_res, cost_source, min_supply_source, max_supply_source, carbon_source = gp.multidict({
    ("S1", "R1"): [10, 0, 1000, 20],
    ("S2", "R2"): [15, 0, 1000, 15],
    ("S3", "R3"): [2, 0, 1000, 10],
    ("S3", "R4"): [2, 0, 1000, 10],
    ("S4", "R3"): [5, 20, 1000, 2],
    ("S4", "R4"): [5, 20, 1000, 2]
})

# Production units (short: prods)
prods, cost_prod, max_through_prod, carbon_prod = gp.multidict({
    "P1": [10, -1, 20],
    "P2": [12, 500, 15]
})

prods_res = [("P1", "R1"), ("P1", "R3"), ("P1", "R4"), ("P1", "R5"), ("P2", "R2"), ("P2", "R3"), ("P2", "R4"), ("P2", "R6")]

# Yield table (production units as rows, resources as columns)
yields_array = np.array([
    #  R1  R2    R3     R4  R5 R6
    [-0.65, 0, -0.25, -0.10, 1, 0],  # P1
    [0, -0.70, -0.15, -0.15, 0, 1]   # P2
])
yields = pd.DataFrame(yields_array, columns=resources, index=prods)

# Depots
depots, min_through_depot, fixed_cost_depot, fixed_carbon_depot = gp.multidict({
    "D1": [0, 500, 1500],
    "D2": [50, 750, 450]
})

depots_res, max_through_depot = gp.multidict({
    ("D1", "R5"): [-1],
    ("D1", "R6"): [-1],
    ("D2", "R5"): [300],
    ("D2", "R6"): [300]
})

# Targets
targets = ["T1", "T2", "T3", "T4"]

targets_res, price_target, min_demand_target, max_demand_target, carbon_target, max_CI_target = gp.multidict({
    ("T1", "R5"): [40, 70, 80, 0, -1],
    ("T1", "R6"): [50, 80, 90, 0, -1],
    ("T2", "R5"): [40, 85, 95, 0, -1],
    ("T2", "R6"): [50, 65, 75, 0, -1],
    ("T3", "R5"): [40, 70, 80, 0, -1],
    ("T4", "R6"): [50, 85, 95, 0, -1],
})

# Subset of prods on which CCS is defined
prods_CCS = ["P1", "P2"]

prods_levels_CCS, fixed_cost_install_CCS, max_capture_CCS = gp.multidict({
    ("P1", "L1"): [1000, 2000],  # can capture 2 units of CO2 with 1 dollar
    ("P1", "L2"): [2000, 5000],  # can capture 2.5 units of CO2 with 1 dollar
    ("P1", "L3"): [3000, 9000],  # can capture 3 units of CO2 with 1 dollar
    ("P2", "L1"): [1000, 2000],
    ("P2", "L2"): [2000, 5000]
})

""" Edges """

modes = ["M1", "M2", "M3", "M4"]

# Arc = from node i to node j using mode m for resource r
arcs, cost_arc, carbon_arc, max_through_arc = gp.multidict({
    # Sources-to-pools
    ("S1", "P1", "M1", "R1"): [1, 1, -1],
    ("S2", "P2", "M1", "R2"): [1, 1, -1],
    ("S3", "P1", "M1", "R3"): [1, 1, -1],
    ("S3", "P1", "M1", "R4"): [1, 1, -1],
    ("S3", "P2", "M1", "R3"): [1, 1, -1],
    ("S3", "P2", "M1", "R4"): [1, 1, -1],
    ("S4", "P1", "M1", "R3"): [1, 1, -1],
    ("S4", "P1", "M1", "R4"): [1, 1, -1],
    ("S4", "P2", "M1", "R3"): [1, 1, -1],
    ("S4", "P2", "M1", "R4"): [1, 1, -1],
    # Pools-to-depots
    ("P1", "D1", "M2", "R5"): [1, 3, -1],
    ("P1", "D1", "M3", "R5"): [2, 1, 200],
    ("P1", "D2", "M2", "R5"): [1, 3, -1],
    ("P1", "D2", "M3", "R5"): [2, 1, 200],
    ("P2", "D1", "M2", "R6"): [1, 3, -1],
    ("P2", "D1", "M3", "R6"): [2, 1, 300],
    ("P2", "D2", "M2", "R6"): [1, 3, -1],
    ("P2", "D2", "M3", "R6"): [2, 1, 300],
    # Depots-to-targets
    ("D1", "T1", "M4", "R5"): [1, 1, -1],
    ("D1", "T1", "M4", "R6"): [1, 1, -1],
    ("D1", "T2", "M4", "R5"): [1, 1, -1],
    ("D1", "T2", "M4", "R6"): [1, 1, -1],
    ("D1", "T3", "M4", "R5"): [1, 1, -1],
    ("D1", "T4", "M4", "R6"): [1, 1, -1],
    ("D2", "T1", "M4", "R5"): [1, 1, -1],
    ("D2", "T1", "M4", "R6"): [1, 1, -1],
    ("D2", "T2", "M4", "R5"): [1, 1, -1],
    ("D2", "T2", "M4", "R6"): [1, 1, -1],
    ("D2", "T3", "M4", "R5"): [1, 1, -1],
    ("D2", "T4", "M4", "R6"): [1, 1, -1]
})
