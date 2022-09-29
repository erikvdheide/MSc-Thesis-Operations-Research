""" Data_SRMT_3

Data for set-up 3: 4 pools (in 2 layers)
NOTE: -1 = infinity.

@author: Erik van der Heide

"""

import gurobipy as gp

""" Nodes """

sources, cost_source, min_supply_source, max_supply_source, CO2_source = gp.multidict({
    "S1": [10, 10, 100, 5],
    "S2": [15, 10, 120, 4],
    "S3": [1, 10, 1000, 100]
})

pools, cost_pool, min_supply_pool, max_supply_pool, CO2_pool = gp.multidict({
    "P1": [5, 10, -1, 1],
    "P2": [1, 10, -1, 20],
    "P3": [2, 10, -1, 1],
    "P4": [1, 10, -1, 2]
})

targets, price_target, min_demand_target, max_demand_target, CO2_target, max_CI_target = gp.multidict({
    "T1": [15, 50, 80, 0, -1],
    "T2": [18, 60, 90, 0, -1],
    "T3": [20, 80, 100, 0, -1]
})

pools_CCS = []

""" Edges """

arcs, cost_arc, CO2_arc, max_arc = gp.multidict({
    # Source-to-pool
    ("S1", "P1", "M1"): [3, 1, -1],
    ("S1", "P2", "M1"): [3, 1, -1],
    ("S2", "P1", "M1"): [2, 1, -1],
    ("S2", "P2", "M1"): [2, 1, -1],
    ("S3", "P1", "M1"): [4, 3, -1],
    ("S3", "P2", "M1"): [4, 3, -1],
    # Pool-to-pool
    ("P1", "P3", "M2"): [1, 1, -1],
    ("P1", "P4", "M2"): [1, 1, -1],
    ("P2", "P3", "M2"): [1, 1, -1],
    ("P2", "P4", "M2"): [1, 1, -1],
    # Pool-to-targe
    ("P3", "T1", "M3"): [1, 1, -1],
    ("P3", "T2", "M3"): [1, 1, -1],
    ("P3", "T3", "M3"): [1, 1, -1],
    ("P4", "T1", "M3"): [1, 1, -1],
    ("P4", "T2", "M3"): [1, 1, -1],
    ("P4", "T3", "M3"): [1, 1, -1]
})