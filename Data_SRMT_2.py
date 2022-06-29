""" Data_SRMT_2

Data for set-up 2: 2 pools
NOTE: -1 = infinity.

@author: Erik van der Heide

"""

import gurobipy as gp

""" Nodes """

sources, cost_source, min_supply_source, max_supply_source, CO2_source = gp.multidict({
    "S1": [10, 20, 100, 5],
    "S2": [15, 20, 120, 4],
    "S3": [1, 0, 1000, 100]
})

pools, cost_pool, min_supply_pool, max_supply_pool, CO2_pool = gp.multidict({
    "P1": [5, 20, 150, 1],
    "P2": [1, 0, -1, 20]
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
    # Pool-to-target
    ("P1", "T1", "M2"): [1, 1, -1],
    ("P1", "T2", "M2"): [1, 1, -1],
    ("P1", "T3", "M2"): [1, 1, -1],
    ("P2", "T1", "M2"): [1, 1, -1],
    ("P2", "T2", "M2"): [1, 1, -1],
    ("P2", "T3", "M2"): [1, 1, -1]
})