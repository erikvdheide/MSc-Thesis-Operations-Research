""" Model_MRMT

The Model for the Multi-Resource, Multi-Transport (MRMT) model.

This file contains the following components:
- Exact model formulation, without any carbon regulations (except pre-implementation CCS);
- Carbon flow and carbon intensity constraints;
- Two ways to calculate carbon-intensities post-optimization.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

""" Import packages & data """
import networkx as nx
from gurobipy import GRB
import time
import sys
from matplotlib import pyplot as plt

# Import the data and the parameters
from Params_MRMT import *
if useScaledData:
    from Data_MRMT_Scaled import *
else:
    from Data_MRMT import *

""" Create set of nodes, arcs and a graph """

# List of purely the nodes
nodes = []
nodes.extend(sources)
nodes.extend(prods)
nodes.extend(depots)
nodes.extend(targets)

# List of resources x nodes (if that resource is defined on that node)
nodes_res = []
nodes_res.extend(sources_res)
nodes_res.extend(prods_res)
nodes_res.extend(depots_res)
nodes_res.extend(targets_res)

# Set of all the arcs [(start, end)]
arcs2 = set()
for edge in arcs:
    arcs2.add((edge[0], edge[1]))

# Set of all the arcs [(start, end, resource)]
arcs3 = set()
for edge in arcs:
    arcs3.add((edge[0], edge[1], edge[3]))

# Create a graph
G = nx.DiGraph()
G.add_nodes_from(nodes)
G.add_edges_from(arcs2)

print("GRAPH: ")
print("Nodes: ", G.nodes)
print("Edges: ", G.edges)
print()

""" Define dictionaries, parameters and booleans from the data """

# Dictionary with the resources defined on each node
node_resourceDict = {}
for node_res in nodes_res:
    # If the node not yet in dictionary, add the first resource
    if node_res[0] not in node_resourceDict:
        node_resourceDict[node_res[0]] = [node_res[1]]
    # If the node already has mode in dictionary, add another resource
    else:
        node_resourceDict[node_res[0]].append(node_res[1])

# Dictionary with the mode defined on the arcs
arc_modeDict = {}
for arc in arcs:
    # If the arc not yet in dictionary, add the first mode
    if (arc[0], arc[1]) not in arc_modeDict:
        arc_modeDict[arc[0], arc[1]] = [arc[2]]
    # If the arc already has mode in dictionary, add another mode
    else:
        if arc[2] not in arc_modeDict[arc[0], arc[1]]:
            arc_modeDict[arc[0], arc[1]].append(arc[2])

# Calculate bigM as the maximum of all products for all customers
bigM = sum(max_demand_target[t] for t in targets_res)

# Check if we have CCS options
includedCCS = False
if len(prods_CCS) > 0:
    includedCCS = True

# In case you don't want to optimize over CCS
if disableCCS:
    includedCCS = False

# Create a dictionary of the investment options per pool where CCS is defined on
if includedCCS:
    CCSDict = {}
    for prod_level in prods_levels_CCS:
        if prod_level[0] not in CCSDict:
            CCSDict[prod_level[0]] = [prod_level[1]]
        else:
            CCSDict[prod_level[0]].append(prod_level[1])

""" Update CI constraints """
max_CI_target[("T1", "R5")] = maxCI_T1_R5
max_CI_target[("T1", "R6")] = maxCI_T1_R6
max_CI_target[("T2", "R5")] = maxCI_T2_R5
max_CI_target[("T2", "R6")] = maxCI_T2_R6
max_CI_target[("T3", "R5")] = maxCI_T3_R5
max_CI_target[("T4", "R6")] = maxCI_T4_R6
max_CI_total = maxCI_Total

############################## Set up model & decision variables ######################################################

""" Initiate model """
Model = gp.Model()
Model.params.nonConvex = 2
Model.params.timelimit = timeLimit
if not printStatus:
    Model.Params.LogToConsole = 0

""" Decision variables """

# Material flow variables [(arc=(start, end, mode, resource)]
x_arc = Model.addVars(arcs, lb=0.0, name="Resource flow variables on the arcs")

# Material flows through a node per resource [(node, resource)]
x = Model.addVars(nodes_res, lb=0.0, name="Material flow through a node")

if optimizeOption1:
    # Carbon intensity per resource at a node [(node, resource)]
    CI_var = Model.addVars(nodes_res, lb=0.0, name="Carbon intensity of resources at nodes")

    # Fraction of each resource that goes through a depot [(depot, resource)]
    frac = Model.addVars(depots_res, lb=0.0, name="Fraction of the resource at depot compared to other resources")

if optimizeOption2:
    # CO2 flow variables defined on the arcs [arc=(start, end, resource)]
    y_arc = Model.addVars(arcs3, lb=0.0, name="CO2 flow variables on the arcs")

    # Total carbon flows per resources at a node [(node, resource)]
    y = Model.addVars(nodes_res, lb=0.0, name="Carbon flows of resources at nodes")

    # Total flow coming into the production unit [(prod)]
    y_prod = Model.addVars(prods, lb=0.0, name="Carbon flow total in production unit")

""" Integer variables """

# Binary variables if source supplies resource [(source, resource)]
z_source_res = Model.addVars(sources_res, vtype=GRB.BINARY, name="Binary variables if source supplies a resource")

# Binary variables if depot is utilized [(depot)]
z_depot = Model.addVars(depots, vtype=GRB.BINARY, name="Binary variables depot is utilized")

""" CCS variables """

# One if investment level l is done at production unit p with CCS defined [(prod_CCS)]
if includedCCS:
    # Capture of CO2 at production unit p with CCS defined [(prod_CCS)]
    c_p = Model.addVars(prods_CCS, lb=0.0, name="Amount of captured CO2 at pool p")

    # One if investment level l is done at production unit p with CCS defined [(prod_CCS)]
    z_pl = Model.addVars(prods_levels_CCS, vtype=GRB.BINARY, name="Binary variable if CO2 is captured using that level of investment")

    # For balance constraints of CI, we need a variable for the capture per unit
    c_p_unit = Model.addVars(prods_CCS, lb=0.0, name="CO2 capture per unit")

############################## Model objective & constraints ##########################################################

""" MODEL as presented in the report """

# (4.2) Resource flows from sources
Model.addConstrs(x[sr] == x_arc.sum(sr[0], '*', '*', sr[1]) for sr in sources_res)

# (4.3) Resource flows into production unit
Model.addConstrs(x[pr] == x_arc.sum('*', pr[0], '*', pr[1]) for pr in prods_res if yields.at[pr] < 0)  # .at to get to the element

# (4.4) Resource flows out of production units
Model.addConstrs(x[pr] == x_arc.sum(pr[0], '*', '*', pr[1]) for pr in prods_res if yields.at[pr] > 0)

# (4.5) Resource flows through depots
Model.addConstrs(x[dr] == x_arc.sum('*', dr[0], '*', dr[1]) for dr in depots_res)

# (4.6) Resource flows into targets
Model.addConstrs(x[tr] == x_arc.sum('*', tr[0], '*', tr[1]) for tr in targets_res)

# (4.7) Revenue calculations
revenue = gp.quicksum(price_target[tr] * x[tr] for tr in targets_res)

# (4.8) Costs sources
cost_sources = gp.quicksum(cost_source[sr] * x[sr] for sr in sources_res)

# (4.9) Costs production units
cost_prods = gp.quicksum(cost_prod[pr[0]] * x[pr] for pr in prods_res if yields.at[pr] < 0)

# (4.10) Costs depots
cost_depots = gp.quicksum(fixed_cost_depot[d] * z_depot[d] for d in depots)

# (4.11) Costs on arcs from transportation
cost_transport = gp.quicksum(cost_arc[arc] * x_arc[arc] for arc in arcs)

# (4.58) Cost CCS
if includedCCS:
    cost_CCS = gp.quicksum(fixed_cost_install_CCS[(p, l)] * z_pl[(p, l)] for p in prods_CCS for l in CCSDict[p])

# (4.12) Total costs
if includedCCS:
    costs_total = cost_sources + cost_prods + cost_depots + cost_transport + cost_CCS
else:
    costs_total = cost_sources + cost_prods + cost_depots + cost_transport

# (4.1) Profit
profit = revenue - costs_total

# (4.64) Total source emissions
emissions_sources = gp.quicksum(carbon_source[sr] * x[sr] for sr in sources_res)

# (4.65) Total production emissions
emissions_prods = gp.quicksum(carbon_prod[pr[0]] * x[pr] for pr in prods_res if yields.at[pr] < 0)

# (4.66) Total emissions saved by CCS
emissions_CCS = -1 * gp.quicksum(c_p[p] for p in prods_CCS)

# (4.67) Total depot emissions
emissions_depots = gp.quicksum(fixed_carbon_depot[d] * z_depot[d] for d in depots)

# (4.68) Total transport emissions
emissions_transport = gp.quicksum(carbon_arc[arc] * x_arc[arc] for arc in arcs)

# (4.60) Total emissions (without upper bound)
if includedCCS:
    emissions_total = emissions_sources + emissions_prods + emissions_depots + emissions_transport + emissions_CCS
else:
    emissions_total = emissions_sources + emissions_prods + emissions_depots + emissions_transport

# Emissions of the T3, R5 subgraph:
if useT3_R5:
    em_sources = carbon_source[("S1", "R1")] * x[("S1", "R1")] + carbon_source[("S3", "R3")] * x[("S3", "R3")] + carbon_source[("S3", "R4")] * x[("S3", "R4")] \
                 + carbon_source[("S4", "R3")] * x[("S4", "R3")] + carbon_source[("S4", "R4")] * x[("S4", "R4")]
    em_prods = carbon_prod["P1"] * x[("P1", "R5")]
    em_depots = fixed_carbon_depot["D1"] * z_depot["D1"] + fixed_carbon_depot["D2"] * z_depot["D2"]
    em_transport = carbon_arc[("S1", "P1", "M1", "R1")] * x_arc[("S1", "P1", "M1", "R1")] \
                     + carbon_arc[("S3", "P1", "M1", "R3")] * x_arc[("S3", "P1", "M1", "R3")] \
                     + carbon_arc[("S3", "P1", "M1", "R4")] * x_arc[("S3", "P1", "M1", "R4")] \
                     + carbon_arc[("S4", "P1", "M1", "R3")] * x_arc[("S4", "P1", "M1", "R3")] \
                     + carbon_arc[("S4", "P1", "M1", "R4")] * x_arc[("S4", "P1", "M1", "R4")] \
                     + carbon_arc[("P1", "D1", "M2", "R5")] * x_arc[("P1", "D1", "M2", "R5")] \
                     + carbon_arc[("P1", "D1", "M3", "R5")] * x_arc[("P1", "D1", "M3", "R5")] \
                     + carbon_arc[("P1", "D2", "M2", "R5")] * x_arc[("P1", "D2", "M2", "R5")] \
                     + carbon_arc[("P1", "D2", "M3", "R5")] * x_arc[("P1", "D2", "M3", "R5")] \
                     + carbon_arc[("D1", "T3", "M4", "R5")] * x_arc[("D1", "T3", "M4", "R5")] \
                     + carbon_arc[("D2", "T3", "M4", "R5")] * x_arc[("D2", "T3", "M4", "R5")]
    em_CCS_T3_R5 = -1 * c_p["P1"]

    if includedCCS:
        em_total_T3_R5 = em_sources + em_prods + em_depots + em_transport + em_CCS_T3_R5
    else:
        em_total_T3_R5 = em_sources + em_prods + em_depots + em_transport

# Emissions of the T4, R6 subgraph:
if useT4_R6:
    em_sources = carbon_source[("S2", "R2")] * x[("S2", "R2")] + carbon_source[("S3", "R3")] * x[("S3", "R3")] + carbon_source[("S3", "R4")] * x[("S3", "R4")] \
                 + carbon_source[("S4", "R3")] * x[("S4", "R3")] + carbon_source[("S4", "R4")] * x[("S4", "R4")]
    em_prods = carbon_prod["P2"] * x[("P2", "R6")]
    em_depots = fixed_carbon_depot["D1"] * z_depot["D1"] + fixed_carbon_depot["D2"] * z_depot["D2"]
    em_transport = carbon_arc[("S2", "P2", "M1", "R2")] * x_arc[("S2", "P2", "M1", "R2")] \
                     + carbon_arc[("S3", "P2", "M1", "R3")] * x_arc[("S3", "P2", "M1", "R3")] \
                     + carbon_arc[("S3", "P2", "M1", "R4")] * x_arc[("S3", "P2", "M1", "R4")] \
                     + carbon_arc[("S4", "P2", "M1", "R3")] * x_arc[("S4", "P2", "M1", "R3")] \
                     + carbon_arc[("S4", "P2", "M1", "R4")] * x_arc[("S4", "P2", "M1", "R4")] \
                     + carbon_arc[("P2", "D1", "M2", "R6")] * x_arc[("P2", "D1", "M2", "R6")] \
                     + carbon_arc[("P2", "D1", "M3", "R6")] * x_arc[("P2", "D1", "M3", "R6")] \
                     + carbon_arc[("P2", "D2", "M2", "R6")] * x_arc[("P2", "D2", "M2", "R6")] \
                     + carbon_arc[("P2", "D2", "M3", "R6")] * x_arc[("P2", "D2", "M3", "R6")] \
                     + carbon_arc[("D1", "T4", "M4", "R6")] * x_arc[("D1", "T4", "M4", "R6")] \
                     + carbon_arc[("D2", "T4", "M4", "R6")] * x_arc[("D2", "T4", "M4", "R6")]
    em_CCS_T4_R6 = -1 * c_p["P2"]

    if includedCCS:
        em_total_T4_R6 = em_sources + em_prods + em_depots + em_transport + em_CCS_T4_R6
    else:
        em_total_T4_R6 = em_sources + em_prods + em_depots + em_transport

# OBJECTIVE: 100% profit, 100% emissions, or everything in between
if optimizeT3_R5:
    objective = -1 * em_total_T3_R5
elif optimizeT4_R6:
    objective = -1 * em_total_T4_R6
else:
    objective = profitImportance * profit + -1 * carbonImportance * emissions_total

""" Constraints """

# (4.13) Minimum supply + (4.14) Maximum supply
for sr in sources_res:
    # If no min supply restriction, then just not larger than maximum flow
    if min_supply_source[sr] == 0:
        if max_supply_source[sr] != -1:
            Model.addConstr(x[sr] <= max_supply_source[sr])
    # If min supply restriction, use of binary variable is necessary
    else:
        Model.addConstr(x[sr] >= min_supply_source[sr] * z_source_res[sr])
        if max_supply_source[sr] != -1:
            Model.addConstr(x[sr] <= max_supply_source[sr] * z_source_res[sr])
        else:
            Model.addConstr(x[sr] <= bigM * z_source_res[sr])

# (4.15) Maximum throughput of production unit
for p in prods:
    if max_through_prod[p] != -1:
        Model.addConstr(gp.quicksum(x[(p, r)] for r in node_resourceDict[p] if yields.at[(p, r)] < 0) <= max_through_prod[p])  # TODO does the if work here?

# (4.16), (4.17): Yield corrections at production units
for pr in prods_res:
    # (4.16) Yield correction of inflowing resources
    if yields.at[pr] < 0:
        Model.addConstr(x[pr] == -1 * yields.at[pr] * gp.quicksum(x[(pr[0], r)] for r in node_resourceDict[pr[0]] if yields.at[(pr[0], r)] > 0))
    # (4.17) Yield correction of outflowing resources
    if yields.at[pr] > 0:
        Model.addConstr(x[pr] == yields.at[pr] * gp.quicksum(x[(pr[0], r)] for r in node_resourceDict[pr[0]] if yields.at[(pr[0], r)] > 0))

# (4.18): Zero flow from resources not defined on the prod units
pass  # not necessary, because of how sets and arcs are defined

# (4.19) Inflow depot is also outflow depot
Model.addConstrs(x[dr] == x_arc.sum(dr[0], '*', '*', dr[1]) for dr in depots_res)

# (4.20) Minimum total flow through a depot
for d in depots:
    if min_through_depot[d] > 0:
        Model.addConstr(gp.quicksum(x[(d, r)] for r in node_resourceDict[d]) >= min_through_depot[d] * z_depot[d])

# (4.21) Maximum flow through depot
for dr in depots_res:
    if max_through_depot[dr] != -1:
        Model.addConstr(x[dr] <= max_through_depot[dr] * z_depot[dr[0]])
    else:
        if min_through_depot[dr[0]] != -1:
            Model.addConstr(x[dr] <= bigM * z_depot[dr[0]])

# (4.22) Demand satisfied at the targets
Model.addConstrs(x[tr] >= min_demand_target[tr] for tr in targets_res)
Model.addConstrs(x[tr] <= max_demand_target[tr] for tr in targets_res)

# (4.23) Maximum throughput on arcs
Model.addConstrs(x_arc[arc] <= max_through_arc[arc] for arc in arcs if max_through_arc[arc] != -1)

######################################## Carbon Constraints #######################################################

# Option 1: use balance constraints on CI values
if optimizeOption1:

    # (4.38): Carbon intensities at the sources
    for sr in sources_res:
        Model.addConstr(CI_var[sr] == carbon_source[sr])

    # (4.39), (4.40): Carbon intensities of incoming and outgoing resources production units
    for pr in prods_res:
        if yields.at[pr] < 0:
            Model.addConstr(CI_var[pr] * x[pr]
                            == gp.quicksum((CI_var[(arc[0], arc[3])] + carbon_arc[arc]) * x_arc[arc] for arc in arcs if arc[1] == pr[0] and arc[3] == pr[1]))
        if yields.at[pr] > 0:
            # (4.61), (4.62): CCS option
            if includedCCS and pr[0] in prods_CCS:
                Model.addConstr(c_p_unit[pr[0]] * gp.quicksum(x[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] < 0)
                                == c_p[pr[0]])
                Model.addConstr(CI_var[pr]
                                == gp.quicksum(CI_var[(pr[0], r)] * -1 * yields.at[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] < 0)
                                + carbon_prod[pr[0]] - c_p_unit[pr[0]])
            # (4.40)
            else:
                Model.addConstr(CI_var[pr]
                                == gp.quicksum(CI_var[(pr[0], r)] * -1 * yields.at[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] < 0) + carbon_prod[pr[0]])

    # (4.41), (4.42): Carbon intensities at depots
    for dr in depots_res:
        Model.addConstr(frac[dr] * gp.quicksum(x[(dr[0], r)] for r in resources if (dr[0], r) in depots_res) == x[dr])
        Model.addConstr(CI_var[dr] * x[dr]
                        == gp.quicksum((CI_var[(arc[0], arc[3])] + carbon_arc[arc]) * x_arc[arc] for arc in arcs if arc[1] == dr[0] and arc[3] == dr[1])
                        + fixed_carbon_depot[dr[0]] * frac[dr])

    # (4.43): Carbon intensities at targets
    for tr in targets_res:
        Model.addConstr(CI_var[tr] * x[tr]
                        == gp.quicksum((CI_var[(arc[0], arc[3])] + carbon_arc[arc]) * x_arc[arc] for arc in arcs if arc[1] == tr[0] and arc[3] == tr[1]))

    # (4.44): Carbon intensity restrictions
    for tr in targets_res:
        if max_CI_target[tr] != -1:
            Model.addConstr(CI_var[tr] <= max_CI_target[tr])

# Option two: use carbon flows directly and balance them over the products
if optimizeOption2:

    # (4.45): Carbon flows leaving the source
    for sr in sources_res:
        for j in G.successors(sr[0]):
            Model.addConstr(y_arc[(sr[0], j, sr[1])]
                            == gp.quicksum((carbon_source[sr] + carbon_arc[(sr[0], j, m, sr[1])]) * x_arc[(sr[0], j, m, sr[1])] for m in arc_modeDict[(sr[0], j)]))

    # (4.46): Total carbon flow entering production unit
    for p in prods:
        if includedCCS and p in prods_CCS:
            Model.addConstr(y_prod[p] == gp.quicksum(y_arc[(i, p, r)] for i in G.predecessors(p) for r in resources if (i, r) in nodes_res and yields.at[(p, r)] < 0)
                            + carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] < 0) - c_p[p])
        else:
            Model.addConstr(y_prod[p] == gp.quicksum(y_arc[(i, p, r)] for i in G.predecessors(p) for r in resources if (i, r) in nodes_res and yields.at[(p, r)] < 0)
                            + carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] < 0))

    # (4.47): Carbon flows leaving production units
    for pr in prods_res:
        for j in G.successors(pr[0]):
            if yields.at[pr] > 0 and (pr[0], j, pr[1]) in arcs3:
                Model.addConstr(y_arc[(pr[0], j, pr[1])] * gp.quicksum(x[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] > 0)
                                == y_prod[pr[0]] * gp.quicksum(x_arc[(pr[0], j, m, pr[1])] for m in arc_modeDict[(pr[0], j)])
                                + gp.quicksum(carbon_arc[(pr[0], j, m, pr[1])] * x_arc[(pr[0], j, m, pr[1])] for m in arc_modeDict[(pr[0], j)])
                                * gp.quicksum(x[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] > 0))

    # (4.48): Carbon flow entering the depot
    for dr in depots_res:
        Model.addConstr(y[dr] * gp.quicksum(x[(dr[0], r)] for r in resources if (dr[0], r) in depots_res)
                        == gp.quicksum(y_arc[(i, dr[0], dr[1])] for i in G.predecessors(dr[0]) if (i, dr[0], dr[1]) in arcs3) * gp.quicksum(
            x[(dr[0], r)] for r in resources if (dr[0], r) in depots_res)
                        + fixed_carbon_depot[dr[0]] * x[dr])

    # (4.49): Carbon flows leaving the depot + if-depot-not-used-correction
    for dr in depots_res:
        for j in G.successors(dr[0]):
            if (dr[0], j, dr[1]) in arcs3:
                Model.addConstr(y_arc[(dr[0], j, dr[1])] * x[dr]
                                == y[dr] * gp.quicksum(x_arc[(dr[0], j, m, dr[1])] for m in arc_modeDict[(dr[0], j)])
                                + gp.quicksum(carbon_arc[(dr[0], j, m, dr[1])] * x_arc[(dr[0], j, m, dr[1])] * x[dr] for m in arc_modeDict[(dr[0], j)]))
                # OLD constraint: manually putting flows from the depots to 0 if they are not used.
                #Model.addConstr(y_arc[(dr[0], j, dr[1])] == y_arc[(dr[0], j, dr[1])] * z_depot[dr[0]])

    # (4.50), (4.51): Flow, carbon intensity and max CI at the targets
    for tr in targets_res:
        # Flow into targets
        Model.addConstr(y[tr] == gp.quicksum(y_arc[(i, tr[0], tr[1])] for i in G.predecessors(tr[0]) if (i, tr[0], tr[1]) in arcs3))
        if max_CI_target[tr] != -1:
            Model.addConstr(y[tr] <= max_CI_target[tr] * x[tr])

################################## CCS & Linear emission constraint #######################################################

# CCS constraints
if includedCCS:

    for p in prods_CCS:
        # (4.52) The amount of CO2 captured is not more than the amount combusted
        Model.addConstr(c_p[p] <= carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] < 0))

        # (4.53) If something is captured, constraint on the corresponding maximum capture capacity
        Model.addConstr(c_p[p] <= gp.quicksum(max_capture_CCS[(p, l)] * z_pl[(p, l)] for l in CCSDict[p]))

        # (4.54) Add that for each pool, at most 1 investment can be chosen
        Model.addConstr(gp.quicksum(z_pl[(p, l)] for l in CCSDict[p]) <= 1)

# (4.69) Put a linear restriction on all emissions
if max_CI_total != -1:
    Model.addConstr(emissions_total <= max_CI_total)

# Put a restriction on the linear amount
if maxCI_Total_T3_R5 != -1:
    Model.addConstr(em_total_T3_R5 <= maxCI_Total_T3_R5)

# Put a restriction on the linear amount
if maxCI_Total_T4_R6 != -1:
    Model.addConstr(em_total_T4_R6 <= maxCI_Total_T4_R6)

# Let the model produce the maximum for R5:
if produceMaximumSubgraphT3_R5:
    Model.addConstr(x[("T1", "R5")] == max_demand_target[("T1", "R5")])
    Model.addConstr(x[("T2", "R5")] == max_demand_target[("T2", "R5")])
    Model.addConstr(x[("T3", "R5")] == max_demand_target[("T3", "R5")])

# Let the model produce the maximum for R6:
if produceMaximumSubgraphT4_R6:
    Model.addConstr(x[("T1", "R6")] == max_demand_target[("T1", "R6")])
    Model.addConstr(x[("T2", "R6")] == max_demand_target[("T2", "R6")])
    Model.addConstr(x[("T4", "R6")] == max_demand_target[("T4", "R6")])

############################## Execute & print functions ################################################################


def executeModel():
    """ Execute the model as specified in Params. """
    start_time = time.time()
    Model.setObjective(objective, GRB.MAXIMIZE)
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted(CI_Total):
    """ Execute the model with a new restriction on total emissions. """
    Model.addConstr(emissions_total <= CI_Total)
    Model.setObjective(objective, GRB.MAXIMIZE)
    Model.optimize()


def executeModelRestricted_R5_Option1(CI_R5):
    """ Execute the model with CI restricted on R5 using Option 1. """
    start_time = time.time()
    Model.addConstr(CI_var[("T3", "R5")] <= CI_R5)
    Model.setObjective(objective, GRB.MAXIMIZE)
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted_R6_Option1(CI_R6):
    """ Execute the model with CI restricted on R6 using Option 1. """
    start_time = time.time()
    Model.addConstr(CI_var[("T4", "R6")] <= CI_R6)
    Model.setObjective(objective, GRB.MAXIMIZE)
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted_R5_Option2(CI_R5):
    """ Execute the model with CI restricted on R5 using Option 2. """
    start_time = time.time()
    Model.addConstr(y[("T3", "R5")] <= CI_R5 * x[("T3", "R5")])
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted_R6_Option2(CI_R6):
    """ Execute the model with CI restricted on R6 using Option 2. """
    start_time = time.time()
    Model.addConstr(y[("T4", "R6")] <= CI_R6 * x[("T4", "R6")])
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted_R5_R6_Option1(CI_R5, CI_R6):
    """ Execute the model with CI restricted on R5 and R6 using Option 1. """
    start_time = time.time()
    Model.addConstr(CI_var[("T3", "R5")] <= CI_R5)
    Model.addConstr(CI_var[("T4", "R6")] <= CI_R6)
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def executeModelRestricted_R5_R6_Option2(CI_R5, CI_R6):
    """ Execute the model with CI restricted on R5 and R6 using Option 2. """
    start_time = time.time()
    Model.addConstr(y[("T3", "R5")] <= CI_R5 * x[("T3", "R5")])
    Model.addConstr(y[("T4", "R6")] <= CI_R6 * x[("T4", "R6")])
    Model.optimize()
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Tot. time   : {round(total_time, 3)} sec")


def getObjective():
    return objective.getValue()


def getTotalEmissions():
    return emissions_total.getValue()


def printObjective():
    """ Print the output of the model """
    print("Num. variables   : ", Model.numVars)  # (including those that are by definition 0)
    print("Num. constraints : ", Model.numConstrs)  # (including those that make variables 0 and redundant constraints for 0 variables)
    print()
    print(Model.printStats())
    print(f"OBJECTIVE:")
    print(f"Tot. Profit    : {round((revenue - costs_total).getValue(), 3)} [rev={round(revenue.getValue(), 3)}, cost={round(costs_total.getValue(), 3)}]")
    print(f"Tot. Emissions : {round(emissions_total.getValue(), 3)}")
    print(f"Optimality gap :", Model.MIPGap)
    print()
    #print(Model.printStats())  # prints statistics of the model
    #print()
    # for arc in arcs:
    #     print("arc:", arc, " - value: ", x_arc[arc].X)


def printVariables():
    """ Print the variables of the model """
    print(f"=== VARIABLES: ===")
    print()

    # Print outflow from each source, resource pair node
    print("SOURCES: ")
    for sr in sources_res:
        outflow_source = x[sr].X  # x_arc.sum(sr[0], '*', '*', sr[1]).getValue()
        print(f" Outflow from source {sr[0]} of resource {sr[1]}: {round(outflow_source, 3)}")
    print()

    # Print inflow and outflow for each production unit
    print("PRODUCTION UNITS: ")
    for p in prods:
        total_inflow_prod = 0
        total_outflow_prod = 0
        for r in node_resourceDict[p]:
            flow_prod = x[(p, r)].X
            if yields.at[(p, r)] < 0:
                total_inflow_prod += flow_prod
                print(f" inflow from prod. unit {p} of resource {r}: {round(flow_prod, 3)}")
            if yields.at[(p, r)] > 0:
                total_outflow_prod += flow_prod
                print(f" outflow from prod. unit {p} of resource {r}: {round(flow_prod, 3)}")
        print(f"  Total flow into prod. unit {p}: {round(total_inflow_prod, 3)}")
        print(f"  Total flow out prod. unit {p}: {round(total_inflow_prod, 3)}")
    print()

    # Print CCS captures at production units
    if includedCCS:
        print("CCS CAPTURES: ")
        for p in prods_CCS:
            print(f" Amount of CO2 captured at unit {p}: {round(c_p[p].X, 2)}")
            print(f"   Total combustion at unit {p}: "
                  f"{round(carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] > 0).getValue(), 2)}")
            for l in CCSDict[p]:
                if z_pl[(p, l)].X == 1:
                    print(f"   Capacity of investment made at unit {p}: {max_capture_CCS[(p, l)]}")
        print()

    # Print flow of the resources through depots
    print("DEPOTS: ")
    for dr in depots_res:
        flow_depot = x[dr].X
        print(f" Flow from source {dr[0]} of resource {dr[1]}: {round(flow_depot, 3)}")
    print()

    # Print total flow over different modes (specific for this dataset)
    print("MODES: ")
    total_modes = [0] * len(modes)
    for arc in arcs:
        for mode in modes:
            if arc[2] == mode:
                total_modes[modes.index(mode)] += x_arc[arc]
    for mode in modes:
        print(f" Total flow using mode {mode}: ", round(total_modes[modes.index(mode)].getValue(), 3))
    print()

    # Print demand satisfied at each target node
    print("TARGETS: ")
    for tr in targets_res:
        inflow_target = x[tr].X  # x_arc.sum('*', tr[0], '*', tr[1]).getValue()
        print(f" Demand satisfied at target {tr[0]} of resource {tr[1]}: {round(inflow_target, 3)}")
    print()

    # CI during optimization - Option 1
    if optimizeOption1:
        print("Carbon Intensities DURING-optimization way 1:")
        if printAllCI:
            for nr in nodes_res:
                print(f" CI of resource {nr[1]} at location {nr[0]}: {round(CI_var[nr].X, 2)}")
        else:
            for tr in targets_res:
                print(f" CI of resource {tr[1]} at location {tr[0]}: {round(CI_var[tr].X, 2)}")
        print()

    # CI during optimization - Option 2
    if optimizeOption2:
        # Final carbon intensities from optimization
        CI_flow_post = {}
        for tr in targets_res:
            CI_flow_post[tr] = y[tr].X / (x[tr].X + 10e-8)

        if printAllCI:
            print("Carbon flows DURING-optimization way 2:")
            for arc in arcs3:
                print(f" Flow on arc {arc}: {round(y_arc[arc].X, 2)}")
            for nr in nodes_res:
                if nr[0] not in prods:
                    print(f" Total carbon flow of node {nr[0]} of resource {nr[1]}: {round(y[nr].X, 2)}")
            for p in prods:
                print(f" Total flow in prod. unit {p}: {round(y_prod[p].X, 2)}")
            for tr in targets_res:
                print(f" CI of target {tr[0]} at resource {tr[1]}: {round(CI_flow_post[tr], 2)}")
        else:
            for tr in targets_res:
                print(f" CI of resource {tr[1]} at location {tr[0]}: {round(CI_flow_post[tr], 2)}")
        print()


def sourceChecker():
    """ checks if the lower-carbon source (S4) is activated. """
    if x[("S4", "R3")].X > 0.01 or x[("S4", "R4")].X > 0.01:
        return 1
    else:
        return 0


def CCSChecker():
    """ checks how many `levels' of CCS are activated (0, 1, 2, 3, 4, 5). """
    # small node: 100% carbon minimization will choose investment level L3 for P1, which is never useful when considering costs.
    CCS_counter_P1 = 0
    CCS_counter_P2 = 0
    if z_pl[("P1", "L1")].X > 0.99:
        CCS_counter_P1 = 1
    if z_pl[("P1", "L2")].X > 0.99:
        CCS_counter_P1 = 2
    if z_pl[("P1", "L3")].X > 0.99:
        CCS_counter_P1 = 3
    if z_pl[("P2", "L1")].X > 0.99:
        CCS_counter_P2 = 1
    if z_pl[("P2", "L2")].X > 0.99:
        CCS_counter_P2 = 2
    return CCS_counter_P1 + CCS_counter_P2


def modeChecker():
    """ checks if the lower-carbon mode (M3) is activated. """
    if x_arc[("P1", "D1", "M3", "R5")].X > 0.01 or x_arc[("P1", "D2", "M3", "R5")].X > 0.01 \
            or x_arc[("P2", "D1", "M3", "R6")].X > 0.01 or x_arc[("P2", "D2", "M3", "R6")].X > 0.01:
        return 1
    else:
        return 0


def depotChecker():
    """ checks if the lower-carbon depot (D2) is activated. """
    if x[("D2", "R5")].X > 0.01 or x[("D2", "R6")].X > 0.01:
        return 1
    else:
        return 0


def demandChecker():
    """ checks if not the maximum demand is satisfied at one of the target nodes. """
    demand_checker = 0
    if x[("T1", "R5")].X < (max_demand_target[("T1", "R5")]-0.01):
        demand_checker = 1
    if x[("T1", "R6")].X < (max_demand_target[("T1", "R6")]-0.01):
        demand_checker = 1
    if x[("T2", "R5")].X < (max_demand_target[("T2", "R5")]-0.01):
        demand_checker = 1
    if x[("T2", "R6")].X < (max_demand_target[("T2", "R6")]-0.01):
        demand_checker = 1
    if x[("T3", "R5")].X < (max_demand_target[("T3", "R5")]-0.01):
        demand_checker = 1
    if x[("T4", "R6")].X < (max_demand_target[("T4", "R6")]-0.01):
        demand_checker = 1
    return demand_checker


def postCalcOption1():
    """ Post calculate CI for Option 1: CI directly. """
    # Dictionary that stores CI values
    CI = {}
    for nr in nodes_res:
        CI[nr] = 0

    # (4.26): Carbon intensity at the sources
    for sr in sources_res:
        CI[sr] = carbon_source[sr]

    # (4.27): Carbon intensity of incoming resources of production units
    for pr in prods_res:
        if yields.at[pr] < 0:
            CI[pr] = gp.quicksum((CI[(i, pr[1])] + carbon_arc[(i, pr[0], m, pr[1])]) * (x_arc[(i, pr[0], m, pr[1])]) / (x[pr].X + 10e-8)
                                 for i in G.predecessors(pr[0]) if (i, pr[1]) in nodes_res for m in arc_modeDict[(i, pr[0])]).getValue()

    # (4.28): Carbon intensity of outgoing resources of production unit (for safety splitted)
    for pr in prods_res:
        if yields.at[pr] > 0:
            # (4.59) CCS option
            if includedCCS and pr[0] in prods_CCS:
                CI[pr] = gp.quicksum(CI[(pr[0], r)] * -1 * yields.at[(pr[0], r)] for r in resources if yields.at[(pr[0], r)]).getValue() + carbon_prod[pr[0]] \
                         - c_p[pr[0]].X / (gp.quicksum(x[(pr[0], r)] for r in resources if (pr[0], r) in prods_res and yields.at[(pr[0], r)] < 0).getValue() + 10e-8)
            else:
                CI[pr] = gp.quicksum(CI[(pr[0], r)] * -1 * yields.at[(pr[0], r)] for r in resources if yields.at[(pr[0], r)]).getValue() + carbon_prod[pr[0]]

    # (4.29): Carbon intensity of resources at the depots
    for dr in depots_res:
        if x[dr].X > 0:  # Must be at least some inflow
            CI[dr] = gp.quicksum((CI[(i, dr[1])] + carbon_arc[(i, dr[0], m, dr[1])]) * (x_arc[(i, dr[0], m, dr[1])]) / (x[dr].X + 10e-8)
                                 for i in G.predecessors(dr[0]) if (i, dr[1]) in nodes_res for m in arc_modeDict[(i, dr[0])]).getValue() \
                     + fixed_carbon_depot[dr[0]] / (gp.quicksum(x[(dr[0], r)] for r in resources if (dr[0], r) in depots_res).getValue() + 10e-8)

    # (4.30): Carbon intensity of the targets
    for tr in targets_res:
        CI[tr] = gp.quicksum((CI[(i, tr[1])] + carbon_arc[(i, tr[0], m, tr[1])]) * (x_arc[(i, tr[0], m, tr[1])]) / (x[tr].X + 10e-8)
                             for i in G.predecessors(tr[0]) if (i, tr[1]) in nodes_res for m in arc_modeDict[(i, tr[0])]).getValue()

    # Print option for debugging
    if printAllCI:
        for nr in nodes_res:
            print(f" CI of resource {nr[1]} at location {nr[0]}: {round(CI[nr], 2)}")

    return CI


def printCIOption1():
    """ Print CI for Option 1: CI directly """
    print("Carbon Intensities POST-optimization way 1:")
    total_emissions = 0
    for tr in targets_res:
        total_emissions += postCalcOption1()[tr] * x[tr]
        print(f" CI of resource {tr[1]} at location {tr[0]}: {round(postCalcOption1()[tr], 2)}")
    print(f" CHECK: CO2 emitted in total: {round(total_emissions.getValue(), 2)}")
    print()


def postCalcOption2():
    """ Post calculate CI for Option 1: CI directly """
    # Carbon flows on arcs (post-optimization)
    y_arc_post = {}
    for arc in arcs:
        y_arc_post[arc] = 0

    # Carbon flows
    y_nr_post = {}
    for nr in nodes_res:
        if nr[0] not in prods:
            y_nr_post[nr] = 0

    # Total carbon flow at the production units
    y_prod_post = {}
    for p in prods:
        y_p_post = 0

    # Final carbon intensities
    CI_post = {}
    for tr in targets_res:
        CI_post[tr] = 0

    # (4.31): Carbon flows leaving the source
    for arc in arcs:
        if arc[0] in sources:
            y_arc_post[arc] = (carbon_source[(arc[0], arc[3])] + carbon_arc[arc]) * x_arc[arc].X

    # (4.32): Total carbon flow through production unit
    for p in prods:
        # (4.60): CCS option
        if includedCCS and p in prods_CCS:
            y_prod_post[p] = (gp.quicksum(y_arc_post[arc] for arc in arcs if arc[1] == p and yields.at[(p, arc[3])] < 0)
                              + carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] < 0)).getValue() - c_p[p]
        else:
            y_prod_post[p] = (gp.quicksum(y_arc_post[arc] for arc in arcs if arc[1] == p and yields.at[(p, arc[3])] < 0)
                              + carbon_prod[p] * gp.quicksum(x[(p, r)] for r in resources if (p, r) in prods_res and yields.at[(p, r)] < 0)).getValue()

    # (4.33): Carbon flows leaving the production units
    for arc in arcs:
        if arc[0] in prods:
            y_arc_post[arc] = (y_prod_post[arc[0]] * (x_arc[arc]) / (gp.quicksum(x[(arc[0], r)] for r in resources if (arc[0], r) in prods_res and yields.at[(arc[0], r)] < 0).getValue() + 10e-8)
                               + carbon_arc[arc] * x_arc[arc]).getValue()

    # (4.34): Total carbon flows entering the depots
    for dr in depots_res:
        y_nr_post[dr] = (gp.quicksum(y_arc_post[arc] for arc in arcs if arc[1] == dr[0] and arc[3] == dr[1])
                         + fixed_carbon_depot[dr[0]] * (x[dr]) / (gp.quicksum(x[(dr[0], r)] for r in resources if (dr[0], r) in depots_res).getValue() + 10e-8)).getValue()

    # (4.35): Total carbon flows on arcs leaving the depots
    for arc in arcs:
        if arc[0] in depots:
            y_arc_post[arc] = (y_nr_post[(arc[0], arc[3])] * (x_arc[arc]) / (x[(arc[0], arc[3])].X + 10e-8)
                               + carbon_arc[arc] * x_arc[arc]).getValue()

    # (4.36): Total carbon flow coming into a depot
    for tr in targets_res:
        y_nr_post[tr] = (gp.quicksum(y_arc_post[arc] for arc in arcs if arc[1] == tr[0] and arc[3] == tr[1])).getValue()

    # (4.37): Print the CI values post-optimization
    for tr in targets_res:
        CI_post[tr] = y_nr_post[tr] / (x[tr].X + 10e-8)

    # Print option for debugging
    if printAllCI:
        print("Carbon flows POST-optimization way 2:")
        for arc in arcs:
            print(f" Flow on arc {arc}: {round(y_arc_post[arc], 2)}")
        for nr in nodes_res:
            if nr[0] not in prods:
                print(f" Total carbon flow of node {nr[0]} of resource {nr[1]}: {round(y_nr_post[nr], 2)}")
        for p in prods:
            print(f" Total flow in prod. unit {p}: {round(y_prod_post[p], 2)}")
        for tr in targets_res:
            print(f" CI of target {tr[0]} at resource {tr[1]}: {round(CI_post[tr], 2)}")

    return CI_post


def printCIOption2():
    """ Print CI for Option 2: Carbon flows """
    print("Carbon Intensities POST-optimization way 2:")
    total_emissions = 0
    for tr in targets_res:
        total_emissions += postCalcOption2()[tr] * x[tr]
        print(f" CI of resource {tr[1]} at location {tr[0]}: {round(postCalcOption2()[tr], 2)}")
    print(f" CHECK: CO2 emitted in total: {round(total_emissions.getValue(), 2)}")
    print()


def getCI_R5():
    return postCalcOption1()[("T3", "R5")]


def getCI_R6():
    return postCalcOption1()[("T4", "R6")]
