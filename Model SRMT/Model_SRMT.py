""" Model_SRMT

The single-resource, multi-transport model (with special cases of single-transport)

This model can be used for the distribution of one resource (production). It includes options for:
- multiple pools;
- pool-to-pool connections;
- source-to-target connections;
- multiple modes of transport;
- CCS investments.

The model can do the following:
- Optimize material flows without any carbon accounting;
- Have explicit constraints on carbon intensities (2 options);
- Post-calculate carbon intensities (2 options);
- Optimize over carbon intensities (though optimizing over emissions linearly is much faster!).

@author: Erik van der Heide

"""

############################## Decisions #############################################################################

usedData = 5  # Data 1, 2, 3, 4 of 5
timeLimit = 100  # Time limit in seconds
printStatus = True  # Print optimization status updates
minimizeOnlyCI = False  # Minimize all CIs for positive profit
profitImportance = 1.0  # Profit importance in weighted objective
correctCCSPost = True  # Correct CI post-optimization if not capacity captured from CCS, but was possible
balanceConstraintsWay1 = False  # Way 1 considers inflow and outflow, Way 2 only uses inflow + CI_d directly
atMostOneMode = False  # True to add the option that only 1 mode of transport can be used between A and B

optimizeOption1 = True  # Optimize using balance constraints
optimizeOption2 = False  # Optimize using carbon flows

postCalcOption1 = True  # Post-calc using CI vars
postCalcOption2 = False  # Post-calc using flow vars

# Restrictions on CIs:
maxCI_T1 = -1
maxCI_T2 = -1
maxCI_T3 = -1

############################## Imports & Data #########################################################################

""" Imports """
import networkx as nx
from gurobipy import GRB
import time
import sys

if usedData == 1:
    from Data_SRMT_1 import *
if usedData == 2:
    from Data_SRMT_2 import *
if usedData == 3:
    from Data_SRMT_3 import *
if usedData == 4:
    from Data_SRMT_4 import *
if usedData == 5:
    from Data_SRMT_5 import *

""" Create graph """

# List of all the nodes
nodes = []
nodes.extend(sources)
nodes.extend(pools)
nodes.extend(targets)

# Set of all the arcs without modes of transport
arcs2 = set()
for edge in arcs:
    arcs2.add((edge[0], edge[1]))

# Create a graph
G = nx.MultiDiGraph()
G.add_nodes_from(sources)
G.add_nodes_from(pools)
G.add_nodes_from(targets)
G.add_edges_from(arcs)

# For the arcs, we should define which modes of transport are defined on them
modesDict = {}
for edge in G.edges:
    # If the arc not yet in dictionary, add the first mode
    if (edge[0], edge[1]) not in modesDict:
        modesDict[edge[0], edge[1]] = [edge[2]]
    # If arc already has mode in dictionary, add another mode
    else:
        modesDict[edge[0], edge[1]].append(edge[2])

# Determine if there are multiple modes on SOME arcs
multiModesPerArc = False
for value in modesDict.values():
    if len(value) > 1:
        multiModesPerArc = True
        break

# Calculate bigM over the arcs
bigM = sum(max_demand_target[t] for t in targets)

# Check if we have CCS options
includedCCS = False
if len(pools_CCS) > 0:
    includedCCS = True

# Extract the unique pools on which CCS is defined
uniquePoolsCCS = set()
if includedCCS:
    for entry in pools_CCS:
        uniquePoolsCCS.add(entry[0])

# Create a dictionary of the investment options per pool where CCS is defined on
if includedCCS:
    CCSDict = {}
    for pl in pools_CCS:
        if pl[0] not in CCSDict:
            CCSDict[pl[0]] = [pl[1]]
        else:
            CCSDict[pl[0]].append(pl[1])

""" Update CI constraints """
max_CI_target['T1'] = maxCI_T1
max_CI_target['T2'] = maxCI_T2
max_CI_target['T3'] = maxCI_T3

############################## Set up model & decision variables ######################################################

""" Initiate model """
Model = gp.Model()
Model.params.nonConvex = 2
Model.params.timelimit = timeLimit
if not printStatus:
    Model.Params.LogToConsole = 0

""" Decision variables """

# Material flow variables [(arc=(start, finish, mode)]
x_arc = Model.addVars(arcs, lb=0.0, name="Material flow variables on the arcs")

# CO2 flow variables defined on the arcs [arc=(start, finish)]
y_arc = Model.addVars(arcs2, lb=0.0, name="CO2 flow variables on the arcs")

# Material flows through a node
x = Model.addVars(nodes, lb=0.0, name="Material flow through a node")

# Carbon intensity (total or per product) at a node
y = Model.addVars(nodes, lb=0.0, name="Carbon intensity of any node")

# Binary variables on source and pool nodes (only used if minimum supply is provided)
I = Model.addVars(nodes, vtype=GRB.BINARY, name="Binary variables if node is used")  # later also a "z" used

# Indicator whether arc mode m is used on arc i,j
if atMostOneMode:
    I_arc = Model.addVars(arcs, vtype=GRB.BINARY, name="Indicators for modes of transport")

# CCS variables
if includedCCS:
    c_p = Model.addVars(uniquePoolsCCS, lb=0.0, name="Amount of captured CO2 at pool p")
    I_pl = Model.addVars(pools_CCS, vtype=GRB.BINARY, name="Binary variable if CO2 is captured using that level of investment")

############################## Model objective & constraints ##########################################################

""" TEST: Fix if you fix CI variables (did not work!) """
# Model.addConstr(y['P1'] == 7.000000099067903 * 4.786)
# Model.addConstr(y['P2'] == 122.99999965945723 * 265.214)
# Model.addConstr(y['P3'] == 9.000000099034116 * 4.786)
# Model.addConstr(y['P4'] == 125.99999965927057 * 265.214)
# Model.addConstr(y['T1'] == 119.9999996799399 * 120.0)
# Model.addConstr(y['T2'] == 126.99999965884611 * 127.0)
# Model.addConstr(y['T3'] == 126.9999996588786 * 127.0)

""" Objective function """

# Revenue
revenue = gp.quicksum(price_target[t] * x_arc.sum('*', t, '*') for t in targets)

# Costs
total_cost_sources = gp.quicksum(cost_source[s] * x_arc.sum(s, '*', '*') for s in sources)
total_cost_pools = gp.quicksum(cost_pool[p] * x_arc.sum('*', p, '*') for p in pools)
total_cost_arc = gp.quicksum(cost_arc[arc] * x_arc[arc] for arc in arcs)
total_costs = total_cost_sources + total_cost_pools + total_cost_arc
if includedCCS:
    total_costs += gp.quicksum(cost_install_pool[pl] * I_pl[pl] for pl in pools_CCS)

# Profit
profit = revenue - total_costs

# Total CI of all final products combined
CI_targets = gp.quicksum(y_arc[(i, t)] for t in targets for i in G.predecessors(t))

# Objective
if not minimizeOnlyCI:
    # Profit maximization
    if profitImportance == 1.0:
        objective = profit
    # Balancing profit and sum of CI's
    else:
        objective = profitImportance * profit + (1 - profitImportance) * -1 * CI_targets
else:
    # Minimize CI with 0 profit
    objective = -1 * CI_targets
    Model.addConstr(profit >= 0)

""" Material constraints """

# The minimum and maximum amount of flow leaving a source
for s in sources:
    # If no min supply restriction, then just not larger than maximum flow
    if min_supply_source[s] == 0:
        if max_supply_source[s] != -1:
            Model.addConstr(x_arc.sum(s, '*', '*') <= max_supply_source[s])
    # If min supply restriction, use of binary variable is necessary
    else:
        Model.addConstr(x_arc.sum(s, '*', '*') >= min_supply_source[s] * I[s])
        if max_supply_source[s] != -1:
            Model.addConstr(x_arc.sum(s, '*', '*') <= max_supply_source[s] * I[s])
        else:
            Model.addConstr(x_arc.sum(s, '*', '*') <= bigM * I[s])

# The flow into a pool equals the flow out of a pool
Model.addConstrs(x_arc.sum('*', p, '*') == x_arc.sum(p, '*', '*') for p in pools)

# The minimum and maximum amount of flow of through a pool
for p in pools:
    # If no min supply restriction, then just not larger than maximum flow
    if min_supply_pool[p] == 0:
        if max_supply_pool[p] != -1:
            Model.addConstr(x_arc.sum(p, '*', '*') <= max_supply_pool[p])
    # If min supply restriction, use of binary variable is necessary
    else:
        Model.addConstr(x_arc.sum(p, '*', '*') >= min_supply_pool[p] * I[p])
        if max_supply_pool[p] != -1:
            Model.addConstr(x_arc.sum(p, '*', '*') <= max_supply_pool[p] * I[p])
        else:
            Model.addConstr(x_arc.sum(p, '*', '*') <= bigM * I[p])

# The satisfied demand at the target node is at least its minimum and at most its maximum
Model.addConstrs(x_arc.sum('*', t, '*') >= min_demand_target[t] for t in targets)
Model.addConstrs(x_arc.sum('*', t, '*') <= max_demand_target[t] for t in targets)

# Arc flows below maximum
Model.addConstrs(x_arc[arc] <= max_arc[arc] for arc in arcs if max_arc[arc] != -1)

# Only one mode of transport can be used between node i and j
if multiModesPerArc and atMostOneMode:
    Model.addConstrs(I_arc.sum(arc[0], arc[1], '*') <= 1 for arc in arcs2)  # TODO try if there are performance issues for == 1
    for arc in arcs:
        if max_arc[arc] == -1:
            Model.addConstr(x_arc[arc] <= I_arc[arc] * bigM)
        # Not sure if this will really improve the speed
        else:
            Model.addConstr(x_arc[arc] <= I_arc[arc] * max_arc[arc])

######################################## Carbon Constraints #######################################################

""" Option 1: standard balance constraints """

if optimizeOption1:

    for s in sources:
        # Carbon intensities at the source at given
        Model.addConstr(y[s] == CO2_source[s])

    for p in pools:
        # For CCS, don't distribute captured carbon over the outgoing arcs
        if p in uniquePoolsCCS:
            # Carbon intensity balance constraints
            if balanceConstraintsWay1:
                # Way 1: multiply the total sum OUTflow with y[p] on left-hand-side
                Model.addConstr(gp.quicksum((y[i] + CO2_arc[(i, p, m)] + CO2_pool[p]) * x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)])
                                == y[p] * gp.quicksum(x_arc[(p, j, m)] for j in G.successors(p) for m in modesDict[(p, j)]) + c_p[p])
            else:
                # Way 2: multiply the total sum INFLOW with y[p] on left-hand-side
                Model.addConstr(gp.quicksum((y[i] + CO2_arc[(i, p, m)] + CO2_pool[p]) * x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)])
                                == y[p] * gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]) + c_p[p])

        else:
            # Carbon intensity balance constraints
            if balanceConstraintsWay1:
                # Way 1: multiply the total sum OUTflow with y[p] on left-hand-side
                Model.addConstr(gp.quicksum((y[i] + CO2_arc[(i, p, m)] + CO2_pool[p]) * x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)])
                                == y[p] * gp.quicksum(x_arc[(p, j, m)] for j in G.successors(p) for m in modesDict[(p, j)]))
            else:
                # Way 2: multiply the total sum INFLOW with y[p] on left-hand-side
                Model.addConstr(gp.quicksum((y[i] + CO2_arc[(i, p, m)] + CO2_pool[p]) * x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)])
                                == y[p] * gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]))

    for t in targets:
        # Carbon intensity at the targets and constraints on it
        if balanceConstraintsWay1:
            # Way 1: calculate total inflow of carbon at the target and restrict on that proportionally
            Model.addConstr(y[t] == gp.quicksum((y[i] + CO2_arc[(i, t, m)]) * x_arc[(i, t, m)] for i in G.predecessors(t) for m in modesDict[(i, t)]))
            Model.addConstr(x[t] == x_arc.sum('*', t, '*'))
            if max_CI_target[t] != -1:
                Model.addConstr(y[t] <= max_CI_target[t] * x[t])
        else:
            # Way 2: use a balance constraint at the target nodes and restrict on the CI values directly
            Model.addConstr(gp.quicksum((y[i] + CO2_arc[(i, t, m)]) * x_arc[(i, t, m)] for i in G.predecessors(t) for m in modesDict[(i, t)])
                            == y[t] * gp.quicksum(x_arc[(i, t, m)] for i in G.predecessors(t) for m in modesDict[(i, t)]))
            if max_CI_target[t] != -1:
                Model.addConstr(y[t] <= max_CI_target[t])

""" Option 2: my way """

if optimizeOption2:

    for s in sources:
        for j in G.successors(s):
            # Determine carbon flows leaving a source
            Model.addConstr(y_arc[(s, j)] == gp.quicksum((CO2_source[s] + CO2_arc[(s, j, m)]) * x_arc[(s, j, m)] for m in modesDict[(s, j)]))

    for p in pools:
        # Determine material flows through pools
        Model.addConstr(x[p] == gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]))

        # Determine material flow OUT OF a pool
        if p in uniquePoolsCCS:  # Total carbon at the pool does not include captured carbon
            Model.addConstr(y[p] == gp.quicksum(y_arc[(i, p)] for i in G.predecessors(p)) + CO2_pool[p] * x[p] - c_p[p])
        else:
            Model.addConstr(y[p] == gp.quicksum(y_arc[(i, p)] for i in G.predecessors(p)) + CO2_pool[p] * x[p])

        # Determine carbon flows per pool
        for j in G.successors(p):
            Model.addConstr(y_arc[(p, j)] * x[p] == y[p] * gp.quicksum(x_arc[(p, j, m)] for m in modesDict[(p, j)])
                            + gp.quicksum(CO2_arc[(p, j, m)] * x_arc[(p, j, m)] * x[p] for m in modesDict[(p, j)]))

    for t in targets:
        # Determine material flow and CI of the target nodes
        Model.addConstr(x[t] == x_arc.sum('*', t, '*'))
        Model.addConstr(y[t] == y_arc.sum('*', t, '*'))

        # CI restrictions
        if max_CI_target[t] != -1:
            Model.addConstr(y[t] <= max_CI_target[t] * x[t])

""" Some final CCS constraints """
if includedCCS:
    for p in uniquePoolsCCS:
        # Constrained on actual combustion (save option avoiding x[p])
        Model.addConstr(c_p[p] <= CO2_pool[p] * gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]))
        # Constrained on maximum capture IF something captured, otherwise 0 is captured
        Model.addConstr(c_p[p] <= gp.quicksum(max_capture_pool[(p, l)] * I_pl[(p, l)] for l in CCSDict[p]))
        # Add that for each pool, at most 1 investment can be chosen
        Model.addConstr(gp.quicksum(I_pl[(p, l)] for l in CCSDict[p]) <= 1)

############################## Execute & print results ################################################################

""" Execute the model """

start_time = time.time()
Model.setObjective(objective, GRB.MAXIMIZE)
Model.optimize()
end_time = time.time()
total_time = end_time - start_time
print()

""" Print the output """

# Print objective value & time
if Model.Status == GRB.OPTIMAL:
    print('Model is feasible!')
elif Model.Status == GRB.INF_OR_UNBD:
    print('Model is infeasible or unbounded')
    sys.exit(0)
elif Model.Status == GRB.INFEASIBLE:
    print('Model is infeasible.')
    sys.exit(0)
elif Model.Status == GRB.UNBOUNDED:
    print('Model is unbounded.')
    sys.exit(0)
else:
    print('Optimization ended with status %d' % Model.Status)

print(f"Tot. Profit : {round((revenue-total_costs).getValue(),3)} [rev={round(revenue.getValue(),3)}, cost={round(total_costs.getValue(),3)}]")
print(f"Tot. CI     : {round(CI_targets.getValue(),3)}")
print(f"Tot. time   : {round(total_time, 3)} sec")
print()

# Print outflow from each source node
for s in sources:
    outflow_source = x_arc.sum(s, '*').getValue()
    print(f"Outflow from source {s}: {round(outflow_source,3)}")
print()

# Print material going through each pool node
for p in pools:
    inflow_pool = x_arc.sum('*', p).getValue()
    print(f"Flow through pool {p}: {round(inflow_pool, 3)}")
print()

# Print demand satisfied at each target node
CI_per_unit = []
for t in targets:
    inflow_target = x_arc.sum('*', t).getValue()
    if optimizeOption1:
        inflow_CO2_target = y[t].X
    if optimizeOption2:
        inflow_CO2_target = y_arc.sum('*', t).getValue()
    if (optimizeOption1 and balanceConstraintsWay1) or optimizeOption2:
        CI_per_unit.append(inflow_CO2_target / inflow_target)
    if optimizeOption1 and balanceConstraintsWay1==False:
        CI_per_unit.append(inflow_CO2_target)
    print(f"Demand satisfied at target {t}: {round(inflow_target,3)}")
print()

# Print CI values at the targets (post-optimization)
for CI_val in CI_per_unit:
    print(f"DURING-optimization Carbon intensity at target {t}: {round(CI_val,3)}")
print()

############################## Calculate CI post-optimization ###########################################################

for i in G.nodes:
    if i in sources:
        G.nodes[i]['CO2'] = CO2_source[i]
    if i in pools:
        G.nodes[i]['CO2'] = CO2_pool[i]
    if i in sources:
        G.nodes[i]['y_post'] = CO2_source[i]
    else:
        G.nodes[i]['y_post'] = 0

for arc in G.edges:
    G.edges[arc]['y_post'] = 0

""" Way 1: per-unit (before the pool) """

if postCalcOption1:

    for p in pools:

        # Calculate material flow into pools
        G.nodes[p]['x_post'] = gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]).getValue()

        # Calculate the CI at the pools
        for i in G.predecessors(p):
            for m in modesDict[(i, p)]:
                if i in sources:
                    G.nodes[p]['y_post'] += (G.nodes[i]['CO2'] + CO2_arc[(i, p, m)]) * (x_arc[(i, p, m)] / (G.nodes[p]['x_post']+10e-8))
                if i in pools:
                    G.nodes[p]['y_post'] += (G.nodes[i]['y_post'] + CO2_arc[(i, p, m)]) * (x_arc[(i, p, m)] / (G.nodes[p]['x_post']+10e-8))
        G.nodes[p]['y_post'] += CO2_pool[p]
        if p in uniquePoolsCCS:  # Remove captured CO2 per unit in post-calculation
            if correctCCSPost:
                for l in CCSDict[p]:
                    G.nodes[p]['y_post'] -= min(CO2_pool[p], (max_capture_pool[(p, l)]/G.nodes[p]['x_post'])) * I_pl[(p, l)]  # or first check if any I_pl is non-zero
            else:
                G.nodes[p]['y_post'] -= (c_p[p] / (G.nodes[p]['x_post'] + 10e-8))
        G.nodes[p]['y_post'] = G.nodes[p]['y_post'].getValue()

    for t in targets:

        # Calculate the material flow into targets
        G.nodes[t]['x_post'] = gp.quicksum(x_arc[(i, t, m)] for i in G.predecessors(t) for m in modesDict[(i, t)]).getValue()

        # Calculate the CI at the targets
        G.nodes[t]['y_post'] = gp.quicksum((G.nodes[i]['y_post'] + CO2_arc[(i, t, m)]) * (x_arc[(i, t, m)] / (G.nodes[t]['x_post']+10e-8))
                                            for i in G.predecessors(t) for m in modesDict[(i, t)]) + CO2_target[t]

    # Print the CI values of each target node (post-optimization)
    for t in targets:
        if balanceConstraintsWay1==False and optimizeOption1:
            CI_per_unit = y[t].X
        else:
            CI_per_unit = G.nodes[t]['y_post'].getValue()
        print(f"POST-optimization Carbon intensity at target {t}: {round(CI_per_unit, 3)}")
    print()


""" Way 2: total carbon flows (after the pool) """

if postCalcOption2:

    # Calculate flows leaving the source
    # NOTE: for post-calculations, it is fine to have multiple carbon flows between 2 nodes
    for s in sources:
        for j in G.successors(s):
            for m in modesDict[(s, j)]:
                G.edges[(s, j, m)]['y_post'] = (CO2_source[s] + CO2_arc[(s, j, m)]) * x_arc[(s, j, m)]

    # Calculate material flows into the pools
    for p in pools:
        G.nodes[p]['x_post'] = gp.quicksum(x_arc[(i, p, m)] for i in G.predecessors(p) for m in modesDict[(i, p)]).getValue()

    # Calculate material flows into the target
    for t in targets:
        G.nodes[t]['x_post'] = gp.quicksum(x_arc[(i, t, m)] for i in G.predecessors(t) for m in modesDict[(i, t)]).getValue()

    # Calculate flows in and out pools
    for p in pools:
        # Calculate the TOTAL flow that goes out of a node
        G.nodes[p]['y_post'] = (gp.quicksum(G.edges[(i, p, m)]['y_post'] for i in G.predecessors(p) for m in modesDict[(i, p)]) + (CO2_pool[p] * G.nodes[p]['x_post'])).getValue()
        if p in uniquePoolsCCS:  # Remove captured flow from total flow
            if correctCCSPost:
                for l in CCSDict[p]:
                    G.nodes[p]['y_post'] -= min(CO2_pool[p] * G.nodes[p]['x_post'], max_capture_pool[(p, l)]) * I_pl[(p, l)]
            else:
                G.nodes[p]['y_post'] -= c_p[p]

        # Distribute the CO2 flows over the outgoing arcs
        for j in G.successors(p):
            for m in modesDict[(p, j)]:
                G.edges[(p, j, m)]['y_post'] = G.nodes[p]['y_post'] * (x_arc[(p, j, m)] / (G.nodes[p]['x_post']+10e-8)) + CO2_arc[(p, j, m)] * x_arc[(p, j, m)]

    # Calculate CI of the targets
    for t in targets:
        for i in G.predecessors(t):
            for m in modesDict[(i, t)]:
                G.nodes[t]['y_post'] += G.edges[(i, t, m)]['y_post']

    # Print the CI values of each target node (post-optimization)
    for t in targets:
        CI_per_unit = G.nodes[t]['y_post'].getValue() / G.nodes[t]['x_post']
        print(f"POST-optimization Carbon intensity at target {t}: {round(CI_per_unit,3)}")
    print()

# Print CIs
# for n in nodes:
#     if x[n].X > 0.0001:
#         print("node: ", n, " with CI: ", (y[n].X/x[n].X))

