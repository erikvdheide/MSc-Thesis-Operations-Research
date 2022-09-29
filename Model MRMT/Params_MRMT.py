""" Params_MRMT

The Parameters of the Multi-Resource, Multi-Transport (MRMT) model.

All parameters are documented - see below.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

# General
timeLimit = 100  # maximum time limit of the exact model in seconds
printStatus = True  # print progress while running - mainly for best-found objective and optimality gap
profitImportance = 1.0  # how much you want to optimize over profit (between 0.0 and 1.0)
carbonImportance = 0.0  # how much you want to optimize of carbon emissions in the chain (between 0.0 and 1.0)
printAllCI = False  # put on True to also print CI of nodes other than targets
disableCCS = False  # put on True to not included CCS in the model (even if there is CCS data in the dataset)
useScaledData = False  # put on True to use Data_MRMT_Scaled instead of Data_MRMT

# CI optimization options
optimizeOption1 = True  # put on True to use the formulation based on CI balance constraints
optimizeOption2 = False  # put on True to use the formulation based on arc flow variable distribution

# Subgraphs
useT3_R5 = False  # put on True to define the subgraph emissions of T3, R5
useT4_R6 = False  # put on True to define the subgraph emissions of T4, R6
optimizeT3_R5 = False  # put on True to find the minimum emissions of T3, R5
optimizeT4_R6 = False  # put on True to find the minimum emissions of T4, R6
produceMaximumSubgraphT3_R5 = False  # put on True to impose producing the maximum resource of  R5
produceMaximumSubgraphT4_R6 = False  # put on True to impose producing the maximum resource of R6

# Restrictions on CI for (target, resource) pairs [-1 = no restriction]
maxCI_T1_R5 = -1
maxCI_T1_R6 = -1
maxCI_T2_R5 = -1
maxCI_T2_R6 = -1
maxCI_T3_R5 = -1
maxCI_T4_R6 = -1
# Restrictions on LINEAR emissions^max for the total chain and two subgraphs [-1 = no restriction]
maxCI_Total = -1
maxCI_Total_T3_R5 = -1
maxCI_Total_T4_R6 = -1
