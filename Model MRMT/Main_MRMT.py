""" Main_MRMT

The Main of the Multi-Resource, Multi-Transport (MRMT) model.

This file performs runs on the MRMT model, based on the parameters that are chosen in Params_MRMT.
The model interacts with Model_MRMT, which on its turn interacts with Data_MRMT or Data_MRMT_Scaled.

The Main file can do the following, which have to be un-commented to be activated:
- Basic execution of the model; printing objective; printing variables; printing CI post-optimization with option 1 or 2;
- Perform progressively stricter linear restrictions and plot the objective, CI, and some model indicators for these restrictions;
- Perform progressively stricter CI restrictions for R5 and/or R6 for Option 1 and Option 2;
- Some hardcoded data to ease the runs.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

from Model_MRMT import *
from matplotlib import pyplot as plt

############################## Execute & print results ################################################################

""" Basic execution of the model """

executeModel()
printObjective()
printVariables()

printCIOption1()
printCIOption2()


"""
===========================================================================
Section: 5.1.3. - Absolute Emission Restrictions 
Perform a stricter linear restriction in steps, until no longer achievable
===========================================================================
"""

""" 
# Run the linear models and save all the results
objectives = []
x_axis = []
total_emissions = []
CI_R5 = []
CI_R6 = []
source_list = []
CCS_list = []
mode_list = []
depot_list = []
demand_list = []
rangeUsed = range(20793, 7500, -20)  # 7500
for i in rangeUsed:
    executeModelRestricted(i)
    x_axis.append(i)
    objectives.append(getObjective())
    total_emissions.append(getTotalEmissions())
    CI_R5.append(getCI_R5())
    CI_R6.append(getCI_R6())
    source_list.append(sourceChecker())
    CCS_list.append(CCSChecker())
    mode_list.append(modeChecker())
    depot_list.append(depotChecker())
    demand_list.append(demandChecker())

print("Objectives     :", objectives)
print("Tot. emissions :", total_emissions)
print("CI R5          :", CI_R5)
print("CI R6          :", CI_R6)
print("Source list :", source_list)
print("CCS list    :", CCS_list)
print("Mode list   :", mode_list)
print("Depot list  :", depot_list)
print("Demand list :", demand_list)
"""

"""
# Print profit objective
plt.plot(rangeUsed, objectives, 'blue', label='profit')
plt.xlabel('Total chain emissions maximum')
plt.ylabel('Objective')
plt.grid(axis='y', alpha=0.75)
plt.gca().invert_xaxis()
plt.legend()
plt.show()
"""

"""
# Print carbon intensities
CCS_check = 1
firstCCS_check = False
for j in CCS_list:
    if j == CCS_check:
        if firstCCS_check == False:
            plt.axvline(x=x_axis[CCS_list.index(j)], linewidth=2, color='lightgray', label='CCS investment')
            firstCCS_check = True
        else:
            plt.axvline(x=x_axis[CCS_list.index(j)], linewidth=2, color='lightgray')
        CCS_check = CCS_check + 1
plt.plot(rangeUsed, CI_R5, linestyle='--', color='red', label='Carbon Intensity R5')
plt.plot(rangeUsed, CI_R6, linestyle='-.', color='green', label='Carbon Intensity R6')
plt.xlabel('Total chain emissions maximum')
plt.ylabel('Carbon intensity')
plt.grid(axis='y', alpha=0.75)
plt.gca().invert_xaxis()
plt.legend(loc='upper right')
plt.show()
"""

"""
# Print when sources, modes, depots, demand and CCS are used
for source_el in source_list:
    source_list[source_list.index(source_el)] = source_el * 4
for mode_el in mode_list:
    mode_list[mode_list.index(mode_el)] = mode_el * 3
for depot_el in depot_list:
    depot_list[depot_list.index(depot_el)] = depot_el * 2

CCS_check = 1
firstCCS_check = False
for j in CCS_list:
    if j == CCS_check:
        if firstCCS_check == False:
            plt.axvline(x=x_axis[CCS_list.index(j)], linewidth=2, color='lightgray', label='CCS investment')
            firstCCS_check = True
        else:
            plt.axvline(x=x_axis[CCS_list.index(j)], linewidth=2, color='lightgray')
        CCS_check = CCS_check + 1

plt.plot(rangeUsed, source_list, 's', markersize=2)
plt.text(x=7350, y=3.95, s="sources", fontsize=10)
plt.plot(rangeUsed, mode_list, 's', markersize=2)
plt.text(x=7350, y=2.95, s="modes", fontsize=10)
plt.plot(rangeUsed, depot_list, 's', markersize=2)
plt.text(x=7350, y=1.95, s="depots", fontsize=10)
plt.plot(rangeUsed, demand_list, 's', markersize=2)
plt.text(x=7350, y=0.95, s="demand", fontsize=10)
plt.ylim(0.5, 4.5)
plt.xlabel('Total chain emissions maximum')
plt.gca().invert_xaxis()
plt.xlim(x_axis[0]+500, 5400)
plt.show()
"""

""" 
===================================================================
 Section 5.1.3. - Comparing CI formulations
 Perform CI restrictions with step-by-step stricter CI restrictions
===================================================================
"""

"""
# Restriction on CI of R5 in steps of one
rangeCI_R5 = range(44, 15, -1)
for i in rangeCI_R5:
    print("====RESTRICTION R5: ", i, "====")
    executeModelRestricted_R5_Option2(i)
    printObjective()
    printCIOption2()
"""

"""
# Restriction on CI of R6 in steps of one
rangeCI_R6 = range(36, 12, -1)
for i in rangeCI_R6:
    print("====RESTRICTION R6: ", i, "====")
    executeModelRestricted_R6_Option2(i)
    printObjective()
    printCIOption2()
"""

"""
# Perform some CI restrictions on both R5 and R6 simultaneously (5 options in total)
# MANUALLY change if you perform Option 1 or Option 2
print("===== RESTRICTION CI_R5=22, CI_R6=19: =====")
executeModelRestricted_R5_R6_Option2(22, 19)
printObjective()
printCIOption2()

print("===== RESTRICTION CI_R5=21, CI_R6=18: =====")
executeModelRestricted_R5_R6_Option2(21, 18)
printObjective()
printCIOption2()

print("===== RESTRICTION CI_R5=20, CI_R6=17: =====")
executeModelRestricted_R5_R6_Option2(20, 17)
printObjective()
printCIOption2()

print("===== RESTRICTION CI_R5=19, CI_R6=16: =====")
executeModelRestricted_R5_R6_Option2(19, 16)
printObjective()
printCIOption2()

print("===== RESTRICTION CI_R5=18, CI_R6=15: =====")
executeModelRestricted_R5_R6_Option2(18, 15)
printObjective()
printCIOption2()
"""


""" Hardcoded information on Min-Max Data for Data_MRMT """
"""
# Highest profit solution: 10763.0
# Total emissions: 20792.5
# CI R5: 44.41
# CI R6: 36.41

# Lowest total emissions solution: 7500.5
# Profit: 4328.75 (3328.75 in minimization, but correct for P1, L3)
# CI R5: 17.91
# CI R6: 15.09

# Minimum emissions if always maximum demand: 8584.5 (profit 5551.25)
# CI R5: 18.4
# CI R6: 14.97
"""

