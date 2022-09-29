""" Params_case

The Parameters and choices of the case study model.

All parameters are documented - see below.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

# Some settings
import pandas as pd
import numpy as np
import math
pd.options.display.show_dimensions = False

""" Choices """
# Main model type
solveTimePeriod = 4  # For the STATIC model, solve time period 1, 2, 3 or 4
timeLimit = 10*60  # Time limit in minutes
onlyOneImportProduct = False  # True to only be able to import one hydrogen product in a grid
costImportance = 1.0  # between 0.0 and 1.0
CO2Importance = 0.0  # between 0.0 and 1.0
useStartSolution = False  # True to use a (CO2-minimizing) starting solution

# Prints
printStatus = True  # True to print the status while running
printData = True  # True to print the data that is used
printVariables = True  # True to print the variables of the model
printBinaryVariables = False  # True to print the binary variables X, Y, Z
printTests = False  # True to print some tests

# Settings in the back
numGrids = 25  # this still has to be hardcoded
numProds = 2  # this still has to be hardcoded
numPeriods = 4  # this still has to be hardcoded
avgCapSmall = math.floor(20+(99-20)/2)  # Where TCC is based on for small units
avgCapMedium = math.floor(100+(499-100)/2)  # Where TCC is based on for medium units
avgCapLarge = math.floor(500+(1000-500)/2)  # Where TCC is based on for large units
distanceWithinGrid = 5  # 0 to assume no transport cost in grid, otherwise take avg. transport distance at location
CCSCostInput = 25  # $/ton CO2 produced; variable between 5 and 55

# Extra choices which are usually fixed
optimizeCI = True  # True to optimize carbon intensity variables
includeCCS = True  # True to include CCS as option
multipleSizes = True  # True to include plants with sizes small, medium and large
underEstimateNTU = True  # True to calculate the number of units as they can take each others trips
overEstimateNTU = False  # True to calculate more transportation units than might be needed
minCapZero = True  # True to put the minimum capacity of plants to 0 in the dynamic case
fixPlantsStatic = False  # True to fix some plants in the static model (or specify not what is fixed)
newEmissionData = True  # True to reduce the standard biomass emissions with 80%.
groupByTechnology = True  # True to group plants in the static model by their technology

""" Carbon restrictions """

# Total CO2 emissions
TotalCO2Max = -1
TotalCO2Max_t = [-1, -1, -1, -1]

# Individual CI restrictions (CI^{max}_{ig} and CI^{max}_{g})
maxCIStatic = -1 * np.ones(shape=((numProds+1), numGrids))
maxCIDynamic = -1 * np.ones(shape=((numProds+1), numGrids, numPeriods))


def setRestrictionStatic(gridNr, restrictionValue, productType):
    """ Function to set a restriction a carbon intensity
    :param gridNr: The grid for which you want to put a restriction;
    :param restrictionValue: The value of the CI restriction;
    :param productType: Put restriction on "CH2", "LH2", or "ALL".  """

    # Exceptions
    if gridNr > numGrids or gridNr <= 0:
        raise ValueError("Grid number is too high (or zero/negative).")
    elif restrictionValue < 0:
        raise ValueError("CI restriction must be non-negative.")
    elif productType.upper() not in ["CH2", "LH2", "ALL"]:
        raise ValueError("ProductType arguments must be 'CH2', 'LH2' or 'ALL'.")
    # Set the values
    elif productType.upper() == "CH2":
        maxCIStatic[(0, gridNr-1)] = restrictionValue
    elif productType.upper() == "LH2":
        maxCIStatic[(1, gridNr-1)] = restrictionValue
    elif productType.upper() == "ALL":
        maxCIStatic[(2, gridNr-1)] = restrictionValue


def setRestrictionDynamic(gridNr, restrictionValue, productType, timePeriod):
    """ Function to set a restriction a carbon intensity
    :param gridNr: The grid for which you want to put a restriction;
    :param restrictionValue: The value of the CI restriction;
    :param productType: Put restriction on "CH2", "LH2", or "ALL";
    :param timePeriod: The time period for which the restriction should hold. """

    # Exceptions
    if gridNr > numGrids or gridNr <= 0:
        raise ValueError("Grid number is too high (or zero/negative).")
    elif restrictionValue < 0:
        raise ValueError("CI restriction must be non-negative.")
    elif productType.upper() not in ["CH2", "LH2", "ALL"]:
        raise ValueError("ProductType arguments must be 'CH2', 'LH2' or 'ALL'.")
    elif timePeriod > numPeriods or timePeriod <= 0:
        raise ValueError("Time period is too high (or zero/negative).")
    # Set the values
    elif productType.upper() == "CH2":
        maxCIDynamic[(0, gridNr-1, timePeriod-1)] = restrictionValue
    elif productType.upper() == "LH2":
        maxCIDynamic[(1, gridNr-1, timePeriod-1)] = restrictionValue
    elif productType.upper() == "ALL":
        maxCIDynamic[(2, gridNr-1, timePeriod-1)] = restrictionValue


""" Set the actual CI restrictions """

""" STATIC restrictions """

## NON-PRODUCING CITIES
# setRestrictionStatic(2, 10, "all")
# setRestrictionStatic(6, 10, "all")
# setRestrictionStatic(11, 10, "all")
# setRestrictionStatic(20, 10, "all")
# setRestrictionStatic(24, 10, "all")

# setRestrictionStatic(2, 1, "all")
# setRestrictionStatic(6, 5, "all")
# setRestrictionStatic(11, 5, "all")
# setRestrictionStatic(20, 5, "all")
# setRestrictionStatic(24, 5, "all")

# setRestrictionStatic(2, 5, "CH2")
# setRestrictionStatic(2, 5, "LH2")

## PRODUCING CITIES
# setRestrictionStatic(1, 10, "all")
# setRestrictionStatic(5, 10, "all")
# setRestrictionStatic(10, 10, "all")
# setRestrictionStatic(19, 10, "all")
# setRestrictionStatic(23, 10, "all")
# setRestrictionStatic(1, 0, "all")
# setRestrictionStatic(5, 5, "all")
# setRestrictionStatic(10, 5, "all")
# setRestrictionStatic(19, 5, "all")
# setRestrictionStatic(23, 5, "all")

## COMBINATION OF CITIES
# setRestrictionStatic(2, 10, "all")
# setRestrictionStatic(24, 10, "all")

# setRestrictionStatic(1, 10, "all")
# setRestrictionStatic(24, 10, "all")

# setRestrictionStatic(2, 10, "all")
# setRestrictionStatic(23, 10, "all")

# setRestrictionStatic(1, 10, "all")
# setRestrictionStatic(23, 10, "all")

# setRestrictionStatic(1, 10, "all")
# setRestrictionStatic(2, 10, "all")
# setRestrictionStatic(3, 10, "all")
# setRestrictionStatic(4, 10, "all")
# setRestrictionStatic(5, 10, "all")
# setRestrictionStatic(6, 10, "all")
# setRestrictionStatic(7, 10, "all")
# setRestrictionStatic(8, 10, "all")
# setRestrictionStatic(9, 10, "all")
# setRestrictionStatic(10, 10, "all")
# setRestrictionStatic(11, 10, "all")
# setRestrictionStatic(12, 10, "all")
# setRestrictionStatic(13, 10, "all")
# setRestrictionStatic(14, 10, "all")
# setRestrictionStatic(15, 10, "all")
# setRestrictionStatic(16, 10, "all")
# setRestrictionStatic(17, 10, "all")
# setRestrictionStatic(18, 10, "all")
# setRestrictionStatic(19, 10, "all")
# setRestrictionStatic(20, 10, "all")
# setRestrictionStatic(21, 10, "all")
# setRestrictionStatic(22, 10, "all")
# setRestrictionStatic(23, 10, "all")
# setRestrictionStatic(24, 10, "all")
# setRestrictionStatic(25, 10, "all")

# setRestrictionStatic(2, 5, "all")
# setRestrictionStatic(24, 5, "all")

# setRestrictionStatic(1, 5, "all")
# setRestrictionStatic(24, 5, "all")

# setRestrictionStatic(2, 5, "all")
# setRestrictionStatic(23, 5, "all")

# setRestrictionStatic(1, 5, "all")
# setRestrictionStatic(23, 5, "all")

# setRestrictionStatic(1, 5, "all")
# setRestrictionStatic(2, 5, "all")
# setRestrictionStatic(3, 5, "all")
# setRestrictionStatic(4, 5, "all")
# setRestrictionStatic(5, 5, "all")
# setRestrictionStatic(6, 5, "all")
# setRestrictionStatic(7, 5, "all")
# setRestrictionStatic(8, 5, "all")
# setRestrictionStatic(9, 5, "all")
# setRestrictionStatic(10, 5, "all")
# setRestrictionStatic(11, 5, "all")
# setRestrictionStatic(12, 5, "all")
# setRestrictionStatic(13, 5, "all")
# setRestrictionStatic(14, 5, "all")
# setRestrictionStatic(15, 5, "all")
# setRestrictionStatic(16, 5, "all")
# setRestrictionStatic(17, 5, "all")
# setRestrictionStatic(18, 5, "all")
# setRestrictionStatic(19, 5, "all")
# setRestrictionStatic(20, 5, "all")
# setRestrictionStatic(21, 5, "all")
# setRestrictionStatic(22, 5, "all")
# setRestrictionStatic(23, 5, "all")
# setRestrictionStatic(24, 5, "all")
# setRestrictionStatic(25, 5, "all")

""" DYNAMIC restrictions """

# setRestrictionDynamic(2, 10, "all", 1)
# setRestrictionDynamic(2, 8, "all", 2)
# setRestrictionDynamic(2, 6, "all", 3)
# setRestrictionDynamic(2, 4, "all", 4)

# setRestrictionDynamic(1, 10, "all", 1)
# setRestrictionDynamic(1, 8, "all", 2)
# setRestrictionDynamic(1, 6, "all", 3)
# setRestrictionDynamic(1, 4, "all", 4)

# setRestrictionDynamic(2, 8, "all", 1)
# setRestrictionDynamic(2, 6, "all", 2)
# setRestrictionDynamic(2, 4, "all", 3)
# setRestrictionDynamic(2, 0, "all", 4)

# setRestrictionDynamic(1, 8, "all", 1)
# setRestrictionDynamic(1, 6, "all", 2)
# setRestrictionDynamic(1, 4, "all", 3)
# setRestrictionDynamic(1, 2, "all", 4)

# # Period 1
# setRestrictionDynamic(1, 15, "all", 1)
# setRestrictionDynamic(2, 15, "all", 1)
# setRestrictionDynamic(3, 15, "all", 1)
# setRestrictionDynamic(4, 15, "all", 1)
# setRestrictionDynamic(5, 15, "all", 1)
# setRestrictionDynamic(6, 15, "all", 1)
# setRestrictionDynamic(7, 15, "all", 1)
# setRestrictionDynamic(8, 15, "all", 1)
# setRestrictionDynamic(9, 15, "all", 1)
# setRestrictionDynamic(10, 15, "all", 1)
# setRestrictionDynamic(11, 15, "all", 1)
# setRestrictionDynamic(12, 15, "all", 1)
# setRestrictionDynamic(13, 15, "all", 1)
# setRestrictionDynamic(14, 15, "all", 1)
# setRestrictionDynamic(15, 15, "all", 1)
# setRestrictionDynamic(16, 15, "all", 1)
# setRestrictionDynamic(17, 15, "all", 1)
# setRestrictionDynamic(18, 15, "all", 1)
# setRestrictionDynamic(19, 15, "all", 1)
# setRestrictionDynamic(20, 15, "all", 1)
# setRestrictionDynamic(21, 15, "all", 1)
# setRestrictionDynamic(22, 15, "all", 1)
# setRestrictionDynamic(23, 15, "all", 1)
# setRestrictionDynamic(24, 15, "all", 1)
# setRestrictionDynamic(25, 15, "all", 1)

# # Period 2
# setRestrictionDynamic(1, 10, "all", 2)
# setRestrictionDynamic(2, 10, "all", 2)
# setRestrictionDynamic(3, 10, "all", 2)
# setRestrictionDynamic(4, 10, "all", 2)
# setRestrictionDynamic(5, 10, "all", 2)
# setRestrictionDynamic(6, 10, "all", 2)
# setRestrictionDynamic(7, 10, "all", 2)
# setRestrictionDynamic(8, 10, "all", 2)
# setRestrictionDynamic(9, 10, "all", 2)
# setRestrictionDynamic(10, 10, "all", 2)
# setRestrictionDynamic(11, 10, "all", 2)
# setRestrictionDynamic(12, 10, "all", 2)
# setRestrictionDynamic(13, 10, "all", 2)
# setRestrictionDynamic(14, 10, "all", 2)
# setRestrictionDynamic(15, 10, "all", 2)
# setRestrictionDynamic(16, 10, "all", 2)
# setRestrictionDynamic(17, 10, "all", 2)
# setRestrictionDynamic(18, 10, "all", 2)
# setRestrictionDynamic(19, 10, "all", 2)
# setRestrictionDynamic(20, 10, "all", 2)
# setRestrictionDynamic(21, 10, "all", 2)
# setRestrictionDynamic(22, 10, "all", 2)
# setRestrictionDynamic(23, 10, "all", 2)
# setRestrictionDynamic(24, 10, "all", 2)
# setRestrictionDynamic(25, 10, "all", 2)

# # Period 3
# setRestrictionDynamic(1, 5, "all", 3)
# setRestrictionDynamic(2, 5, "all", 3)
# setRestrictionDynamic(3, 5, "all", 3)
# setRestrictionDynamic(4, 5, "all", 3)
# setRestrictionDynamic(5, 5, "all", 3)
# setRestrictionDynamic(6, 5, "all", 3)
# setRestrictionDynamic(7, 5, "all", 3)
# setRestrictionDynamic(8, 5, "all", 3)
# setRestrictionDynamic(9, 5, "all", 3)
# setRestrictionDynamic(10, 5, "all", 3)
# setRestrictionDynamic(11, 5, "all", 3)
# setRestrictionDynamic(12, 5, "all", 3)
# setRestrictionDynamic(13, 5, "all", 3)
# setRestrictionDynamic(14, 5, "all", 3)
# setRestrictionDynamic(15, 5, "all", 3)
# setRestrictionDynamic(16, 5, "all", 3)
# setRestrictionDynamic(17, 5, "all", 3)
# setRestrictionDynamic(18, 5, "all", 3)
# setRestrictionDynamic(19, 5, "all", 3)
# setRestrictionDynamic(20, 5, "all", 3)
# setRestrictionDynamic(21, 5, "all", 3)
# setRestrictionDynamic(22, 5, "all", 3)
# setRestrictionDynamic(23, 5, "all", 3)
# setRestrictionDynamic(24, 5, "all", 3)
# setRestrictionDynamic(25, 5, "all", 3)

# # Period 4
# setRestrictionDynamic(1, 0, "all", 4)
# setRestrictionDynamic(2, 0, "all", 4)
# setRestrictionDynamic(3, 0, "all", 4)
# setRestrictionDynamic(4, 0, "all", 4)
# setRestrictionDynamic(5, 0, "all", 4)
# setRestrictionDynamic(6, 0, "all", 4)
# setRestrictionDynamic(7, 0, "all", 4)
# setRestrictionDynamic(8, 0, "all", 4)
# setRestrictionDynamic(9, 0, "all", 4)
# setRestrictionDynamic(10, 0, "all", 4)
# setRestrictionDynamic(11, 0, "all", 4)
# setRestrictionDynamic(12, 0, "all", 4)
# setRestrictionDynamic(13, 0, "all", 4)
# setRestrictionDynamic(14, 0, "all", 4)
# setRestrictionDynamic(15, 0, "all", 4)
# setRestrictionDynamic(16, 0, "all", 4)
# setRestrictionDynamic(17, 0, "all", 4)
# setRestrictionDynamic(18, 0, "all", 4)
# setRestrictionDynamic(19, 0, "all", 4)
# setRestrictionDynamic(20, 0, "all", 4)
# setRestrictionDynamic(21, 0, "all", 4)
# setRestrictionDynamic(22, 0, "all", 4)
# setRestrictionDynamic(23, 0, "all", 4)
# setRestrictionDynamic(24, 0, "all", 4)
# setRestrictionDynamic(25, 0, "all", 4)

# split 'all' in 2 constraints

# setRestrictionDynamic(2, 10, "CH2", 1)
# setRestrictionDynamic(2, 8, "CH2", 2)
# setRestrictionDynamic(2, 6, "CH2", 3)
# setRestrictionDynamic(2, 4, "CH2", 4)
# setRestrictionDynamic(2, 10, "LH2", 1)
# setRestrictionDynamic(2, 8, "LH2", 2)
# setRestrictionDynamic(2, 6, "LH2", 3)
# setRestrictionDynamic(2, 4, "LH2", 4)

# setRestrictionDynamic(1, 10, "CH2", 1)
# setRestrictionDynamic(1, 8, "CH2", 2)
# setRestrictionDynamic(1, 6, "CH2", 3)
# setRestrictionDynamic(1, 4, "CH2", 4)
# setRestrictionDynamic(1, 10, "LH2", 1)
# setRestrictionDynamic(1, 8, "LH2", 2)
# setRestrictionDynamic(1, 6, "LH2", 3)
# setRestrictionDynamic(1, 4, "LH2", 4)

# setRestrictionDynamic(2, 8, "CH2", 1)
# setRestrictionDynamic(2, 6, "CH2", 2)
# setRestrictionDynamic(2, 4, "CH2", 3)
# setRestrictionDynamic(2, 2, "CH2", 4)
# setRestrictionDynamic(2, 8, "LH2", 1)
# setRestrictionDynamic(2, 6, "LH2", 2)
# setRestrictionDynamic(2, 4, "LH2", 3)
# setRestrictionDynamic(2, 2, "LH2", 4)

# setRestrictionDynamic(1, 8, "CH2", 1)
# setRestrictionDynamic(1, 6, "CH2", 2)
# setRestrictionDynamic(1, 4, "CH2", 3)
# setRestrictionDynamic(1, 2, "CH2", 4)
# setRestrictionDynamic(1, 8, "LH2", 1)
# setRestrictionDynamic(1, 6, "LH2", 2)
# setRestrictionDynamic(1, 4, "LH2", 3)
# setRestrictionDynamic(1, 2, "LH2", 4)

# setRestrictionDynamic(1, 5, "all", 4)
# setRestrictionDynamic(1, 5, "all", 3)
# setRestrictionDynamic(1, 5, "all", 2)
# setRestrictionDynamic(1, 5, "all", 1)
