""" Params_case

The Parameters of the case study model.

All parameters are documented - see below.

Thesis OR&QL - EUR 2022
@author: Erik van der Heide

"""

# Some settings
import pandas as pd
import numpy as np
import math
pd.options.display.show_dimensions = False
numGrids = 25  # this still has to be hardcoded
numProds = 2  # this still has to be hardcoded

""" Choices """
printData = True  # True to print the data that is used
printVariables = True  # True to print the variables of the model
printBinaryVariables = False  # True to print the binary variables X, Y, Z
timeLimit = 10*60  # Time limit in seconds
printStatus = True  # True to print the status while running
runStaticModel = True  # Run the model without time periods
if runStaticModel:
    solveTimePeriod = 4  # Solve time period 1, 2, 3 or 4
distanceWithinGrid = 5  # 0 to assume no transport cost in grid, otherwise take avg. transport distance at location
onlyOneImportProduct = True  # True to only be able to import one hydrogen product in a grid
costImportance = 1.0  # between 0.0 and 1.0
CO2Importance = 0.0  # between 0.0 and 1.0
printTests = False  # True to print some tests
CCSCostInput = 25  # $/ton CO2 produced; variable between 5 and 55
optimizeCI = True  # True to optimize carbon intensity variables
includeCCS = True  # True to include CCS as option
multipleSizes = True  # True to include plants with sizes small, medium and large
avgCapSmall = math.floor(20+(99-20)/2)  # Where TCC is based on for small units
avgCapMedium = math.floor(100+(499-100)/2)  # Where TCC is based on for medium units
avgCapLarge = math.floor(500+(1000-500)/2)  # Where TCC is based on for large units

""" Carbon restrictions """

# Total CO2 emissions
TotalCO2Max = -1  # 21915.1

# Individual CI restrictions (CI^{max}_{ig} and CI^{max}_{g})
maxCI = -1 * np.ones(shape=((numProds+1), numGrids))


def setRestriction(gridNr, restrictionValue, productType):
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
        maxCI[(0, gridNr-1)] = restrictionValue
    elif productType.upper() == "LH2":
        maxCI[(1, gridNr-1)] = restrictionValue
    elif productType.upper() == "ALL":
        maxCI[(2, gridNr-1)] = restrictionValue


# Set the actual restrictions
#setRestriction(1, 10, "all")
#setRestriction(2, 3, "all")

