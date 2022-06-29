# MSc Thesis Operations Research & Quantitative Logistics
### Carbon Intensity Modeling in Energy Value Chain Optimization
*Erik van der Heide. Erasmus University Rotterdam 2022.*

This repository contains all files used for my master thesis "Carbon Intensity Modeling in Energy Value Chain Optimization".

The code split into into two parts. The first part considers the single-resource, multi-transport model and the second part the multi-resource, multi-transport model. All modeling is done in the Python language using the Gurobi solver.

## Single-resource, multi-transport model (SRMT Model)
The model which considers only one resource (product) is split into 5 separate example setups. This model considers *pools*: all materials that come into an intermediate pool are linearly blended. This model was a build-up towards the multi-resource model, no results were stated in the report. It containts of the following files:
* ```Model_SRMT.py```. Serves at the **main** for this model. Containts the model which can run all 5 data setups. Also contains all parameters and choices. 
* ```Data_SRMT_1.py```. Consists of 3 source nodes, 1 pool and 3 target nodes.
* ```Data_SRMT_2.py```. Consists of 3 sources, 2 pools (in 1 layer) and 3 target nodes.
* ```Data_SRMT_3.py```. Consists 
* ```Data_SRMT_4.py```.
* ```Data_SRMT_5.py```.


## Multi-resource, multi-transport model (MRMT Model)
The model which is described extensively in the report. It contains of the following files:
* ```Main_MRMT.py```. Serves as the **main** of this 
* ```Model_MRMT.py```.
* ```Data_SRMT.py```.

The topology of the model: 
![final_setup](https://user-images.githubusercontent.com/75078739/176399010-dc961d8e-abca-42a3-b799-a092e74f72fa.jpg)
