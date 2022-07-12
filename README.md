# MSc Thesis Operations Research & Quantitative Logistics
### Carbon Intensity Modeling in Energy Value Chain Optimization
*Erik van der Heide. Erasmus University Rotterdam 2022.*

This repository contains all files used for my master thesis "Carbon Intensity Modeling in Energy Value Chain Optimization".

The code split into into two parts. The first part considers the single-resource, multi-transport model and the second part the multi-resource, multi-transport model. All modeling is done in the Python language using the Gurobi solver.

## Single-resource, multi-transport model (```Model SRMT```)
The model which considers only one resource (product) is split into 5 separate example setups. This model considers *pools*: all materials that come into an intermediate pool are linearly blended. This model was a build-up towards the multi-resource model, no results were stated in the report. It containts of the following files:
* ```Model_SRMT.py```. Serves at the **main** for this model. Containts the model which can run all 5 data setups. Also contains options to run the model. 
* ```Data_SRMT_1.py```. Consists of 3 source nodes, 1 pool and 3 target nodes.
* ```Data_SRMT_2.py```. Consists of 3 sources, 2 pools (in 1 layer) and 3 target nodes.
* ```Data_SRMT_3.py```. Consists of 3 sources, 4 pools (in 2 layers) and 3 target nodes.
* ```Data_SRMT_4.py```. Consists of 3 sources, 2 pools (in 1 layer), 3 target nodes and one extra model of transport.
* ```Data_SRMT_5.py```. Consists of 3 sources, 2 pools (in 1 layer), 3 target nodes and CCS investment options defined on the second pool.

The topology of the five experiments:
![5setups (2)](https://user-images.githubusercontent.com/75078739/178499245-e3a81d48-622b-4cdb-8bd6-9de10eb12f7d.jpg)


## Multi-resource, multi-transport model (```Model MRMT```)
The model which is described extensively in the report. It contains the following files:
* ```Main_MRMT.py```. Serves as the **main** of this model. To extensively test this model, this main is able to perform single or multiple runs of the model, as well as making insightful plots of the results.
* ```Params_MRMT.py```. Specifies several settings of single runs, such as max. running time and carbon restrictions.
* ```Model_MRMT.py```. Builds the main model, as well as carbon intensity restrictions and calculations. 
* ```Data_SRMT.py```. Consists of the data setup for the model (see topology below).
* ```Data_SRMT_Scaled.py```. Consists of data of the same type of topology as ```Data_SRMT.py```, but scaled up 5 times and random disturbances in the data are applied.

The topology of the MRMT model that is considered in the data: 
![final_setup](https://user-images.githubusercontent.com/75078739/176399010-dc961d8e-abca-42a3-b799-a092e74f72fa.jpg)
