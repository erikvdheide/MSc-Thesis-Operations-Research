# MSc Thesis Operations Research & Quantitative Logistics
### Carbon Intensity Modeling in Energy Value Chain Optimization
*Erik van der Heide. Erasmus University Rotterdam 2022.*

This repository contains all files used for my master thesis "Carbon Intensity Modeling in Energy Value Chain Optimization".

The code is split into 3 parts. The first part considers the single-resource, multi-transport model. The second part the multi-resource, multi-transport model. The third part considers the case study model (both static and dynamic versions). All modeling is done in the Python language using the Gurobi solver.


## Single-resource, multi-transport model (```Model SRMT```)
The model which considers only one resource (product) is split into 5 separate example setups. This model considers *pools*: all materials that come into an intermediate pool are linearly blended. This model was a build-up towards the multi-resource model, no results were stated in the report. It containts of the following files:
* ```Model_SRMT.py```. Run this file to execute the model. Containts the model formulations. In the top region of the page, you can fill in some parameters, including which of the 5 datasets you want to run. Run this file to execute the model. 
* ```Data_SRMT_1.py```. Consists of 3 source nodes, 1 pool and 3 target nodes.
* ```Data_SRMT_2.py```. Consists of 3 sources, 2 pools (in 1 layer) and 3 target nodes.
* ```Data_SRMT_3.py```. Consists of 3 sources, 4 pools (in 2 layers) and 3 target nodes.
* ```Data_SRMT_4.py```. Consists of 3 sources, 2 pools (in 1 layer), 3 target nodes and one extra model of transport.
* ```Data_SRMT_5.py```. Consists of 3 sources, 2 pools (in 1 layer), 3 target nodes and CCS investment options defined on the second pool.

The topology of the five experiments:
![5setups (2)](https://user-images.githubusercontent.com/75078739/178499245-e3a81d48-622b-4cdb-8bd6-9de10eb12f7d.jpg)


## Multi-resource, multi-transport model (```Model MRMT```)
The model which is described extensively in the report. It contains the following files:
* ```Main_MRMT.py```. Run this file to execute the model. Standard, the model stands on 1 model run given the parameters in ```Params_MRMT.py```. If you want to perform multiple runs and/or plot some results, you can uncomment the corresponding boxes. 
* ```Params_MRMT.py```. Specifies several settings of single runs, such as max. running time and carbon restrictions.
* ```Model_MRMT.py```. Contains all model components, as well as carbon intensity restrictions and calculations. 
* ```Data_SRMT.py```. Consists of the data setup for the model (see topology below).
* ```Data_SRMT_Scaled.py```. Consists of data of the same type of topology as ```Data_SRMT.py```, but scaled up 5 times and random disturbances in the data are applied.

The topology of the MRMT model that is considered in the data: 
![final_setup](https://user-images.githubusercontent.com/75078739/176399010-dc961d8e-abca-42a3-b799-a092e74f72fa.jpg)


## Hydrogen supply chain network case study model (```Model HSCN```)
The case study model for an optimal hydrogen configuration in the Netherlands. It contains the following files:
* ```Data_HSCN.py```. Contains the complete (hardcoded) dataset used for both the static and dynamic model (incl. different assumptions on CO2 data).
* ```Model_HSCN_dynamic.py```. Run this file to execute the *dynamic* model. Contains all dynamic model components.
* ```Model_HSCN_static.py```. Run this file to execute the *static* model. Contains all static model components.
* ```Params_HSCN.py```. Consists of the choices that are made in the models, as well as the option to put CI restrictions. Some important choices are which time period you run (in case of static model) and whether you allow to import multiple hydrogen products into a grid.
* ```Plots_HSCN.py```. Code that makes plots for using a starting solution or not (data hardcoded).

The main pathways of the HSCN model are given by:
![pathwaysH2](https://user-images.githubusercontent.com/75078739/188835006-22dd2b93-1954-41ac-b70c-abe42c34df6c.jpg)

