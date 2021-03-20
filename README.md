# 2021_CadhMultisite
 scripts related to our paper on Cadherin11

These scripts are supposed to help understanding and reproducing the results we are currently publishing
(will add reference here)

List of files

*** goglobalfit.py ***  
Main file to perform the global fitting process. The file is divided in sections to allow stepwise
execution of each stage.  
*** hkimports2.py ***  
This file imports all modules which are not included in this package. We intend to keep all those imports
together so that it is easier to keep track of all required modules.  
*** iofunctions.py ***  
Functions related to import and export of data  
*** datamodelfunctions.py ***  
Functions related to the datamodels used for handling experimental data and parameters  
*** globalfitunctions.py ***  
Functions related to the execution of the global fitting process  
*** preprocessing.py ***  
Functions related to transferring / pre-processing / converting experimental data to our data model  
*** RDmath.py ***  
Approximations and functions describing relaxation dispersion experiment. Related to our previous
theoretical publications.  

*** exptl_data ***  
experimental data. The path "path2020" to this folder has to be set in goglobalfitfuntions.py and goglobalfit.py  