# HGSAE
This respository contains code and data to reproduce the results presented in Section 5 of the paper titled "A Hierarchical Gamma Prior for Modeling Random Effects in Small Area Estimation".

* data directory contains the state-level ([data1.txt](data/data1.txt)) and county-level ([data2.txt](data/data2.txt) data. The benchmark for the state-level estimation is in [data1_theta0.txt](data/data1_theta0.txt).
* [util_functions.R](util_functions.R) contains samplers for the proposed hierarchical gamma model and the global-local model proposed by Tang et al. (2018). It also provides functions for computing DIC and deviation measures.
* [state_code.R](state_code.R) is used to produce the state-level estimation results.
* [county_code.R](county_code.R) is used to produce the state-level estimation results.
