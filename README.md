## Supporting Information for "Bayesian nonparametric analysis for the detection of spikes in noisy calcium imaging data"

Finite version of the Common Atom Model of Denti et al. (2020) based on the mixtures of finite mixtures of Frühwirth-Schnatter, Malsiner-Walli e Grün (2020).

Application to calcium imaging data from the Allen Brain Observatory (Allen  Institute  for  Brain  Science,  2016) for the detection of neurons' spiking activity.



The repository contains:
- **fCAM_calcium_imaging.cpp** : cpp code implementing the Gibbs sampler;
- **Simulation_study/** : folder containing the R scripts for generating the synthetic data (sim_data_scenario#.R) and an example code to run the Gibbs sampler (run.R);
- **Allen_Brain_data_analysis/** : folder containing the R script to run the analysis on a real series of calcium imaging data. The folder also contains a csv file (cell_data.csv) of a trace for a neuron, extracted from the Allen brain observatory database. The file group.Rdata contains the information about the times of the visual stimuli for the considered neuron.

#### References
- Allen Institute for Brain Science (2016). Allen brain observatory.  http://observatory.brain-map.org/visualcoding

- F. Denti, F. Camerlenghi, M. Guindani and A. Mira (2020). A Common Atom Model for the Bayesian Nonparametric Analysis of Nested Data. *arXiv preprint arXiv:2008.07077*

- S. Frühwirth-Schnatter, G. Malsiner-Walli and B. Grün (2020). Generalized mixtures of finite mixtures and telescoping sampling. *arXiv preprint arXiv:2005.09918*
