## Supporting Information for "Bayesian nonparametric analysis for the detection of spikes in noisy calcium imaging data"

Finite version of the Common Atom Model of Denti et al. (2020) based on the mixtures of finite mixtures of Frühwirth-Schnatter, Malsiner-Walli e Grün (2020).

Application to calcium imaging data from the Allen Brain Observatory (Allen  Institute  for  Brain  Science,  2016) for the detection of neurons' spiking activity.



The repository contains:
- **fCAM_calcium_imaging.cpp** : cpp code implementing the Gibbs sampler;
- **Simulation_study/** : folder containing the R scripts for generating the synthetic data (sim_data_scenario#.R) and an example code to run the Gibbs sampler (run.R);
- **Allen_Brain_data_analysis/** : folder containing the R script to run the analysis on a real series of calcium imaging data. The folder also contains a csv file (cell_data.csv) of a trace for a neuron, extracted from the Allen brain observatory database. The file group.Rdata contains the information about the times of the visual stimuli for the considered neuron.

#### References
- L. D'Angelo, A. Canale, Z. Yu and M. Guindani (2022). Bayesian nonparametric analysis for the detection of spikes in noisy calcium imaging data. *Biometrics*. doi = 10.1111/biom.13626
  
- Allen Institute MindScope Program (2016). Allen Brain Observatory – 2-photon visual coding
[dataset]. brain-map.org/explore/circuits observatory.brain-map.org/visualcoding .

- de Vries, S., Lecoq, J., Buice, M., Groblewski, P., Ocker, G., Oliver, M. et al. (2020). A large-scale standardized physiological survey reveals functional organization of the mouse visual cortex. *Nature neuroscience* **23**, 138–151.

- Denti, F., Camerlenghi, F., Guindani, M., and Mira, A. (2021). A common atoms model for the Bayesian nonparametric analysis of nested data. *Journal of the American Statistical Association*.

- Frühwirth-Schnatter, S., Malsiner-Walli, G., and Grün, B. (2021). Generalized mixtures of finite mixtures and telescoping sampling. *Bayesian Analysis* **16**, 1279 – 1307.
