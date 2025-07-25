# zIP_modeling
Code for modeling the intersection point of human quiet balance and quantifying control strategy

## Quick Start Guide:
1. Ensure that you have MATLAB R2024a or later installed, along with the following toolboxes:
    > Control System Toolbox, Signal Processing Toolbox, Statistics and Machine Learning Toolbox.
2. Navigate to the zIP_modeling folder in MATLAB and run run_main.m.

## Organization of the repository:

    \bestfitparams\   contains best-fit-model parameters, from parameter 
                      search via search_params_bestfit.m

    \data\            contains pre-processed human zIP data from 
                      dos Santos et al. (2017)

    \model\           contains functions related to the balance model and
                      simulation

    \visuals\         contains functions used for generating figures

    \zIP\             contains functions for computing zIP from data 
                      and from models

## Acknowledgements:
The simulation code is largely based on [Shiozawa et al. (2021)](https://doi.org/10.1186/s12984-021-00907-2), while the [zIP code](https://github.com/rikasd/zIP_spectralmethod) is based on [Sugimoto-Dimitrova et al. (2024)](https://doi.org/10.1152/jn.00084.2024).

## References:

dos Santos, D. A., Fukuchi, C. A., Fukuchi, R. K., & Duarte, M. (2017). A data set with kinematic and ground reaction forces of human balance. PeerJ, 2017(7), 1–17. https://doi.org/10.7717/peerj.3626

Shiozawa, K., Lee, J., Russo, M. et al. (2021). Frequency-dependent force direction elucidates neural control of balance. J NeuroEngineering Rehabil 18, 145. https://doi.org/10.1186/s12984-021-00907-2

Sugimoto-Dimitrova, R., Shiozawa, K., Gruben, K. G., & Hogan, N. (2024). Frequency-domain patterns in foot-force line-of-action: an emergent property of standing balance control. Journal of Neurophysiology. https://doi.org/10.1152/jn.00084.2024