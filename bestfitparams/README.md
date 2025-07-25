# Best-fit-model parameters

Contains results of parameter search from search_params_bestfit.m,
of the models that best reproduce the intersection-point-height (zIP) data
from each data group (unimpaired, paretic, nonparetic).

## File content

Each file contains the following variables:

    BestFitParams - struct containing the parameter values (beta, omega, sigma_r) for the best-fit model for that group

    Input         - struct containing input to computeAnalyticIP.m, including fields:
                    f: frequency at which to compute zIP
                    lumped_params: lumped parameters of double-inverted-pendulum model
                    controller_params: controller parameters of model (Q & R matrices of LQR)
                    motorNoiseRatio: ankle-to-hip motor noise ratio, sigma_r

    NumSubjects   - number of subjects in data set

    RMSE          - root-mean-square error between model zIP and human average zIP

    SearchParams  - struct containing the parameter values (beta, omega, sigma_r) that were searched over

    SubjectType   - details on data group that the model was fit to