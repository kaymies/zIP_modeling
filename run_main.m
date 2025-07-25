% Main script

clear; clc; close all;

run setup

%% Run parameter search for best-fit models 
% May take over 30 min to run; you may wish to skip running this section, 
% and use the search results already available in the "bestfitparams" 
% folder

run search_params_bestfit  % Run global parameter search to find model 
                           % parameters that best fit human zIP

%% Generate figures from paper:

run plot_human_zIP         % Fig. 1 - zIP from human data


run plot_model_zIP         % Fig. 3 - zIP from best-fit model


run plot_stiffness_ellipse % Fig. 4 - ellipse-circle representation of 
                           %          apparent stiffness of best-fit models