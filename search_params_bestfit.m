%% SEARCH_PARAMS_BESTFIT  Run global parameter search
%                         to find model parameters that best fit human zIP
%                         (search may take over 30 min each to run)
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu) & Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2025-07-01

clear;clc;

run setup

%% Loop through data files
for index_group = 1:n_groups

    clear human_struct rmse_results params_bestfit
    
    %% Select data to fit:
    subject_type = data_filenames{index_group};
    filename_bestfit_results = ['bestfitparams_' subject_type];
    
    %% Load human data
    human_struct = load(sprintf('%s.mat',subject_type));
    zIP_human = human_struct.IPDataAverage;
    input_struct.f = human_struct.Frequency;
    
    %% Set model parameters
    totalMass_kg = human_struct.MeanMass_kg;
    totalHeight_m = human_struct.MeanHeight_m;  
    gender = 'M';
    plane = human_struct.Plane;
    pose = human_struct.Pose;
    
    input_struct.lumped_params = getLumpedParams_DIP(totalMass_kg, totalHeight_m, gender, plane, pose);
    input_struct.lumped_params.L_COM = 0.56*totalHeight_m;
    
    input_struct.controller_params.Q = eye(4);
    input_struct.controller_params.R = alpha*eye(2); % alpha = 1e6, initialized in setup.m
    
    %% Parameters to search through
    params.beta = 0.15:0.01:0.45;
    params.sigma_r = [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1];
    % params.omega = -1:0.1:1;
    
    % Create empty array to log errors
    rmse_results = zeros(length(params.beta),length(params.sigma_r));%,length(params.omega));
    
    %% Run Simulation for n trials across different parameter values
    disp(' ')
    disp(['Running parameter search for ' data_types{index_group} ' group.'])
    disp('This may take several minutes . . .')
    disp(' ')
    tic
    for b = 1:length(params.beta)
        for s = 1:length(params.sigma_r)
            % for o = 1:length(params.omega)
                input_struct.controller_params.R = ...
                    alpha*[params.beta(b) 0;...
                           0 1/params.beta(b)];            
                input_struct.motorNoiseRatio = params.sigma_r(s);
    
                % check that R is positive definite (as required by LQR)
                if all(eig(input_struct.controller_params.R)>0.0001) 
                    [~, zIP_model] = predictZIPfromModel(input_struct);
                else
                    zIP_model = NaN(size(input_struct.f));
                end
    
                rmse_results(b,s) = sqrt(mean((zIP_model-zIP_human).^2));
            % end
        end
    end
    toc
    
    %% Get best-fit parameters
    [v_min,loc_min] = min(abs(rmse_results(:)));
    [b_min,s_min] = ind2sub(size(rmse_results),loc_min);
    params_bestfit.beta = params.beta(b_min);
    params_bestfit.sigma_r = params.sigma_r(s_min);
    % params_bestfit.omega = params.omega(o_min);
    
    disp(' ')
    disp('Best-fit-model parameters:')
    disp(['beta: ' num2str(params_bestfit.beta)])
    % disp(['omega: ' num2str(params_bestfit.omega)])
    disp(['sigma_r: ' num2str(params_bestfit.sigma_r)])
    disp(['RMSE: ' num2str(v_min)])
    disp(' ')
    
    %% Save results
    bestfit_results_struct = struct(...
        'Input', input_struct,...
        'SubjectType', subject_type,...
        'NumSubjects', human_struct.NumSubjects,...
        'SearchParams', params,...
        'RMSE', rmse_results,...
        'BestFitParams', params_bestfit);
    save([pwd '\bestfitparams\' filename_bestfit_results], '-struct', 'bestfit_results_struct');
    
    disp(['Saved parameter search results at \data\' filename_bestfit_results '.mat'])

end