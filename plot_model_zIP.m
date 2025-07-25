%% PLOT_MODEL_ZIP  Plot best-fit-model intersection-point-height plot (Fig. 3)
%                  against human data (mean +/- standard deviation)
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu) & Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2025-07-02

run setup

%% Loop through best-fit models for 1) unimpaired, 2) paretic, 3) nonparetic
for index_group = 1:3

    figure;hold on;

    human_struct = load(sprintf('%s.mat',data_filenames{index_group}));
    model_struct = load(sprintf('bestfitparams_%s.mat',data_filenames{index_group}));

    f_Hz = human_struct.Frequency;
    mean_zIP_human = human_struct.IPDataAverage;
    SD_zIP  = human_struct.StandardDeviation;

    input_struct = model_struct.Input;
    beta = model_struct.BestFitParams.beta;
    omega = model_struct.BestFitParams.omega;
    input_struct.controller_params.R = alpha*[beta omega; omega 1/beta];  % alpha = 1e6, initialized in setup.m   
    input_struct.controller_params.Q = eye(4);
    input_struct.motorNoiseRatio = model_struct.BestFitParams.sigma_r;

    [~, zIP_model] = predictZIPfromModel(input_struct);

    errorbar(f_Hz, mean_zIP_human, SD_zIP, 'x', 'Color', [0.2 0.2 0.2], ...
        'LineWidth', 2, 'MarkerSize', 15);
    plot(f_Hz, zIP_model, 'o', 'Color', colors(index_group,:), 'MarkerSize', 10);

    box on;set(gca,'linewidth',2);
    xlabel('Frequency [Hz]');
    ylabel({'Intersection Point Height','(Normalized by CoM height)'});
    xticks(0:6);
    legend('Human Data','Modeled Data');
    title(data_types{index_group});
    
end