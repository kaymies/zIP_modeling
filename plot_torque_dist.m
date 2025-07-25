%% PLOT_TORQUE_DIST  Simulate best-fit models and plot torque distributions (Fig. 5)
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

    input_struct.lumped_params = model_struct.Input.lumped_params;
    beta = model_struct.BestFitParams.beta;
    omega = model_struct.BestFitParams.omega;
    input_struct.controller_params.R = alpha*[beta omega; omega 1/beta];  % alpha = 1e6, initialized in setup.m   
    input_struct.controller_params.Q = eye(4);
    input_struct.motorNoiseRatio = model_struct.BestFitParams.sigma_r;
    
    disp(' ')
    disp(['Simulating best-fit model of ' data_types{index_group} ' group.'])
    disp('This may take a few minutes . . .')
    disp(' ')

    output_struct = simulate_nonlinDIP(input_struct);

    torques_cmd = output_struct.torque_cmd;

    histogram2(torques_cmd(1,:),torques_cmd(2,:),'DisplayStyle','tile','EdgeColor','none');
    a = colorbar;clim([0, 2500]);
    a.Label.String = 'Number of occurences';
    a.FontSize = 16;
    box on;set(gca,'linewidth',2);
    axis square
    xlabel('\tau_a [Nm]');xlim([-.6,.6]);xticks(-0.6:0.2:0.6);
    ylabel('\tau_h [Nm]');ylim([-.6,.6]);yticks(-0.6:0.2:0.6);
    title(data_types{index_group});    
    
end