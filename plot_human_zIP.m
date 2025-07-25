%% PLOT_HUMAN_ZIP  Plot human intersection-point-height plot (Fig. 1)
%                  mean +/- standard error
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu) & Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2025-07-02

run setup

fh = figure;hold on;

%% Loop through data files for 1) unimpaired, 2) paretic, 3) nonparetic
for index_group = 1:n_groups

    human_struct = load(sprintf('%s.mat',data_filenames{index_group}));
    f_Hz = human_struct.Frequency;
    mean_zIP = human_struct.IPDataAverage;
    SE_zIP  = human_struct.StandardDeviation/sqrt(human_struct.NumSubjects);
    lineSpecs.color = colors(index_group,:);
    lineSpecs.width = 2;
    lineSpecs.alpha = 0.5;
    lineSpecs.style = '-';
    [p_mean, p_SE] = ...
        plotSDArea(f_Hz, mean_zIP, SE_zIP, fh.Number, 2, lineSpecs);
    p_mean_all(index_group) = p_mean;
    p_SE_all(index_group) = p_SE;
    
end

box on;set(gca,'linewidth',2);
xlabel('Frequency [Hz]');
ylabel({'Intersection Point Height','(Normalized by CoM height)'});
xticks(0:6);
legend(p_mean_all,data_types);