%% PROCESS_HUMAN_DATA  Process human force plate data to obtain
%                      intersection-point heights for each train of each
%                      subject
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu), 
%          Rika Sugimoto Dimitrova (rikasd@mit.edu),
%          Kreg G. Gruben
% 2025-07-22

run setup

ind_time = 1;

ind_Fx_d = 2;
ind_Fz_d = 3;
ind_CP_d = 4;
ind_Fx_nd = 5;
ind_Fz_nd = 6;
ind_CP_nd = 7;

ind_Fx_p = 2;
ind_Fz_p = 3;
ind_CP_p = 4;
ind_Fx_np = 5;
ind_Fz_np = 6;
ind_CP_np = 7;

params_ip.method = 'gruben_noSDnorm';
params_ip.f_int = 0.5;

VAF_thresh = 0.8;

data_heights = importdata('data\heights.txt');

%% Controls
dosSantos_data_struct = load('data\dosSantos2017_old.mat');
zIP_all = zeros(size(dosSantos_data_struct.IP));

N_subj = dosSantos_data_struct.NumSubjects;
N_trials = 3;

count = 1;

for iSubj = 1:N_subj 

for itrial = 1:N_trials

clear data_ip

data_raw = ...
    importdata(['data\Control ' num2str(iSubj-1) ' trial ' num2str(itrial-1) '.txt']);

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(12 + iSubj)];

% Dominant leg
data_ip.COP = data_raw.data(:,ind_CP_d)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_d)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_d)';

[~,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;
zIP_all(:,count) = zIP;

% Non-dominant leg
data_ip.COP = data_raw.data(:,ind_CP_nd)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_nd)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_nd)';

[~,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;
zIP_all(:, count + N_subj*N_trials) = zIP;

if iSubj == 5
    zIP_all(:, count + N_subj*N_trials) = NaN;
end

count = count + 1;

% if count == 40 || count == 14
%     pause(1);
% end

end

end

%%
zIP_mean = mean(zIP_all,2,'omitnan');
zIP_SD = std(zIP_all,[],2,'omitnan');

% Save results
controls_zIP_struct = struct(...
    'DataInfo', dosSantos_data_struct.DataInfo,...
    'Frequency', dosSantos_data_struct.Frequency,...
    'IP', zIP_all,...
    'IPDataAverage', zIP_mean,...
    'MeanHeight_m', dosSantos_data_struct.MeanHeight_m,...
    'MeanMass_kg', dosSantos_data_struct.MeanMass_kg,...
    'NumSubjects', N_subj,...
    'Plane', 'sgt',...
    'Pose', 'pose_I',...
    'StandardDeviation', zIP_SD,...
   'SubjectType', dosSantos_data_struct.SubjectType );
save([pwd '\data\dosSantos2017_old_v5_rmSubj5nd'], '-struct', 'controls_zIP_struct');

%% Post-stroke
Bartloff_p_struct = load('data\Bartloff2024_paretic.mat');
Bartloff_np_struct = load('data\Bartloff2024_nonparetic.mat');
zIP_p_all = zeros(size(Bartloff_p_struct.IP));
zIP_np_all = zeros(size(Bartloff_np_struct.IP));

N_subj = Bartloff_p_struct.NumSubjects;

count = 1;

for iSubj = 1:N_subj

clear data_ip

data_raw = importdata(['data\Stroke ' num2str(iSubj-1) '.txt']);

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(iSubj)];

% Paretic leg
data_ip.COP = data_raw.data(:,ind_CP_p)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_p)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_p)';

[~,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;
zIP_p_all(:,count) = zIP;

% Non-paretic leg
data_ip.COP = data_raw.data(:,ind_CP_np)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_np)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_np)';

[~,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;
zIP_np_all(:, count) = zIP;

count = count + 1;

end

%%
zIP_p_mean  = mean(zIP_p_all,2,'omitnan');
zIP_p_SD    = std(zIP_p_all,[],2,'omitnan');
zIP_np_mean = mean(zIP_np_all,2,'omitnan');
zIP_np_SD   = std(zIP_np_all,[],2,'omitnan');

% Save results
paretic_zIP_struct = struct(...
    'DataInfo', Bartloff_p_struct.DataInfo,...
    'Frequency', Bartloff_p_struct.Frequency,...
    'IP', zIP_p_all,...
    'IPDataAverage', zIP_p_mean,...
    'MeanHeight_m', Bartloff_p_struct.MeanHeight_m,...
    'MeanMass_kg', Bartloff_p_struct.MeanMass_kg,...
    'NumSubjects', N_subj,...
    'Plane', 'sgt',...
    'Pose', 'pose_I',...
    'StandardDeviation', zIP_p_SD,...
   'SubjectType', Bartloff_p_struct.SubjectType );
save([pwd '\data\Bartloff2024_paretic_v5'], '-struct', 'paretic_zIP_struct');

nonparetic_zIP_struct = struct(...
    'DataInfo', Bartloff_np_struct.DataInfo,...
    'Frequency', Bartloff_np_struct.Frequency,...
    'IP', zIP_np_all,...
    'IPDataAverage', zIP_np_mean,...
    'MeanHeight_m', Bartloff_np_struct.MeanHeight_m,...
    'MeanMass_kg', Bartloff_np_struct.MeanMass_kg,...
    'NumSubjects', N_subj,...
    'Plane', 'sgt',...
    'Pose', 'pose_I',...
    'StandardDeviation', zIP_np_SD,...
   'SubjectType', Bartloff_np_struct.SubjectType );
save([pwd '\data\Bartloff2024_nonparetic_v5'], '-struct', 'nonparetic_zIP_struct');