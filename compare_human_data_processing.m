%% COMPARE_HUMAN_DATA_PROCESSING  Process human force plate data to obtain
%                      intersection-point heights for each train of each
%                      subject, and compare the different methods of doing
%                      so
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu), 
%          Rika Sugimoto Dimitrova (rikasd@mit.edu),
%          Kreg G. Gruben
% 2025-07-22

run setup

warning('off','MATLAB:table:ModifiedAndSavedVarnames')

ind_time = 1;
ind_Fx_d = 2;
ind_Fz_d = 3;
ind_CP_d = 4;
ind_Fx_nd = 5;
ind_Fz_nd = 6;
ind_CP_nd = 7;

params_ip.method = 'gruben_noSDnorm';
params_ip.f_int = 0.5;

VAF_thresh = 0.8;

data_zIP_gruben = readtable('data\zIP lab-frame.txt');

f_Hz_gruben = data_zIP_gruben.lab_frame;

%% Controls

for iSubj = 1:22
figure('Position',[50,0,1100,640]);colors = get(gca,'colororder');

for itrial = 1:3

clear data_ip

data_raw = ...
    importdata(['data\Control ' num2str(iSubj-1) ' trial ' num2str(itrial-1) '.txt']);
data_heights = importdata('data\heights.txt');

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(12 + iSubj)];

% Dominant leg
zIP_gruben = data_zIP_gruben.(['D' num2str(iSubj) '_' num2str(itrial)])/0.56;

data_ip.COP = data_raw.data(:,ind_CP_d)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_d)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_d)';

[f,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;

subplot(221);hold on;
plot(f_Hz_gruben,zIP_gruben,'-','Color',colors(itrial,:))
plot(f,zIP,'--','Color',colors(itrial,:));

subplot(223);hold on;
plot(f,zIP-zIP_gruben,'-','Color',colors(itrial,:));

% Non-dominant leg
zIP_gruben = data_zIP_gruben.(['ND' num2str(iSubj) '_' num2str(itrial)])/0.56;

data_ip.COP = data_raw.data(:,ind_CP_nd)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_nd)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_nd)';

[f,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;

subplot(222);hold on;
plot(f_Hz_gruben,zIP_gruben,'-','Color',colors(itrial,:))
plot(f,zIP,'--','Color',colors(itrial,:));

subplot(224);hold on;
plot(f,zIP-zIP_gruben,'-','Color',colors(itrial,:));

end

subplot(221);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
title('Dominant Leg');

subplot(223);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

subplot(222);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
legend('Original','non-norm, 2D pca');
title('Non-dominant Leg');

subplot(224);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);
legend('Trial 1','Trial 2','Trial 3');

sgtitle(['Control Subject ' num2str(iSubj)]);

end

%% Stroke patients
itrial = 1;

for iSubj = 1:12

clear data_ip

figure('Position',[50,0,1100,640]);colors = get(gca,'colororder');

data_raw = ...
    importdata(['data\Stroke ' num2str(iSubj-1) '.txt']);
data_heights = importdata('data\heights.txt');

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(iSubj)];

% Paretic leg
zIP_gruben = data_zIP_gruben.(['P' num2str(iSubj)])/0.56;

data_ip.COP = data_raw.data(:,ind_CP_d)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_d)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_d)';

[f,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;

subplot(221);hold on;
plot(f_Hz_gruben,zIP_gruben,'-','Color',colors(itrial,:))
plot(f,zIP,'--','Color',colors(itrial,:));

subplot(223);hold on;
plot(f,zIP-zIP_gruben,'-','Color',colors(itrial,:));

% Non-paretic leg
zIP_gruben = data_zIP_gruben.(['NP' num2str(iSubj)])/0.56;

data_ip.COP = data_raw.data(:,ind_CP_nd)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_nd)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_nd)';

[f,zIP,VAF] = getZIPfromData(data_ip,params_ip);
zIP(VAF < VAF_thresh) = NaN;

subplot(222);hold on;
plot(f_Hz_gruben,zIP_gruben,'-','Color',colors(itrial,:))
plot(f,zIP,'--','Color',colors(itrial,:));

subplot(224);hold on;
plot(f,zIP-zIP_gruben,'-','Color',colors(itrial,:));

subplot(221);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
title('Paretic Leg');

subplot(223);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

subplot(222);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
legend('Original','non-norm, 2D pca');
title('Non-paretic Leg');

subplot(224);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

sgtitle(['Post-Stroke Subject ' num2str(iSubj)]);

end

%%
%% Compare rectangular window vs. Hann window

%% Controls

for iSubj = 1:22
figure('Position',[50,0,1100,640]);colors = get(gca,'colororder');

for itrial = 1:3

clear data_ip

data_raw = ...
    importdata(['data\Control ' num2str(iSubj-1) ' trial ' num2str(itrial-1) '.txt']);
data_heights = importdata('data\heights.txt');

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(12 + iSubj)];

% Dominant leg
data_ip.COP = data_raw.data(:,ind_CP_d)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_d)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_d)';

params_ip.method = 'gruben_noSDnorm';
[f,zIP_rect,VAF] = getZIPfromData(data_ip,params_ip);
zIP_rect(VAF < VAF_thresh) = NaN;
params_ip.method = 'bpf';
[~,zIP_hann,VAF] = getZIPfromData(data_ip,params_ip);
zIP_hann(VAF < VAF_thresh) = NaN;

subplot(221);hold on;
plot(f,zIP_hann,'-','Color',colors(itrial,:))
plot(f,zIP_rect,'--','Color',colors(itrial,:));

subplot(223);hold on;
plot(f,zIP_rect-zIP_hann,'-','Color',colors(itrial,:));

% Non-dominant leg
data_ip.COP = data_raw.data(:,ind_CP_nd)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_nd)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_nd)';

params_ip.method = 'gruben_noSDnorm';
[f,zIP_rect,VAF] = getZIPfromData(data_ip,params_ip);
zIP_rect(VAF < VAF_thresh) = NaN;
params_ip.method = 'bpf';
[~,zIP_hann,VAF] = getZIPfromData(data_ip,params_ip);
zIP_hann(VAF < VAF_thresh) = NaN;

subplot(222);hold on;
plot(f,zIP_hann,'-','Color',colors(itrial,:))
plot(f,zIP_rect,'--','Color',colors(itrial,:));

subplot(224);hold on;
plot(f,zIP_rect-zIP_hann,'-','Color',colors(itrial,:));

end

subplot(221);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
title('Dominant Leg');

subplot(223);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

subplot(222);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
legend('hann','rect');
title('Non-dominant Leg');

subplot(224);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);
legend('Trial 1','Trial 2','Trial 3');

sgtitle(['Control Subject ' num2str(iSubj)]);

end

%% Stroke patients
itrial = 1;

for iSubj = 1:12

clear data_ip

figure('Position',[50,0,1100,640]);colors = get(gca,'colororder');

data_raw = ...
    importdata(['data\Stroke ' num2str(iSubj-1) '.txt']);
data_heights = importdata('data\heights.txt');

time_s = data_raw.data(:,ind_time)';
data_ip.sampFreq_Hz = 1/(time_s(2) - time_s(1));
data_ip.COM = zeros([2,length(time_s)]) + [0; 0.56*data_heights.data(iSubj)];

% Paretic leg
data_ip.COP = data_raw.data(:,ind_CP_d)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_d)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_d)';

params_ip.method = 'gruben_noSDnorm';
[f,zIP_rect,VAF] = getZIPfromData(data_ip,params_ip);
zIP_rect(VAF < VAF_thresh) = NaN;
params_ip.method = 'bpf';
[~,zIP_hann,VAF] = getZIPfromData(data_ip,params_ip);
zIP_hann(VAF < VAF_thresh) = NaN;

subplot(221);hold on;
plot(f,zIP_hann,'-','Color',colors(itrial,:))
plot(f,zIP_rect,'--','Color',colors(itrial,:));

subplot(223);hold on;
plot(f,zIP_rect-zIP_hann,'-','Color',colors(itrial,:));

% Non-paretic leg
data_ip.COP = data_raw.data(:,ind_CP_nd)';
data_ip.FootForce(1,:) = -data_raw.data(:,ind_Fx_nd)';
data_ip.FootForce(2,:) =  data_raw.data(:,ind_Fz_nd)';

params_ip.method = 'gruben_noSDnorm';
[f,zIP_rect,VAF] = getZIPfromData(data_ip,params_ip);
zIP_rect(VAF < VAF_thresh) = NaN;
params_ip.method = 'bpf';
[~,zIP_hann,VAF] = getZIPfromData(data_ip,params_ip);
zIP_hann(VAF < VAF_thresh) = NaN;

subplot(222);hold on;
plot(f,zIP_hann,'-','Color',colors(itrial,:))
plot(f,zIP_rect,'--','Color',colors(itrial,:));

subplot(224);hold on;
plot(f,zIP_rect-zIP_hann,'-','Color',colors(itrial,:));

subplot(221);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
title('Paretic Leg');

subplot(223);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

subplot(222);
yline(1,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM}');
xlim([0,6]);
ylim([0,3]);
legend('hann','rect');
title('Non-paretic Leg');

subplot(224);
yline(0,'k--');
xlabel('Frequency (Hz)');
ylabel('z_{IP}/z_{COM} difference');
xlim([0,6]);
ylim([-0.5,0.5]);

sgtitle(['Post-Stroke Subject ' num2str(iSubj)]);

end
