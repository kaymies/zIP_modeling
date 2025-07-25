function [f_zIP, zIP_ratio, VAF] = getZIPfromData(data,params)
%% GETZIPFROMDATA Get zIP over freq from data
% Obtain intersection-point height (zIP) at different frequencies from
% standing data
%
% Inputs:
% data = a struct containing a single trial of data, or
%        a cell with M structs (one struct per trial, for M trials of data)
%        containing the following info:
%      - COM = 2xN COM position in the 2D plane of interest
%              with the first row corresponding to the horizontal direction
%              and the second row to the vertical direction.
%      - COP = 1xN COP position in the 2D plane of interest
%      - FootForce = 2xN Foot-force vector
%              with first row corresponding to the horizontal component
%              and the second row to the vertical component.
%      - sampFreq_Hz = sampling frequency of data, in Hz
% params = struct containing parameters for zIP analysis
%      - method = string that specifies the method for computing zIP from
%                 experimental data; choose from:
%                 'cpsd' (default) = co-spectral density method from
%                              Sugimoto-Dimitrova et al. 2024
%                 'bpf'      = band-pass-filter (BPF) method from
%                              Boehm et al. 2019                 
%                 'cpsd2bpf' = smeared version of cpsd, to reproduce the
%                              bpf results based on cpsd results averaged
%                              over a bandwidth (weighted average based on 
%                              the bpf magnitude response)
%      - window_size = CPSD window size (default: 2^9)
%      - detrend_option = 1: no detrend, 2: zeromean, 3: linear detrend (default)
%      - window_taper_option = 1: rectangular, 2: Hann (default), 3: Hamming
%      - f_int = BPF bandwidth (default: 0.2)
%
% Outputs:
% f_zIP     = frequencies at which zIP is evaluated
% zIP_ratio = zIP/zCOM; intersection-point height normalized by COM height
%
% Rika Sugimoto Dimitrova
% 2024-02-06
% Last updated: 2025-07-23

% CPSD detrend settings
iRAW = 1;
iZEROMEAN = 2;
iDETREND = 3;

% CPSD window taper settings
iRECTWIN = 1;
iHANNWIN = 2;
iHAMMWIN = 3;

iX = 1; % index of horizontal direction
iZ = 2; % index of vertical direction

% Default parameters
method = 'cpsd';
window_size = 2^9;
detrend_option = iDETREND;
window_taper_option = iHANNWIN;
f_int = 0.2; % BPF frequency interval

if iscell(data)
    data_cell = data;
elseif isstruct(data)
    data_cell{1} = data;
else
    ME = MException('getZIPfromData:incorrectInputType', ...
        'data argument should be of type cell or struct');
    throw(ME)
end

if nargin > 1 && isfield(params,'method')
    method = params.method;
end

if strcmp(method,'bpf') || strcmp(method,'cpsd2bpf')
    window_size = length(data_cell{1}.COP);   % Use full data length
elseif nargin > 1 && isfield(params,'window_size')
    window_size = params.window_size;
end

if nargin > 1 && isfield(params,'detrend_option')
    detrend_option = params.detrend_option;
end

if nargin > 1 && isfield(params,'window_taper_option')
    window_taper_option = params.window_taper_option;
end

if nargin > 1 && isfield(params,'f_int')
    f_int = params.f_int;
end

switch window_taper_option
    case iRECTWIN
        windowtaper = rectwin(window_size);
    case iHANNWIN
        windowtaper = hann(window_size);
    case iHAMMWIN
        windowtaper = hamming(window_size);
end

nfft = window_size*2;
noverlap = 0.5*window_size;

% BPF frequency bins
f_i = f_int/2; f_f = 6 - f_i; % starting freq, frequency range
f_bpf = (f_i:f_int:f_f)'; % frequency bins
        
Fs_Hz = data_cell{1}.sampFreq_Hz;
dt = 1/Fs_Hz;

switch method
    case 'cpsd'

        for iTrial = 1:length(data_cell)
            COM_z = data_cell{iTrial}.COM(iZ,:);    % vertical CoM
            COP_x = data_cell{iTrial}.COP(iX,:);      % 1D CoP
            Fx = data_cell{iTrial}.FootForce(iX,:); % horizontal ground reaction force
            Fz = data_cell{iTrial}.FootForce(iZ,:); % vertical ground reaction force
            theta_F = -Fx./Fz;
        
            switch detrend_option
                case iRAW
                    y = [theta_F' COP_x']; 
                case iZEROMEAN
                    y = [detrend(theta_F,0)' detrend(COP_x,0)']; 
                case iDETREND
                    y = [detrend(theta_F,1)' detrend(COP_x,1)']; 
            end 
            
            [Gyy_temp, f_zIP] = ...
                cpsd_custom(y,y,windowtaper,noverlap,nfft,'mimo',Fs_Hz,...
                detrend_option);
    
            Gyy11(:,iTrial) = Gyy_temp(:,1,1);
            Gyy12(:,iTrial) = Gyy_temp(:,1,2);
            Gyy21(:,iTrial) = Gyy_temp(:,2,1);
            Gyy22(:,iTrial) = Gyy_temp(:,2,2);
        end
            
        zIP_ratio = zeros(window_size+1,1);
        VAF = zeros(window_size+1,1);
        for i = 1:(window_size+1)
            Gyy(1,1) = mean(Gyy11(i,:),2);
            Gyy(1,2) = mean(Gyy12(i,:),2);
            Gyy(2,1) = mean(Gyy21(i,:),2);
            Gyy(2,2) = mean(Gyy22(i,:),2);
            [V,D] = eig(real(Gyy(:,:)));
            D_vec = diag(D);
            [~,i_max] = max(D_vec);
            zIP_ratio(i) = V(2,i_max)/V(1,i_max) / mean(COM_z);
            VAF(i) = D_vec(i_max)/sum(D_vec);
        end

    case 'cpsd2bpf'

        for iTrial = 1:length(data_cell)
            COM_z = data_cell{iTrial}.COM(iZ,:);
            COP_x = data_cell{iTrial}.COP(iX,:);
            Fx = data_cell{iTrial}.FootForce(iX,:);
            Fz = data_cell{iTrial}.FootForce(iZ,:);
            theta_F = -Fx./Fz;

            switch detrend_option
                case iRAW
                    y = [theta_F' COP_x']; 
                case iZEROMEAN
                    y = [detrend(theta_F,0)' detrend(COP_x,0)']; 
                case iDETREND
                    y = [detrend(theta_F,1)' detrend(COP_x,1)']; 
            end 

            [Gyy, f_cpsd] = ...
                cpsd_custom(y,y,windowtaper,noverlap,nfft,'mimo',Fs_Hz,...
                detrend_option);

            N_smeared = length(f_bpf);
            for i = 1:N_smeared
                [B,A] = butter(2,[f_bpf(i) f_bpf(i)+f_int]/(1/dt/2));
                h = freqz(B,A,f_cpsd,Fs_Hz);
                Coyy_smeared_temp(i,:,:) = ... % weighted average based on BPF mag. response
                    sum((h.*conj(h)).^2.*real(Gyy))./sum((h.*conj(h)).^2);
            end % i

            Coyy_smeared11(:,iTrial) = Coyy_smeared_temp(:,1,1);
            Coyy_smeared12(:,iTrial) = Coyy_smeared_temp(:,1,2);
            Coyy_smeared21(:,iTrial) = Coyy_smeared_temp(:,2,1);
            Coyy_smeared22(:,iTrial) = Coyy_smeared_temp(:,2,2);
        end
            
        zIP_ratio = zeros(N_smeared,1);
        VAF = zeros(N_smeared,1);
        for i = 1:N_smeared
            Coyy_smeared(1,1) = mean(Coyy_smeared11(i,:),2);
            Coyy_smeared(1,2) = mean(Coyy_smeared12(i,:),2);
            Coyy_smeared(2,1) = mean(Coyy_smeared21(i,:),2);
            Coyy_smeared(2,2) = mean(Coyy_smeared22(i,:),2);
            [V,D] = eig(real(Coyy_smeared(:,:)));
            D_vec = diag(D);
            [~,i_max] = max(D_vec);
            zIP_ratio(i) = V(2,i_max)/V(1,i_max) / mean(COM_z);
            VAF(i) = D_vec(i_max)/sum(D_vec);
        end % i

        f_zIP = f_bpf + f_int/2;

    case 'bpf'
        zIP = zeros(length(f_bpf),1);
        VAF = zeros(length(f_bpf),1);
        
        for iTrial = 1:length(data_cell)
            COM_z = data_cell{iTrial}.COM(iZ,:);
            COP_x = data_cell{iTrial}.COP(iX,:);
            Fx = data_cell{iTrial}.FootForce(iX,:);
            Fz = data_cell{iTrial}.FootForce(iZ,:);
            theta_F = -Fx./Fz;
    
            cop_tapered = detrend(COP_x,1)'.*windowtaper; 
            theta_tapered = detrend(theta_F,1)'.*windowtaper;
            cop_padded = [ zeros(floor(length(cop_tapered)/2),1); 
                   cop_tapered;
                   zeros(floor(length(cop_tapered)/2),1)];
            theta_padded = [ zeros(floor(length(theta_tapered)/2),1); 
                   theta_tapered;
                   zeros(floor(length(theta_tapered)/2),1)];
            for f = 1:length(f_bpf)
                [B,A] = butter(2,[f_bpf(f) f_bpf(f)+f_int]/(1/dt/2)); % BPF
                COP_f(f,iTrial,:) = filtfilt(B,A,cop_padded);
                theta_f(f,iTrial,:) = filtfilt(B,A,theta_padded);
            end % f

            COM_z_allTrials = mean(COM_z);
        end % iTrials
        
        for f = 1:length(f_bpf)
            COP_f_allTrials = [];
            theta_f_allTrials = [];
            for iTrial = 1:length(data_cell)
                COP_f_allTrials = [COP_f_allTrials; squeeze(COP_f(f,iTrial,:))];
                theta_f_allTrials = [theta_f_allTrials; squeeze(theta_f(f,iTrial,:))];
            end % iTrials
            [coeff,~,latent] = pca([theta_f_allTrials COP_f_allTrials]);
            zIP(f) = coeff(2,1)/coeff(1,1);  
            VAF(f) = latent(1)/sum(latent);
        end % f

        zIP_ratio = zIP/mean(COM_z_allTrials);
        f_zIP = f_bpf + f_int/2;

    case 'gruben_SDnorm'

        rm_frac = 0.06;
        zIP = zeros(length(f_bpf),1);
        VAF = zeros(length(f_bpf),1);

        iTrial = 1;

        COM_z = data_cell{iTrial}.COM(iZ,:);
        COP_x = data_cell{iTrial}.COP(iX,:);
        Fx = data_cell{iTrial}.FootForce(iX,:);
        Fz = data_cell{iTrial}.FootForce(iZ,:);
        theta_F = -Fx./Fz;

        cop_zm = detrend(COP_x,0)'; 
        theta_zm = detrend(theta_F,0)';
        for f = 1:length(f_bpf)
            [B,A] = butter(2,[f_bpf(f) f_bpf(f)+f_int]/(1/dt/2)); % BPF
            COP_f_temp = filtfilt(B,A,cop_zm);
            theta_f_temp = filtfilt(B,A,theta_zm);

            N = length(COP_f_temp);
            index_start = round(rm_frac*N);
            index_end = N - index_start;
            COP_f(f,:) = COP_f_temp(index_start:index_end);
            theta_f(f,:) = theta_f_temp(index_start:index_end);

            SD_COP_f(f) = std(COP_f(f,:));
            SD_theta_f(f) = std(theta_f(f,:));
            COP_f_norm(f,:) = COP_f(f,:)/SD_COP_f(f);
            theta_f_norm(f,:) = theta_f(f,:)/SD_theta_f(f);
            [coeff,~,latent] = pca([theta_f_norm(f,:)' COP_f_norm(f,:)']);
            % y = detrend([theta_f_norm(f,:)' COP_f_norm(f,:)'],0);
            % A = y'*y;
            % A1 = [A(1:2); A];
            % A2 = [A1(:,1) A1];
            % [V,D] = eig(A);
            % D_vec = diag(D);
            % [D_max,i_max] = max(abs(D_vec));
            % zIP(f) = V(2,i_max)/V(1,i_max)*SD_COP_f(f)/SD_theta_f(f);
            % VAF(f) = D_max/sum(D_vec);
            zIP(f) = coeff(2,1)/coeff(1,1)*SD_COP_f(f)/SD_theta_f(f);
            % zIP(f) = SD_COP_f/SD_theta_f;
            VAF(f) = latent(1)/sum(latent);
        end % f

        zIP_ratio = zIP/mean(COM_z);
        f_zIP = f_bpf + f_int/2;

    case 'gruben_noSDnorm'

        rm_frac = 0.06;
        zIP = zeros(length(f_bpf),1);
        VAF = zeros(length(f_bpf),1);

        iTrial = 1;

        COM_z = data_cell{iTrial}.COM(iZ,:);
        COP_x = data_cell{iTrial}.COP(iX,:);
        Fx = detrend(data_cell{iTrial}.FootForce(iX,:),0);
        Fz = data_cell{iTrial}.FootForce(iZ,:);
        theta_F = -Fx./Fz;

        cop_zm = detrend(COP_x,0)'; 
        theta_zm = detrend(theta_F,0)';
        for f = 1:length(f_bpf)
            [B,A] = butter(2,[f_bpf(f) f_bpf(f)+f_int]/(1/dt/2)); % BPF
            COP_f_temp = filtfilt(B,A,cop_zm);
            theta_f_temp = filtfilt(B,A,theta_zm);

            N = length(COP_f_temp);
            index_start = round(rm_frac*N);
            index_end = N - index_start;
            COP_f(f,:) = COP_f_temp(index_start:index_end);
            theta_f(f,:) = theta_f_temp(index_start:index_end);

            [coeff,~,latent] = pca([theta_f(f,:)' COP_f(f,:)']);
            zIP(f) = coeff(2,1)/coeff(1,1);
            VAF(f) = latent(1)/sum(latent);
        end % f

        zIP_ratio = zIP/mean(COM_z);
        f_zIP = f_bpf + f_int/2;

end % switch method

end