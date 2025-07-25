%% PLOT_STIFFNESS_ELLIPSE  Plot ellipse-circle visualization of stiffness (Fig. 4)
%
% Authors: Kaymie Shiozawa (kaymies@mit.edu) & Rika Sugimoto Dimitrova (rikasd@mit.edu)
% 2025-07-02

run setup

fh = figure;hold on;

%% Loop through best-fit models for 1) unimpaired, 2) paretic, 3) nonparetic
% Print best-fit gain matrices and K_ratio
disp(' ')
disp('Best-fit gain and K_ratio:')
for index_group = 1:n_groups
    disp(' ')
    disp([data_types{index_group} ' group:'])
    disp(' ')
    model_struct = load(sprintf('bestfitparams_%s.mat',data_filenames{index_group}));

    % Get best-fit-model parameters
    lumped_params = model_struct.Input.lumped_params;

    beta = model_struct.BestFitParams.beta;
    omega = 0;
    % omega = model_struct.BestFitParams.omega;
    R = alpha*[beta omega; omega 1/beta]; % alpha = 1e6, initialized in setup.m
    Q = eye(4);

    % Linearized mass matrix M, gravitational stiffness J_G
    q_eq = [0;0]; Dq_eq = [0;0];
    [M,~,~] = getDynamics_nonlinDIP(lumped_params,q_eq,Dq_eq);
    [~,~,J_G] = getJacobians_nonlinDIP(lumped_params,q_eq,Dq_eq);
    
    % Oopen-loop system state-space matrices
    A_ol = [zeros(2) eye(2); -M\J_G zeros(2)];
    B_ol = [zeros(2); M\eye(2)];
    
    % LQR controller gains
    K_x = lqr(A_ol, B_ol, Q, R);
    disp('Gain:')
    disp(K_x)
    K = K_x(1:2,1:2); % apparent stiffness matrix
    K_s = (K+K')/2; % symmetric stiffness
    K_a = (K-K')/2; % antisymmetric stiffness
    K_ratio = sqrt(abs(det(K_a)))/sqrt(abs(det(K_s)));
    disp('Stiffness:')
    disp(K)
    disp('Symmetric Stiffness:')
    disp(K_s)
    disp('Antisymmetric Stiffness:')
    disp(K_a)
    disp(['K_ratio: ' num2str(K_ratio)])
    disp(' ')

    lineSpecs.color = colors(index_group,:);
    lineSpecs.width = 2;
    lineSpecs.style = '-';
    [~, p_ellipse, p_circ] = drawStiffnessEllipse(K, fh.Number, lineSpecs);
    p_ellipse_all(index_group) = p_ellipse;
    p_circ_all(index_group) = p_circ;
    
end

box on;set(gca,'linewidth',2);
xlabel('Ankle [Nm/rad]');
ylabel('Hip [Nm/rad]');
title('Stiffness');
xticks(-500:500:500);
yticks(-500:500:500);
legend(p_ellipse_all, data_types, 'Location', 'eastoutside');