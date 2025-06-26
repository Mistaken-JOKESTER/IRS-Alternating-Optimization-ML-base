close all
clear all

M = 4;
K = 2;
N = 50;

snr = 0;
Pt = 10.^(snr/10); 
it = 15;

data_hd = csvread('data/Hd.csv'); 
data_hr = csvread('data/hr.csv'); 
data_G = csvread('data/G.csv');
data_omega = csvread('data/omega.csv');

% Initialize arrays to accumulate SNRs
snr_gen_all = zeros(1, it+1);
snr_seperate_all = zeros(1, it+1);
snr_togather_all = zeros(1, it+1);
snr_hybrid_all = zeros(1, it+1);
snr_nonml_all = zeros(1, it+1);

snr_gen_all_vector = [];
snr_seperate_all_vector = [];
snr_togather_all_vector = [];
snr_hybrid_all_vector = [];
snr_nonml_all_vector = [];


function snr_values = optimizeBeamforming(theta, W, M, K, N, Hd, Hr, G, omega)
    snr=0;
    Pt=10.^(snr/10); 
    it=100;
    Theta=diag(theta);
    W_old=zeros(M,K);
    beta=zeros(1,K);
    H=Hd+Hr*Theta*G;
    [ W,grt,f0 ] = init_W_2( H,M,K,Pt,omega, W);
    f1=0;
    W_span=W;
    W_last=W;
    [ ~,L_last ] = Proxlinear_beam_para( H,K,M,beta );
    flag=0;
    t_old=1;
    [ beta ] =upadte_beta( H,W,K,grt);

    snr_values = [f0];
    for con0=1:15
        [ Qx,qx,theta ] = surface_U_v_direct( W,Hd,Hr,Theta,G,N,K,grt,beta );
        theta_old=theta;

        U=-Qx;v=qx;
        x0=theta_old;
        phi0=angle(x0);
        grad=real((2*U*x0-2*v).*(-1j.*conj(x0)));
        dir=-grad;
        [ Ltheta ] = SCA_phi_step_para( -Qx,qx,N, theta );
        [ theta,t3, qqq ] = armijo_theta( Ltheta,dir,f0,phi0, grad,grt,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);

        [f1,grt,beta,W,W_span,t_old,L_last,H,Theta ] = fun_theta_package(grt,theta,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);

        f0=f1;     

        snr_values = [snr_values f0];   
    end
end

for zz = 1:100

    
    Hd = reconstructMatrix(K, M, data_hd(zz,:));
    G = reconstructMatrix(N, M, data_G(zz,:));
    Hr = reconstructMatrix(K, N, data_hr(zz,:)); 
    omega = reconstructMatrix(1, K, data_omega(zz,:));
    
    %% ************************************************************************
    %% Generated Data


    data_theta = csvread('data/theta.csv'); 
    data_W = csvread('data/W.csv'); 
    d_theta=reconstructMatrix(N, 1, data_theta(zz,:));
    d_W = reconstructMatrix(M, K, data_W(zz,:));
    snr_gen = optimizeBeamforming(d_theta, d_W, M, K, N, Hd, Hr, G, omega);


    %% ************************************************************************
    %% Theta sepearte and W seperate

    data_theta = csvread('data/only_theta_pred.csv'); 
    data_W = csvread('data/only_w_pred.csv'); 
    d_theta=reconstructMatrix(N,1, data_theta(zz,:));
    d_W = reconstructMatrix(M, K, data_W(zz,:));
    snr_seperate = optimizeBeamforming(d_theta, d_W, M, K, N, Hd, Hr, G, omega);


    %% ************************************************************************
    %% Theta and W togather

    data_theta_W = csvread('data/common_pred.csv'); 
    d_theta=reconstructMatrix(N,1, data_theta_W(zz,17:end));
    d_W = reconstructMatrix(M, K, data_theta_W(zz,1:16));
    snr_togather = optimizeBeamforming(d_theta, d_W, M, K, N, Hd, Hr, G, omega);


    %% ************************************************************************
    %% W and using that W Theta is predicted

    data_theta = csvread('data/theta_after_w_pred.csv'); 
    data_W = csvread('data/only_w_pred.csv'); 
    d_theta=reconstructMatrix(N, 1, data_theta(zz,:));
    d_W = reconstructMatrix(M, K, data_W(zz,:));
    snr_hybrid = optimizeBeamforming(d_theta, d_W, M, K, N, Hd, Hr, G, omega);


    %% ************************************************************************
    %% non ML

    theta=exp(1j.*rand(N,1).*2.*pi);
    Theta=diag(theta);
    H=Hd+Hr*Theta*G;
    W_old=zeros(M,K);
    beta=zeros(1,K);
    [ W,grt,f0 ] = init_W( H,M,K,Pt,omega);
    f1=0;
    W_span=W;
    W_last=W;
    [ ~,L_last ] = Proxlinear_beam_para( H,K,M,beta );
    flag=0;
    t_old=1;
    [ beta ] =upadte_beta( H,W,K,grt);

    snr_nonml = [f0];
    for con0=1:15
        [ Qx,qx,theta ] = surface_U_v_direct( W,Hd,Hr,Theta,G,N,K,grt,beta );
        theta_old=theta;

        U=-Qx;v=qx;
        x0=theta_old;
        phi0=angle(x0);
        grad=real((2*U*x0-2*v).*(-1j.*conj(x0)));
        dir=-grad;
        [ Ltheta ] = SCA_phi_step_para( -Qx,qx,N, theta );
        [ theta,t3 ] = armijo_theta( Ltheta,dir,f0,phi0, grad,grt,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);

        [f1,grt,beta,W,W_span,t_old,L_last,H,Theta ] = fun_theta_package(grt,theta,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);

        f0=f1;          
        snr_nonml = [snr_nonml f0];   
    end
    

    % Accumulate SNRs
    snr_gen_all_vector = [snr_gen_all_vector; snr_gen];
    snr_gen_all = snr_gen_all + snr_gen;
    snr_seperate_all_vector = [snr_seperate_all_vector; snr_seperate];
    snr_seperate_all = snr_seperate_all + snr_seperate;
    snr_togather_all_vector = [snr_togather_all_vector; snr_togather];
    snr_togather_all = snr_togather_all + snr_togather;
    snr_hybrid_all_vector = [snr_hybrid_all_vector; snr_hybrid];
    snr_hybrid_all = snr_hybrid_all + snr_hybrid;
    snr_nonml_all_vector = [snr_nonml_all_vector; snr_nonml];
    snr_nonml_all = snr_nonml_all + snr_nonml;

end

% Calculate averages
snr_gen_avg = snr_gen_all / 100;
snr_seperate_avg = snr_seperate_all / 100;
snr_togather_avg = snr_togather_all / 100;
snr_hybrid_avg = snr_hybrid_all / 100;
snr_nonml_avg = snr_nonml_all / 100;

figure(1);

% Plotting with different markers and increased line width
plot(snr_gen_avg, 'LineWidth', 1.5, 'Marker', 'o'); % Circle markers
hold on;
plot(snr_seperate_avg, 'LineWidth', 1.5, 'Marker', 'x'); % Cross markers
hold on;
plot(snr_togather_avg, 'LineWidth', 1.5, 'Marker', '^'); % Triangle markers
hold on;
plot(snr_hybrid_avg, 'LineWidth', 1.5, 'Marker', 's'); % Square markers
hold on;
plot(snr_nonml_avg, 'LineWidth', 1.5, 'Marker', 'd'); % Diamond markers
hold off;

% Setting grid and axis properties
grid on;
set(gca, 'GridAlpha', 0.3); % Adjust grid opacity (0 = invisible, 1 = solid)
set(gca, 'FontSize', 12, 'FontWeight', 'bold'); % Default font size for axes (will be overridden for tick labels)

% Increase the font size of only the numbers (tick labels) on the axes
ax = gca; % Get current axes
ax.XAxis.FontSize = 25; % Increase x-axis tick label font size
ax.YAxis.FontSize = 25; % Increase y-axis tick label font size

% Setting font sizes for labels and title
xlabel('Iteration Number', 'FontSize', 20, 'FontWeight', 'bold'); % Label font size and weight
ylabel('Average WSR', 'FontSize', 20, 'FontWeight', 'bold'); % Label font size and weight
title('Average WSR vs Iteration', 'FontSize', 16, 'FontWeight', 'bold'); % Title font size and weight

% Setting legend in the bottom-right corner with larger font size
legend('Base Line', 'Approach 1(Joint Prediciton)', 'Approach 2(Seperate Predicton)', 'Apporach 3(Beaforming prediction followed by phase predction)', 'Radom Initilization', 'FontSize', 14, 'Location', 'southeast', 'FontWeight', 'bold'); 
