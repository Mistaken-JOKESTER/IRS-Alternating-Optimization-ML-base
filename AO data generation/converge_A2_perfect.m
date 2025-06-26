close all
clear all

run('generate_location.m');
run('generate_pathloss.m');
run('generate_channel.m');

%% system parameters
load('user_channel.mat','K','N','M','ite','pd','ps','Hd_w','theta_init','AP_angle','IRS_angle',...
    'G_sig','User_angle','Hr_sig','eb1','eb2','path_d','path_i', 'us');

etai=1;
weight=1./((path_d));
weight=weight./sum(weight);
omega=weight;
%%
snr=0;
Pt=10.^(snr/10); 
it=1e2;
rate_w=zeros(ite,it);
start=zeros(1,ite);
error_Ellipsoid=0;
beamforming_error=0;
armijo_w=zeros(ite,it);
snrzzz = [];

for t0=1:ite
    
    rate_tmp=zeros(1,it);
    W_old=zeros(M,K);
    beta=zeros(1,K);
    %%
    Hd=pd.*Hd_w(:,:,t0);
    theta=theta_init(:,:,t0);
    Theta=diag(theta');
    %%
    G=channel_G(AP_angle,IRS_angle,G_sig(:,:,t0),eb1,eb2,N,M);
    Hr=ps.*channel_Hr(User_angle,Hr_sig(:,:,t0),eb1,eb2,K,N);
    %%
    appendMatrixToCSV("hr.csv", Hr);
    appendMatrixToCSV("G.csv", G);
    appendMatrixToCSV("Hd.csv", Hd);
    appendMatrixToCSV("us.csv", us);
    appendMatrixToCSV("omega.csv", omega);
    H=Hd+Hr*Theta*G;
    %% init
    [ W,grt,f0 ] = init_W( H,M,K,Pt,omega);
    w_intial = W;
     start(t0)=f0;
     f1=0;
     %%
     W_span=W;
     W_last=W;
     [ ~,L_last ] = Proxlinear_beam_para( H,K,M,beta );
     %%
     flag=0;
     t_old=1;
     [ beta ] =upadte_beta( H,W,K,grt);
     snr = [f0];
     for con0=1:it
        [ Qx,qx,theta ] = surface_U_v_direct( W,Hd,Hr,Theta,G,N,K,grt,beta );
        theta_old=theta;
        %%
        U=-Qx;v=qx;
        x0=theta_old;
        phi0=angle(x0);
        grad=real((2*U*x0-2*v).*(-1j.*conj(x0)));
        dir=-grad;
        [ Ltheta ] = SCA_phi_step_para( -Qx,qx,N, theta );
        [ theta,t3,armijo_w(t0,con0) ] = armijo_theta( Ltheta,dir,f0,phi0, grad,grt,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);
        %%
        [f1,grt,beta,W,W_span,t_old,L_last,H,Theta ] = fun_theta_package(grt,theta,W,W_span,t_old,L_last,K,M,Pt,omega,Hd,Hr,G);
         %%
         rate_tmp(con0)=f1;
         f0=f1;     
         snr = [snr f0];
         
     end

     appendMatrixToCSV("W.csv", W);
     appendMatrixToCSV("theta.csv", diag(Theta));
     appendMatrixToCSV("phi.csv", angle(diag(Theta)));
     rate_w(t0,:)=rate_tmp;
end
 
 