clear all;
% close all;clc

%% Parameters
%%Set Grid Resolution
parameters.grid.n_nodes = 1000;      %Horizontal Resolution
parameters.grid.n2_nodes = 40;
parameters.grid.gz_nodes = 203;      %Horizontal Resolution in grounding zone
parameters.grid.sigma_gz = 0.97;
parameters.Dupont_G = 0;          %lateral shear stress

parameters.year = 3600*24*365;                     %length of a year in seconds
parameters.tfinal = 1000e3.*parameters.year;          %total time of integration
parameters.nsteps = 1000e3;                           %number of time steps

parameters.accumrate = 0.3./parameters.year;
% parameters.accum_new = 0.48/parameters.year;
parameters.accum_mean = parameters.accumrate;
parameters.accum_std = 0.1./parameters.year;

parameters.buttress = 0.4;
parameters.buttress_mean = parameters.buttress;
parameters.buttress_std = 0;

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.C_mean   = parameters.C_schoof;
parameters.C_std    = 0*parameters.C_mean;

parameters.C_noise_list = parameters.C_std.*randn(parameters.nsteps,1);
parameters.buttress_noise_list = parameters.buttress_std.*randn(parameters.nsteps,1);
parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);

parameters.icedivide = -100;
parameters.bedslope = -1e-3;
parameters.sill_min = 3000e3;
parameters.sill_max = 3010e3;
parameters.sill_slope = 1e-4;

%%Time step parameters
parameters.dtau = parameters.tfinal/parameters.nsteps; %length of time steps
parameters.dtau_max = parameters.dtau;

%%Newton Parameters
parameters.HS_sensitivity = pi*parameters.year;     %sensitivity of the HS function (as this gets larger, theta approaches the actual HS function)
parameters.uverbose = 1;
parameters.iteration_threshold = 1e-3;
parameters.hiter_max=1e3;
parameters.uiter_max=5e2;
parameters.titer_max=4e1;
parameters.CFL=50;

%%Grid Parameters
parameters.grid.n_elements = parameters.grid.n_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
% parameters.grid.sigma_node = linspace(0,1,parameters.grid.n_nodes)';  %node positions scaled to (0,1)
% parameters.grid.sigma_node = flipud(1-linspace(0,1,parameters.grid.n_nodes)'.^parameters.grid.n_exponent); %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_node = [linspace(0,0.97,parameters.grid.n_nodes-parameters.grid.gz_nodes),linspace(0.97+(.03/parameters.grid.gz_nodes),1,parameters.grid.gz_nodes)]'; %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_element =...
    (parameters.grid.sigma_node(1:parameters.grid.n_nodes-1)+...
    parameters.grid.sigma_node(2:parameters.grid.n_nodes))/2;     %element centres scaled to (0,1)

parameters.grid.n2_elements = parameters.grid.n2_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
parameters.grid.eta_node = linspace(0,1,parameters.grid.n2_nodes)';  %eta node positions scaled to (0,1)
parameters.grid.eta_element =...
    (parameters.grid.eta_node(1:parameters.grid.n2_nodes-1)+...
    parameters.grid.eta_node(2:parameters.grid.n2_nodes))/2;     %eta element centres scaled to (0,1)

%%Glen's Law parameters
parameters.B_Glen = (4.227e-25^(-1/3)).* ones(parameters.grid.n_elements,1);                     %B in Glen's law (vertically averaged if necessary)
parameters.n_Glen = 3;

%%Physical parameters
parameters.rho = 917;  %917                                 %ice density
parameters.rho_w = 1028;  %1028                               %water density
parameters.g = 9.81;                                    %acceleration due to gravity
parameters.D_eps = 1e-10;                               %strain rate regularizer
parameters.u_eps = 1e-9;                %velocity regularizer
parameters.u_in = 0./parameters.year; 

%%Sliding Law Parameters
parameters.frictionlaw = 'Weertman';

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.m_schoof = 1/3;          %See Schoof (2007)

parameters.B_shear = 0;
parameters.width_shear = 1e3;

parameters.float = 1;
n_plots = 1e3;                        %number of plots to make

rho_i = parameters.rho;
rho_w = parameters.rho_w;

g = parameters.g;
n = parameters.n_Glen;
m = parameters.m_schoof;

year = parameters.year;
accum = parameters.accumrate;
% accum_new = parameters.accum_new;

A_glen =(parameters.B_Glen(1).^(-3));
C = parameters.C_schoof;
eta = 2.05614;
h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

gz_frac = 0.15;

%% Accum Noise

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=1000;
xg = 250e3; %initial guess

%run to steady state first
for t = 1:2e5 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    %     Q = (rho_i*g/(C*xg))^n * (h^(2*n + 1));
    Q = (rho_i*g/C)^n * ((h-hg)/(gz_frac*xg))^n * (h^(n + 1));
%     Q = (A_glen/(n+2))*((rho_i*g/xg)^n)*(h^(2*n+2));
%     Q = (rho_i*g/(C*xg*(1-sigma)))^n * (h^(2*n + 1));
%     Q = (rho_i*g/(C*xg))^n * (h^(n + 1)) * (h-hg)^n;
    Q_g = omega*(hg^beta);
    
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
%     dxg_dt = (Q-Q_g)/hg;
%     dh_dt = accum - (Q/(xg*(1-sigma))) - (h/(xg*sigma*(hg - lambda*bx*xg)))*(accum*sigma*xg + Q - Q_g);
%     dxg_dt = (accum*sigma*xg + Q - Q_g)/(sigma*(hg - lambda*bx*xg));
    dh_dt = accum - (Q/xg) - (h/(xg*(hg - lambda*bx*xg)))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/(hg - lambda*bx*xg);
    
    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end
%%
xg0 = xg;
accum0 = accum;
omega0 = omega;

%Change something from steady state

xgs_nl = xg;

for t = 1:nt
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    
    accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
%     theta_tot = (1 - parameters.buttress_mean)^(n/(m+1)) - parameters.buttress_noise_list(t)/sqrt(dt/year);
%     theta = 1-(parameters.buttress+(parameters.buttress_noise_list(t)/sqrt(dt/year)));
%     C = max([parameters.C_mean+(parameters.C_noise_list(t)/sqrt(dt/year)),1e5]);
    theta = theta0;

    omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta^(n/(m+1));
%     omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta_tot;

    hg = -(rho_w/rho_i)*b;
%     Q = (rho_i*g/(C*xg))^n * (h^(2*n + 1));
    Q = (rho_i*g/C)^n * ((h-hg)/(gz_frac*xg))^n * (h^(n + 1));
%     Q = (A_glen/(n+2))*((rho_i*g/xg)^n)*(h^(2*n+2));
%     Q = (rho_i*g/(C*xg*(1-sigma)))^n * (h^(2*n + 1));
%     Q = (rho_i*g/(C*xg))^n * (h^(n + 1)) * (h-hg)^n;
    Q_g = omega*(hg^beta);
    
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
%     dxg_dt = (Q-Q_g)/hg;
%     dh_dt = accum - (Q/(xg*(1-sigma))) - (h/(xg*sigma*(hg - lambda*bx*xg)))*(accum*sigma*xg + Q - Q_g);
%     dxg_dt = (accum*sigma*xg + Q - Q_g)/(sigma*(hg - lambda*bx*xg));
    dh_dt = accum - (Q/xg) - (h/(xg*(hg - lambda*bx*xg)))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/(hg - lambda*bx*xg);
    
    h  = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    xgs_nl(t) = xg;
    
    Qs(t) = Q;
    Qgs(t) = Q_g;
    
    dxgs(t) = dxg_dt;
    
end

% figure(3);plot(linspace(0,tf,length(xgs_nl))/year./1e3,Qs,'k','markersize',20);hold on
% figure(3);plot(linspace(0,tf,length(xgs_nl))/year./1e3,Qgs,'r','markersize',20);hold on
% figure(2);plot(linspace(0,tf,length(xgs_nl))/year./1e3,dxgs,'r','markersize',20);hold on

%Linearized Simple Model (1)
% tau_n = 5;
% 
% bfun0 = Base(xg0,parameters);
% bx0 = dBasedx(xg0,parameters);
% 
% % hbar = -lambda*(7/11)*(bx0*xg0 + bfun0) + (11/7)*((4*7*eta*C*(accum0^(1/n)) / (11*4*rho_i*g))^(3/7)) * (xg0^(4/7));
% hbar = tau_n*h0 - lambda*(bfun0 + xg0*bx0);
% tau =-(hbar/(tau_n+1))*(accum0 + beta*omega0*(rho_w/rho_i) * ((-(rho_w/rho_i)*bfun0)^(beta-1)) * bx0)^(-1);
% 
% P_coeff = (tau_n+1)*(1/hbar)*xg0;
% omega_coeff = -(tau_n+1)*(1/hbar)*(-(rho_w/rho_i)*bfun0)^beta;
% 
% figure(1);plot(linspace(0,tf/year)./1e3,(xg0 + (omega_coeff*tau*(omega-omega0)).*(1-exp(-linspace(0,tf)/tau)))./1e3,'k','linewidth',2);hold on;
figure(1);
% subplot(2,1,1)
plot(linspace(0,tf,length(xgs_nl))/year./1e3,xgs_nl/1e3,'k','linewidth',2);hold on
xlabel('Time (kyr)','fontsize',20);ylabel('x_g (km)','fontsize',20);set(gca,'fontsize',24);

%% Flowline Model
splot=1;
warning('off')
plotting_on=0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
small_save_on=0;
large_save_on=0;
all_save_on = 1;

%%Pre-allocate storage
h_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
xg_all = nan*ones(1,round(parameters.nsteps/5));
time_all = nan*ones(1,round(parameters.nsteps/5));
% ub_all = nan*ones(parameters.grid.n_nodes,round(parameters.nsteps/5));

%%Initialize variables
disp('Initializing Variables...')
load ConvergedWeertman_linearbed_ac0.3_but40_T-12_N1000_208.mat   %load converged solution from earlier run

h = interp1([0;sigma_elems;1],[2*h(1)-h(2);h;-(parameters.rho_w/parameters.rho).*Base(x_g,parameters)],parameters.grid.sigma_element);
% x_g = xg0;
% h_g = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
% h = h - (h-h_g).*parameters.grid.sigma_element.^3;
u = 100.*parameters.grid.sigma_node./parameters.year;

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_fix(u,h,x_g,parameters); %initialize velocity
% parameters.hx_g_old = [h;x_g];
% v = mass_continuity(h,x_g,u_old_full,parameters);

%%Run integration
Flowline_driver_noisy

