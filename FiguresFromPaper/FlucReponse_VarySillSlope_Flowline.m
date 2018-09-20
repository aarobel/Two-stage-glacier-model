clear all;
close all;clc
% load Verif_Flowline_AccumNoise_LinearBed_StrongButtress.mat

accumrate0 = 0.3./(3600*24*365);
accumrate1 = 0.29./(3600*24*365);
%% Flowline Model Run - Step
%%Set Grid Resolution
parameters.grid.n_nodes = 1000;      %Horizontal Resolution
parameters.grid.n2_nodes = 40;
parameters.grid.gz_nodes = 188;      %Horizontal Resolution in grounding zone
parameters.grid.sigma_gz = 0.97;
parameters.Dupont_G = 0;          %lateral shear stress

parameters.year = 3600*24*365;                     %length of a year in seconds
parameters.tfinal = 40e3.*parameters.year;          %total time of integration
parameters.nsteps = 40e3;                           %number of time steps

parameters.accumrate = accumrate0;
% parameters.accum_new = 0.48/parameters.year;
parameters.accum_mean = parameters.accumrate;
parameters.accum_std = 0./parameters.year;

parameters.buttress = 0.4;
parameters.buttress_mean = parameters.buttress;
parameters.buttress_std = 0;

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.C_mean   = parameters.C_schoof;
parameters.C_std    = 0*parameters.C_mean;

parameters.C_noise_list = parameters.C_std.*randn(parameters.nsteps,1);
parameters.buttress_noise_list = parameters.buttress_std.*randn(parameters.nsteps,1);
parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);

parameters.bedslope = -1e-3;
% parameters.icedivide = -100;
parameters.icedivide = -546.375-(parameters.bedslope*4.46375e5);
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

%%Other params
n_plots = 1e2;                        %number of plots to make

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

gz_frac = 1;


%% Run to S-S
splot=0;
warning('off')
plotting_on=1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
small_save_on=0;
large_save_on=0;
all_save_on = 0;

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

%New steady-state
parameters.tfinal = 20e3.*parameters.year;          %total time of integration
parameters.nsteps = 20e3;                           %number of time steps
Flowline_driver_noisy


xg_orig = x_g;
b_orig = Base(xg_orig,parameters);
b0_orig = parameters.icedivide;
bx_orig = parameters.bedslope;
hg_orig = -(parameters.rho_w/parameters.rho).*Base(xg_orig,parameters);
h_orig = h;

tss_stab_slow = -hg_orig/(beta*lambda*xg_orig);
tss_stab_fast = (hg_orig/(beta*lambda*xg_orig))*((3*n+1) + (2*n+1)*(hg_orig/h(1)));

%% Add Sill
close all
x_g = xg_orig;
h = h_orig;

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_fix(u,h,x_g,parameters); %initialize velocity


%add sill
parameters.sill_max = xg_orig-12e3;
parameters.sill_min = xg_orig-17e3;
parameters.sill_slope = bx_orig;%1.5*tss_stab_slow;
parameters.icedivide = b0_orig-(parameters.sill_slope-bx_orig)*(5e3);

%New steady-state
parameters.accum_mean = accumrate0*0.9;
parameters.tfinal = 30e3.*parameters.year;          %total time of integration
parameters.nsteps = 30e3;                           %number of time steps
Flowline_driver_noisy

% figure(2);
% plot(time_all,xg_all./1e3,'b','linewidth',3)

xg_all_fs_ss = xg_all;
%% Add Sill Unstable Slow MISI

x_g = xg_orig;
h = h_orig;

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_fix(u,h,x_g,parameters); %initialize velocity

%add sill
parameters.sill_max = xg_orig-12e3;
parameters.sill_min = xg_orig-17e3;
parameters.sill_slope = 0.5*tss_stab_slow;
parameters.icedivide = b0_orig-(parameters.sill_slope-bx_orig)*(5e3);

%New steady-state
parameters.accum_mean = accumrate0*0.9;
parameters.tfinal = 30e3.*parameters.year;          %total time of integration
parameters.nsteps = 30e3;                           %number of time steps
Flowline_driver_noisy

% figure(2);
% plot(time_all,xg_all./1e3,'r','linewidth',3)

xg_all_fs_su = xg_all;
%% Add Sill Unstable Fast MISI

x_g = xg_orig;
h = h_orig;

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_fix(u,h,x_g,parameters); %initialize velocity

%add sill
parameters.sill_max = xg_orig-12e3;
parameters.sill_min = xg_orig-17e3;
parameters.sill_slope = 1*tss_stab_fast;
parameters.icedivide = b0_orig-(parameters.sill_slope-bx_orig)*(5e3);

%New steady-state
parameters.accum_mean = accumrate0*0.9;
parameters.tfinal = 20e3.*parameters.year;          %total time of integration
parameters.nsteps = 20e3;                           %number of time steps
Flowline_driver_noisy

xg_all_fu_su = xg_all;

%%
figure(2);
plot(time_all,xg_all_fs_ss./1e3,'b','linewidth',3)
plot(time_all,xg_all_fs_su./1e3,'r','linewidth',3)
plot(time_all,xg_all_fu_su./1e3,'k','linewidth',3)