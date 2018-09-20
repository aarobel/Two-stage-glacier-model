clear all;
close all;clc
% load Verif_Flowline_AccumNoise_LinearBed_StrongButtress.mat
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

parameters.accumrate = 0.5./parameters.year;
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
% parameters.accum_noise_list = parameters.accum_std.*randn(parameters.nsteps,1);

parameters.icedivide = 0;
parameters.bedslope = -3e-3;
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
% h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

gz_frac = 1;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=1000;
xg = 500e3; %initial guess

%run to steady state first
for t = 1:2e5 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end
xg_orig = xg;
b_orig = b;
hg_orig = hg;
h_orig = h;

%% Analytic time scales from linearized two-stage mode
accum0=parameters.accumrate;
accums = linspace(0.1*accum0,4*accum0,40);

for j = 1:length(accums)
    
    accum = accums(j);
%     parameters.bedslope = bedslopes(j);
%     parameters.icedivide = b_orig-(parameters.bedslope*xg_orig);
    
    xg = xg_orig;
    h = h_orig;
    for t = 1:2e4 %100kyr to steady state
        b = Base(xg,parameters);
        bx = dBasedx(xg,parameters);
        

        hg = -(rho_w/rho_i)*b;
        Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));

        Q_g = omega*(hg^beta);

        dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
        dxg_dt = (Q-Q_g)/hg;

        h = h + dh_dt*5*year;
        xg = xg + dxg_dt*5*year;
        xgs_nl(t) = xg;

    end
    T_F = ((3*n+1)*Q_g/(hg*xg) + (2*n+1)*Q_g/(h*xg) - beta*lambda*bx*Q_g/(hg^2))^-1;
%     T_S = -(1/Q_g)*(((1/xg)+(beta*lambda*bx/hg))^(-1)) * (((3+(1/n))/(2+(1/n)))*h + hg - ((1/n)/(2+(1/n)))*beta*lambda*bx*xg*h/hg);
    
    T_S = -((1/n)/(2+(1/n)))*(h*hg*(xg^2)/((Q_g^2)*T_F))*(1+beta*lambda*bx*xg/hg)^-1;
    
    
    accums(j) = accum;
    tfs(j) = T_F;
    tss(j) = T_S;
    
    j
end

%% Plot
figure(1);set(1,'units','normalized','position',[0.3 0.1 0.4 0.4]);
% subplot(1,2,1);
semilogy(accums*year,tss/year,'k','linewidth',6);hold on
set(gca,'fontsize',24)
xlabel('Accum (m/yr)','fontsize',24);
ylabel('Time Scale (yrs)','fontsize',24);
% xlim([-3e-3 -1e-4]);


% subplot(1,2,2);
semilogy(accums*year,tfs/year,'r','linewidth',6);hold on
set(gca,'fontsize',24)
% xlabel('Bed Slope (b_x)','fontsize',24);ylabel('Fast Time Scale (yrs)','fontsize',24);
% xlim([-3e-3 -1e-4])
% ylim([10 2e4])
legend('Slow Time Scale','Fast Time Scale')
