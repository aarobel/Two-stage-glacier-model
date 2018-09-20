clear all;
close all;clc

%% Parameters
%%Set Grid Resolution
parameters.grid.n_nodes = 1000;      %Horizontal Resolution
parameters.grid.n2_nodes = 40;
parameters.grid.gz_nodes = 203;      %Horizontal Resolution in grounding zone
parameters.grid.sigma_gz = 0.97;
parameters.Dupont_G = 0;          %lateral shear stress

parameters.year = 3600*24*365;                     %length of a year in seconds
parameters.tfinal = 60e3.*parameters.year;          %total time of integration
parameters.nsteps = 60e3;                           %number of time steps

parameters.accumrate = 1./parameters.year;
% parameters.accum_new = 0.48/parameters.year;
parameters.accum_mean = parameters.accumrate;
% parameters.accum_std = 0.1./parameters.year;

parameters.buttress = 0.5;
parameters.buttress_mean = parameters.buttress;
parameters.buttress_std = 0;

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.C_mean   = parameters.C_schoof;
% parameters.C_std    = 0.2*parameters.C_mean;

rdn_list = randn(parameters.nsteps,1);

parameters.C_noise_list10 = 0.1*parameters.C_mean.*rdn_list;
parameters.C_noise_list20 = 0.2*parameters.C_mean.*rdn_list;
parameters.buttress_noise_list1 = 0.10*parameters.buttress_mean.*rdn_list;
parameters.buttress_noise_list2 = 0.20*parameters.buttress_mean.*rdn_list;
parameters.accum_noise_list1 =  0.10*parameters.accum_mean.*rdn_list;
parameters.accum_noise_list2 =  0.20*parameters.accum_mean.*rdn_list;

parameters.icedivide = -100;
parameters.bedslope = -1e-3;
parameters.sill_min = 3000e3;
parameters.sill_max = 3010e3;
parameters.sill_slope = 1e-4;
parameters.sin_amp = 0;
parameters.sin_length = 10e3;

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

gz_frac = 1;

%% Noise in Accum
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
% h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega0 = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

gz_frac = 1;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=2000;
xg = 500e3; %initial guess

%run to steady state first
for t = 1:1e5 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end

xg0 = xg;
h0 = h;
hg0 = hg;
accum0 = accum;
omega0 = omega;
Q_g0 = Q_g;

%add noise

xgs_nl = xg;
C_gz = C;

for t = 1:nt
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    
    theta = theta0;
    if(t>nt/3 & t<2*nt/3)
        accum = parameters.accum_mean+(parameters.accum_noise_list1(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list1(t)/sqrt(dt/year));
%             omega = omega0+(parameters.omega_noise_list1(t)/sqrt(dt/year));
    else if(t>2*nt/3)
        accum = parameters.accum_mean+(parameters.accum_noise_list1(t)/sqrt(dt/year));             
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list2(t)/sqrt(dt/year));
%         omega = omega0+(parameters.omega_noise_list2(t)/sqrt(dt/year));
        end;end
    theta0 = 1-parameters.buttress;
    omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C_gz))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));

    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h  = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    xgs_nl(t) = xg;
    
    Qs(t) = Q;
    Qgs(t) = Q_g;
    
    dxgs(t) = dxg_dt;
    
end
figure(1);set(1,'units','normalized','position',[0.3 0.1 0.5 0.65]);
subplot(2,2,1)
plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',2);hold on
plot(linspace(0,nt/3,nt/3)./1e3,smooth((xgs_nl(1:nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,smooth((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,smooth((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(0,nt/3,nt/3)./1e3,mean((xgs_nl(1:nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,mean((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,mean((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
xlabel('Time (kyr)','fontsize',24)
ylabel('Grounding Line Deviation (km)','fontsize',20)
set(gca,'fontsize',24)
ylim([-2 2])
title('Noisy Surface Mass Balance','fontsize',24)
text(0.02,1.0,'a','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)

% %% Noise in GL melt
% rho_i = parameters.rho;
% rho_w = parameters.rho_w;
% 
% g = parameters.g;
% n = parameters.n_Glen;
% m = parameters.m_schoof;
% 
% year = parameters.year;
% accum = parameters.accumrate;
% A_glen =(parameters.B_Glen(1).^(-3));
% C = parameters.C_schoof;
% eta = 2.05614;
% 
% theta0 = 1-parameters.buttress;
% omega0 = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
% beta = (m+n+3)/(m+1);
% lambda = rho_w/rho_i;
% 
% gz_frac = 1;
% 
% M0 = 200/year;
% 
% parameters.M_noise_list1 = 0.20*M0.*rdn_list;
% parameters.M_noise_list2 = 0.40*M0.*rdn_list;
% 
% %NL Model
% tf = parameters.tfinal;
% nt = parameters.nsteps;
% dt = tf/nt;
% 
% h=2000;
% xg = 500e3; %initial guess
% 
% %run to steady state first
% for t = 1:1e5 %100kyr to steady state
%     b = Base(xg,parameters);
%     bx = dBasedx(xg,parameters);
%     omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
%     M = M0;
%     
%     hg = -(rho_w/rho_i)*b;
%     Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
%    
%     Q_g = omega*(hg^beta);
%     
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g-(hg*M));
%     dxg_dt = (Q-Q_g)/hg - M;
% 
%     h = h + dh_dt*4*year;
%     xg = xg + dxg_dt*4*year;
%     xgs_nl(t) = xg;
%     
% end
% 
% xg0 = xg;
% h0 = h;
% hg0 = hg;
% accum0 = accum;
% omega0 = omega;
% 
% %add noise
% 
% xgs_nl = xg;
% C_gz = C;
% 
% for t = 1:nt
%     b = Base(xg,parameters);
%     bx = dBasedx(xg,parameters);
%     
%     theta = theta0;
%     if(t>nt/3 & t<2*nt/3)
% %         accum = parameters.accum_mean+(parameters.accum_noise_list1(t)/sqrt(dt/year));
% %         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list1(t)/sqrt(dt/year));
% %             omega = omega0+(parameters.omega_noise_list1(t)/sqrt(dt/year));
%         M = M0+(parameters.M_noise_list1(t)/sqrt(dt/year));
%     else if(t>2*nt/3)
% %         accum = parameters.accum_mean+(parameters.accum_noise_list1(t)/sqrt(dt/year));             
% %         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list2(t)/sqrt(dt/year));
% %         omega = omega0+(parameters.omega_noise_list2(t)/sqrt(dt/year));
%         M = M0+(parameters.M_noise_list2(t)/sqrt(dt/year));
%         end;end
%     theta0 = 1-parameters.buttress;
%     omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C_gz))^(1/(m+1))) * theta0^(n/(m+1));
%     
%     hg = -(rho_w/rho_i)*b;
%     Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
% 
%     Q_g = omega*(hg^beta);
%     
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g-(hg*M));
%     dxg_dt = (Q-Q_g)/hg - M;
% 
%     h  = h + dh_dt*dt;
%     xg = xg + dxg_dt*dt;
%     xgs_nl(t) = xg;
%     
%     Qs(t) = Q;
%     Qgs(t) = Q_g;
%     
%     dxgs(t) = dxg_dt;
%     
% end
% figure(1);subplot(2,2,2)
% plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',3);hold on
% plot(linspace(0,nt/3,nt/3)./1e3,smooth((xgs_nl(1:nt/3)-xgs_nl(1))./1e3,1000),'r','linewidth',2);hold on
% plot(linspace(nt/3,2*nt/3,nt/3)./1e3,smooth((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3,1000),'r','linewidth',2);hold on
% plot(linspace(2*nt/3,nt,nt/3)./1e3,smooth((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3,1000),'r','linewidth',2);hold on
% xlabel('Time (kyr)','fontsize',24)
% ylabel('Grounding Line Deviation (km)','fontsize',20)
% set(gca,'fontsize',24)
% ylim([-1 1])
% title('Noisy GL Melt Rate (M'')','fontsize',24)
% text(0.02,1.0,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)

%% Noise in Omega
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
% h0 = 2710;
% h0 = 5000;

theta0 = 1-parameters.buttress;
omega0 = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

parameters.omega_noise_list1 = 0.10*omega0.*rdn_list;
parameters.omega_noise_list2 = 0.20*omega0.*rdn_list;

gz_frac = 1;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=2000;
xg = 500e3; %initial guess

%run to steady state first
for t = 1:1e5 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end

xg0 = xg;
h0 = h;
hg0 = hg;
accum0 = accum;
omega0 = omega;

%add noise

xgs_nl = xg;
C_gz = C;

for t = 1:nt
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    
    theta = theta0;
    if(t>nt/3 & t<2*nt/3)
%         accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
%           C_gz = parameters.C_mean+(parameters.C_noise_list10(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list1(t)/sqrt(dt/year));
            omega = omega0+(parameters.omega_noise_list1(t)/sqrt(dt/year));
    else if(t>2*nt/3)
%             C_gz = parameters.C_mean+(parameters.C_noise_list20(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list2(t)/sqrt(dt/year));
            omega = omega0+(parameters.omega_noise_list2(t)/sqrt(dt/year));
        end;end
%     theta0 = 1-parameters.buttress;
%     omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C_gz))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));

    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h  = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    xgs_nl(t) = xg;
    
    Qs(t) = Q;
    Qgs(t) = Q_g;
    
    dxgs(t) = dxg_dt;
    
end
figure(1);figure(1);subplot(2,2,2)
plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',2);hold on
plot(linspace(0,nt/3,nt/3)./1e3,smooth((xgs_nl(1:nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,smooth((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,smooth((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(0,nt/3,nt/3)./1e3,mean((xgs_nl(1:nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,mean((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,mean((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
xlabel('Time (kyr)','fontsize',24)
ylabel('Grounding Line Deviation (km)','fontsize',20)
set(gca,'fontsize',24)
ylim([-2 2])
title('Noisy GL Flux Coefficient','fontsize',24)
text(0.02,1.0,'b','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)

% %% Noise in Buttressing
% rho_i = parameters.rho;
% rho_w = parameters.rho_w;
% 
% g = parameters.g;
% n = parameters.n_Glen;
% m = parameters.m_schoof;
% 
% year = parameters.year;
% accum = parameters.accumrate;
% % accum_new = parameters.accum_new;
% 
% A_glen =(parameters.B_Glen(1).^(-3));
% C = parameters.C_schoof;
% eta = 2.05614;
% % h0 = 2710;
% % h0 = 5000;
% 
% theta0 = 1-parameters.buttress;
% omega0 = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
% beta = (m+n+3)/(m+1);
% lambda = rho_w/rho_i;
% 
% 
% gz_frac = 1;
% 
% %NL Model
% tf = parameters.tfinal;
% nt = parameters.nsteps;
% dt = tf/nt;
% 
% h=2000;
% xg = 500e3; %initial guess
% 
% %run to steady state first
% for t = 1:1e5 %100kyr to steady state
%     b = Base(xg,parameters);
%     bx = dBasedx(xg,parameters);
%     omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
%     
%     hg = -(rho_w/rho_i)*b;
%     Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
%    
%     Q_g = omega*(hg^beta);
%     
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
%     dxg_dt = (Q-Q_g)/hg;
% 
%     h = h + dh_dt*4*year;
%     xg = xg + dxg_dt*4*year;
%     xgs_nl(t) = xg;
%     
% end
% 
% xg0 = xg;
% h0 = h;
% hg0 = hg;
% accum0 = accum;
% omega0 = omega;
% 
% %add noise
% 
% xgs_nl = xg;
% C_gz = C;
% 
% for t = 1:nt
%     b = Base(xg,parameters);
%     bx = dBasedx(xg,parameters);
%     
%     theta = theta0;
%     if(t>nt/3 & t<2*nt/3)
% %         accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
% %           C_gz = parameters.C_mean+(parameters.C_noise_list10(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list1(t)/sqrt(dt/year));
% %             omega = omega0+(parameters.omega_noise_list1(t)/sqrt(dt/year));
%     else if(t>2*nt/3)
% %             C_gz = parameters.C_mean+(parameters.C_noise_list20(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list2(t)/sqrt(dt/year));
% %             omega = omega0+(parameters.omega_noise_list2(t)/sqrt(dt/year));
%         end;end
%     theta0 = 1-parameters.buttress;
%     omega = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C_gz))^(1/(m+1))) * theta0^(n/(m+1));
%     
%     hg = -(rho_w/rho_i)*b;
%     Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
% 
%     Q_g = omega*(hg^beta);
%     
%     dh_dt = accum - (Q/xg) - (h/(xg*hg))*(Q-Q_g);
%     dxg_dt = (Q-Q_g)/hg;
% 
%     h  = h + dh_dt*dt;
%     xg = xg + dxg_dt*dt;
%     xgs_nl(t) = xg;
%     
%     Qs(t) = Q;
%     Qgs(t) = Q_g;
%     
%     dxgs(t) = dxg_dt;
%     
% end
% figure(1);figure(1);subplot(2,2,4)
% plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',3);hold on
% plot(linspace(0,nt/3,nt/3)./1e3,mean((xgs_nl(nt/6+1:nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'b','linewidth',3);hold on
% plot(linspace(nt/3,2*nt/3,nt/3)./1e3,mean((xgs_nl(nt/2+1:2*nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'b','linewidth',3);hold on
% plot(linspace(2*nt/3,nt,nt/3)./1e3,mean((xgs_nl(5*nt/6+1:nt)-xgs_nl(1))./1e3).*ones(1,nt/3),'b','linewidth',3);hold on
% xlabel('Time (kyr)','fontsize',24)
% ylabel('Grounding Line Deviation (km)','fontsize',24)
% set(gca,'fontsize',24)
% ylim([-50 5])
% title('Noisy Buttessing Param (\theta'')','fontsize',24)

%% Noise in Ice Shelf Length
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
L_s = 40e3;
W = 20e3;

theta0 = 1-parameters.buttress;
beta = n+1;
lambda = rho_w/rho_i;
omega0 = (n/2)^n * (n+1)^(-(n+1)) * (rho_i*g*(1-(rho_i/rho_w)))^n * A_glen * L_s^(-n) * W^(n+1);

L_s0 = L_s;
parameters.Ls_noise_list1 = 0.10*L_s0.*rdn_list;
parameters.Ls_noise_list2 = 0.20*L_s0.*rdn_list;


gz_frac = 1;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=2000;
xg = 500e3; %initial guess

%run to steady state first
for t = 1:1e5 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*4*year;
    xg = xg + dxg_dt*4*year;
    xgs_nl(t) = xg;
    
end

xg0 = xg;
h0 = h;
hg0 = hg;
accum0 = accum;
omega0 = omega;
Q0 = Q;
%add noise

xgs_nl = xg;
C_gz = C;

for t = 1:nt
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    
    theta = theta0;
    if(t>nt/3 & t<2*nt/3)
%         accum = parameters.accum_mean+(parameters.accum_noise_list(t)/sqrt(dt/year));
%           C_gz = parameters.C_mean+(parameters.C_noise_list10(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list1(t)/sqrt(dt/year));
%             omega = omega0+(parameters.omega_noise_list1(t)/sqrt(dt/year));
          L_s = L_s0+(parameters.Ls_noise_list1(t)/sqrt(dt/year));
    else if(t>2*nt/3)
%             C_gz = parameters.C_mean+(parameters.C_noise_list20(t)/sqrt(dt/year));
%         parameters.buttress = parameters.buttress_mean+(parameters.buttress_noise_list2(t)/sqrt(dt/year));
%             omega = omega0+(parameters.omega_noise_list2(t)/sqrt(dt/year));
        L_s = L_s0+(parameters.Ls_noise_list2(t)/sqrt(dt/year));
        end;end
%     theta0 = 1-parameters.buttress;
    omega = (n/2)^n * (n+1)^(-(n+1)) * (rho_i*g*(1-(rho_i/rho_w)))^n * A_glen * L_s^(-n) * W^(n+1);
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));

    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h  = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    xgs_nl(t) = xg;
    
    Qs(t) = Q;
    Qgs(t) = Q_g;
    
    dxgs(t) = dxg_dt;
    
end
figure(1);figure(1);subplot(2,2,3)
plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',2);hold on
plot(linspace(0,nt/3,nt/3)./1e3,smooth((xgs_nl(1:nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,smooth((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,smooth((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(0,nt/3,nt/3)./1e3,mean((xgs_nl(1:nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,mean((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,mean((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
xlabel('Time (kyr)','fontsize',24)
ylabel('Grounding Line Deviation (km)','fontsize',20)
set(gca,'fontsize',24)
ylim([-140 30]);xlim([0 nt/1e3])
title('Noisy Ice Shelf Length','fontsize',24)
text(0.02,1.0,'c','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)

%% Noise in Ice Shelf Basal Melt
accum = 1/year;
A_glen =(parameters.B_Glen(1).^(-3));
C = parameters.C_schoof;
basal_melt = -10/year;
L_s = 50e3;
W = 5e3;

% omega0 = (n/2)^n * (n+1)^(-(n+1)) * (rho_i*g*(1-(rho_i/rho_w)))^n * A_glen * L_s^(-n) * W^(n+1);

basal_melt0 = basal_melt;
parameters.basal_melt_noise_list1 = 0.10*basal_melt0.*rdn_list;
parameters.basal_melt_noise_list2 = 0.20*basal_melt0.*rdn_list;

%NL Model
tf = parameters.tfinal;
nt = parameters.nsteps;
dt = tf/nt;

h=h0;
xg = xg0;
hg = hg0;
xgs_nl = xg;

omega0 = n^2 * (n+1)^(-1/n) * (rho_i*g*(1-(rho_i/rho_w))/2) * A_glen^(1/n) * W^((n+1)/n);
qg_last = fsolve(@(q) qg_bm(q,basal_melt,omega0,L_s,hg,n) , Q_g0 ,optimoptions(@fsolve,'algorithm','levenberg-marquardt'));
for t = 1:1e4 %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
%     omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));
        
%     Q_g = omega*(hg^beta);
    Q_g = fsolve(@(q) qg_bm(q,basal_melt,omega0,L_s,hg,n) , qg_last,optimoptions(@fsolve,'algorithm','levenberg-marquardt'));
    qg_last = Q_g;
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*10*year;
    xg = xg + dxg_dt*10*year;
    xgs_nl(t) = xg;
    hs(t) = h;
    
%     if(mod(t,10)==0)
%         figure(2);plot(t,xg/1e3,'k.','markersize',20);hold on;drawnow
%         t
%     end
end

xg0 = xg;
h0 = h;
hg0 = hg;
qg_last = Q_g;

%add noise
xgs_nl = xg;
C_gz = C;

for t = 1:nt
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    
    basal_melt = basal_melt0;
    if(t>nt/3 && t<2*nt/3)
      basal_melt = basal_melt0+(parameters.basal_melt_noise_list1(t)/sqrt(dt/year));
    else if(t>2*nt/3)
      basal_melt = basal_melt0+(parameters.basal_melt_noise_list2(t)/sqrt(dt/year));
    end;end
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*gz_frac*xg))^n * (h^(2*n + 1));

%     Q_g = omega*(hg^beta);
    Q_g = fsolve(@(q) qg_bm(q,basal_melt,omega0,L_s,hg,n) , qg_last,optimoptions(@fsolve,'algorithm','levenberg-marquardt','Display', 'off'));
    qg_last = Q_g;
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h  = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    xgs_nl(t) = xg;
    
    Qs(t) = Q;
    Qgs(t) = Q_g;
    
    dxgs(t) = dxg_dt;
    
    
    if(mod(t,100)==0)
%         figure(2);plot(1:t,xgs_nl/1e3,'k-','markersize',20);hold on;drawnow
        t
    end
    
end
figure(1);figure(1);subplot(2,2,4)
plot(linspace(0,nt,nt)./1e3,(xgs_nl-xgs_nl(1))./1e3,'k','linewidth',2);hold on
plot(linspace(0,nt/3,nt/3)./1e3,smooth((xgs_nl(1:nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,smooth((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,smooth((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3,1000),'b','linewidth',4);hold on
plot(linspace(0,nt/3,nt/3)./1e3,mean((xgs_nl(1:nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(nt/3,2*nt/3,nt/3)./1e3,mean((xgs_nl(nt/3+1:2*nt/3)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
plot(linspace(2*nt/3,nt,nt/3)./1e3,mean((xgs_nl(2*nt/3+1:nt)-xgs_nl(1))./1e3).*ones(1,nt/3),'r','linewidth',4);hold on
xlabel('Time (kyr)','fontsize',24)
ylabel('Grounding Line Deviation (km)','fontsize',20)
set(gca,'fontsize',24)
ylim([-2 2])
title('Noisy Basal Melt','fontsize',24)
text(0.02,1.0,'d','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',40)