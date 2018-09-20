%% Noise in Accum 
rho_i = 917;
rho_w = 1028;

g = 9.81;
n = 3;
m = 1/3;

year = 3600*24*365;
accum = 0.5/year;

A_glen =4.227e-25;
C = 7.624e6;
theta0 = 0.9;

omega0 = ((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
beta = (m+n+3)/(m+1);
lambda = rho_w/rho_i;

parameters.icedivide = -100;
parameters.bedslope = -1e-3;
parameters.sill_min = 3000e3;
parameters.sill_max = 3010e3;
parameters.sill_slope = 1e-4;

%NL Model
tf = 10e3.*year;
nt = 10e3;
dt = tf/nt;

h=2000;
xg = 600e3; %initial condition

%% Run two-stage model
for t = 1:nt %100kyr to steady state
    b = Base(xg,parameters);
    bx = dBasedx(xg,parameters);
    omega = omega0;%((A_glen*(rho_i*g)^(n+1) * (1-(rho_i/rho_w))^n / (4^n * C))^(1/(m+1))) * theta0^(n/(m+1));
    
    hg = -(rho_w/rho_i)*b;
    Q = (rho_i*g/(C*xg))^n * (h^(2*n + 1));
   
    Q_g = omega*(hg^beta);
    
    dh_dt = accum - (Q_g/xg) - (h/(xg*hg))*(Q-Q_g);
    dxg_dt = (Q-Q_g)/hg;

    h = h + dh_dt*dt;
    xg = xg + dxg_dt*dt;
    
    xgs_nl(t) = xg;
    hs_nl(t) = h;
end