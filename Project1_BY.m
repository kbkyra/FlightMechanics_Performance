clc
clear
close all

%% Intitial Comments
%{
AEM-368-001, Spring 2021
Project 1
Kyra Bryan, Kristan Young

Description:
MATLAB program that calculates performance parameters for a UAV
based on the principles discussed in AEM 368
Valid only for the input bounds given in the prompt

Nomenclature and variables:

h_g = flight altitude, geometric
c_r = wing root chord
b = wingspan
taper = wing taper ratio 
sweep = wing quarter-chord sweep 
cl_max = maximum lift coefficient
cd0 = parasitic drag coefficient
W0 = gross weight
c_fuel = fuel capacity
SFC = specific fuel consumption 
PA_SSL = power available at SSL 
eta_p = propeller efficiency

r_earth = radius of earth
g0 = gravitational constant
T_SSL = temperature at SSL
P_SSL = pressure at SSL
rho_SSL = density at SSL
R = gas constant 
a0 = lapse rate for altitudes under 30000ft 
gamma = heat capacity ratio

h = geopotential altitude 
T_alt_R = temperature for a given flight altitude in Rankine
mu = dynamic viscosity
P_alt = pressure for a given flight altitude
rho_alt = density for a given flight altitude
T_alt_F = temperature for a given flight altitude in Fahrenheit 

c_t = tip chord
s = wing planform area
AR = aspect ratio 
e_wing = spanwise efficiency factor
e_oswald = oswald efficiency factor
k = _ constant %FIXME
ld_max = Lift/Drag ratio, max
MAC = mean aerodynamic chord
v_tr_min = min thrust velocity and max L/D velocity 
Re = Reynolds number

W1 = weight with fuel loss
PA = power available (constant)
wing_loading = wing loading
v_stall = stall velocity 
R_max = max range
v_R_max = velocity at max range and max L/D
clthreetwo_cd = max cl^(3/2)/cd
E_max = max endurance
v_E_max = velocity at max endurance and max cl^(3/2)/cd
v_max = max speed 
a_alt = speed of sound at given altitude 
M_max = Mach number based off of max speed and speed of sound at given altitude 
vel_indep = velocity values to iterate to get PR curve;
PR = power required
q = dynamic pressure
cl = coefficient of lift
cd = total coefficient of drag
v_vclimb = max rate of climb (v for vertical)
v_PR_min = horizontal velocity at max rate of climb
ceiling_alt = absolute ceiling at for a given altitude (thus density)
ttc = time to climb to 10000ft 
R_glide_max = max glide range
v_indicate = corresponding velciites for requested power required values
PR_indicate = requested power required values
%}

%% Read Input File, Check If Within Bounds

% Read Input File
data_array = load('specs.txt'); % reads (loads) the data as a column array
h_g = data_array(1); %ft, sets var1 to the value on the first line
c_r = data_array(2); %ft, sets var2 to the value on the second line
b = data_array(3); %ft
taper = data_array(4);
sweep = data_array(5); %deg
cl_max = data_array(6);
cd0 = data_array(7);
W0 = data_array(8); %lb
c_fuel = data_array(9); %lb 
SFC = data_array(10); %lb/(hp-h)
P_shaft = data_array(11); %hp
eta_p = data_array(12);

% Check If Within Bounds, and Modify as Needed 
if h_g < 0
    h_g = 0;
elseif h_g > 30000
    h_g = 30000;
end 
if c_r < 1
    c_r = 1;
elseif c_r > 7
    c_r = 7;
end 
if b < 15
    b = 15;
elseif b > 50
    b = 50;
end 
if taper < 0
    taper = 0;
elseif taper > 1
    taper = 1;
end 
if sweep < 0
    sweep = 0;
elseif sweep > 10
    sweep = 10;
end 
if cl_max < 1
    cl_max = 1;
elseif cl_max > 2
    cl_max = 2;
end 
if cd0 < 0.015
    cd0 = 0.015;
elseif cd0 > 0.04
    cd0 = 0.04;
end 
if W0 < 500
    W0 = 500;
elseif W0 > 5000
    W0 = 5000;
end 
if c_fuel < (.2*W0)
    c_fuel = (.2*W0);
elseif c_fuel > (.4*W0)
    c_fuel = (.4*W0);
end 
if SFC < 0.35
    SFC = 0.35;
elseif SFC > 0.6
    SFC = 0.6;
end 
if P_shaft < 50
    P_shaft = 50;
elseif P_shaft > 500
    P_shaft = 500;
end 
if eta_p < 0.5
    eta_p = 0.5;
elseif eta_p > 0.9
    eta_p = 0.9;
end 

%% Other Inputs and Constants

r_earth = 2.09e7; %ft
g0 = 32.2; %ft/s^2
T_SSL = 518.67; %R
P_SSL = 2116.8; %lb/ft^2
rho_SSL = 0.002377; %slug/ft^3
R = 1716; %ft*lb/(slug*R)
a0 = -0.00357; %R/ft
gamma = 1.4;

%% Standard Atmosphere

h = (h_g*r_earth)/(h_g+r_earth); %ft
T_alt_R = T_SSL + a0*h; %R
T_alt_K = T_alt_R*(5/9); %K
mu_SI = 1.458e-6*((T_alt_K^1.5)/(T_alt_K+110.4)); %kg/(m*s)
mu = (mu_SI*0.0685)/(3.281); %slug/(ft*s)
P_alt = P_SSL*(T_alt_R/T_SSL)^(-g0/(a0*R)); %lb/ft^2
rho_alt = P_alt/(T_alt_R*R); %slug/ft^3
T_alt_F = T_alt_R - 459.67; %F

%% Aerodynamics

c_t = taper*c_r; %ft
s = 0.5*(c_r+c_t)*b; %ft^2
AR = b^2/s; 

%to interpolate e_wing
AR_sample = [20 16 12 8 4 2;
        20 16 12 8 4 2;
        20 16 12 8 4 2;
        20 16 12 8 4 2;
        20 16 12 8 4 2;
        20 16 12 8 4 2];
taper_sample = [1 1 1 1 1 1;
    0.8 0.8 0.8 0.8 0.8 0.8;
    0.6 0.6 0.6 0.6 0.6 0.6;
    0.4 0.4 0.4 0.4 0.4 0.4;
    0.2 0.2 0.2 0.2 0.2 0.2;
    0 0 0 0 0 0];
e_wing_sample = [0.82 0.882 0.907 0.937 0.972 0.990;
                0.891 0.909 0.929 0.952 0.979 0.993;
                0.924 0.938 0.953 0.969 0.987 0.995;
                0.950 0.960 0.970 0.980 0.992 0.997;
                0.937 0.942 0.950 0.962 0.980 0.991;
                0.775 0.783 0.797 0.820 0.865 0.913];
e_wing = interp2(AR_sample,taper_sample,e_wing_sample,AR,taper);

e_oswald = e_wing*0.75;
k = 1/(pi*e_oswald*AR);

ld_max = 0.5*sqrt(1/(cd0*k));

MAC = ((2/3)*c_r)*((taper^2+taper+1)/(taper+1)); %ft
v_tr_min = (k/cd0)^(1/4)*sqrt((2*W0)/(rho_alt*s)); %ft/s
Re = (rho_alt*MAC*v_tr_min)/(mu);

%% Thrust, Power, Performance
%for "50% fuel consumed" parameters 

W1 = W0-(0.5*c_fuel); %lb

PA = P_shaft*eta_p*550; %ft*lb/s

wing_loading = W1/s; %lb/ft^s

v_stall = sqrt((2*W1)/(rho_alt*cl_max*s)); %ft/s

%assume all fuel is burned for max R and max E
W_all = W0-c_fuel; %lb

SFC = SFC/1980000; %1/ft
R_max = (eta_p/SFC)*ld_max*log(W0/W_all); %ft
v_R_max = (k/cd0)^(1/4)*sqrt((2*W1)/(rho_alt*s)); %ft/s

clthreetwo_cd = ((3*(1/k))^(3/4))/(4*(cd0)^(0.25));
E_max = (eta_p/SFC)*clthreetwo_cd*sqrt(2*rho_alt*s)*(W_all^(-1/2)-W0^(-1/2)); %s
v_E_max = (k/(3*cd0))^(1/4)*sqrt((2*W1)/(rho_alt*s)); %ft/s

v_max_inter = roots([0.5*rho_alt*s*cd0 0 0 -PA (2*W1^2*k)/(rho_alt*s)]); %ft/s
v_max = max(real(v_max_inter)); %ft/s
a_alt = sqrt(gamma*R*T_alt_R); %ft/s
M_max = v_max/a_alt;

%velocity dependent values
vel_indep = linspace(v_E_max,v_max);
PR = zeros(1,length(vel_indep));
for i = 1:length(vel_indep) %iterate velocity 
    q = 0.5*rho_alt*vel_indep(i)*vel_indep(i);
    cl = W0/(s*q); 
    cd = cd0 + (cl*cl*k);
    PR(i) = (W0/(cl/cd))*vel_indep(i); %ft*lb/s
end 

v_vclimb = (PA/W0) - 0.8776*sqrt((W0/s)/(rho_SSL*cd0))*(1/(ld_max)^1.5); %max R/C, ft/s
v_PR_min = v_E_max; %max RC occurs at pr_min, and prmin=v_E_max

%time to climb and abs flight ceiling
ceiling_alt = -19867*log((W0/PA)*sqrt((2*W0)/(rho_SSL*s))*((0.7436*(cd0)^(1/4))/((e_oswald*AR)^(3/4)))); %eqn 6.15.8, p498 %ft
a_ttc = v_vclimb; %ft/s
b_ttc = v_vclimb/ceiling_alt; %1/s
ttc_inter = log(a_ttc-(b_ttc*10000)) - log(a_ttc); %10000ft requested 
ttc = ((-1)/b_ttc)*ttc_inter; %s

R_glide_max = h*ld_max; %ft

%% Output

fprintf('Temperature at an altitude of %0.0fft = %0.2f °F\n',h_g, T_alt_F);
fprintf('Pressure at an altitude of %0.0fft = %0.2f psf\n',h_g, P_alt);
fprintf('Density at an altitude of %0.0fft = %0.6f sl/ft^3\n',h_g, rho_alt);
fprintf('Wing aspect ratio = %0.2f\n',AR);
fprintf('Spanwise efficiency factor = %0.4f\n',e_wing);
fprintf('Oswald efficiency factor = %0.4f\n',e_oswald);
fprintf('L/D ratio, max = %0.2f\n',ld_max);
fprintf('Reynolds number from wing MAC and flight speed at max L/D = %0.2f\n',Re);
fprintf('Wing loading at 50%% fuel consumed = %0.2f lb/ft^s\n',wing_loading);
fprintf('Stall speed at 50%% fuel consumed = %0.2f ft/s\n',v_stall);
fprintf('Max range and corresponding speed at 50%% fuel consumed = %0.2f ft = %0.2f mi at %0.2f ft/s\n',R_max,R_max/5280,v_R_max);
fprintf('Max endurance and corresponding speed at 50%% fuel consumed = %0.2f s = %0.2f hr at %0.2f ft/s\n',E_max,E_max/3600,v_E_max);
fprintf('Max speed and corresponding Mach number at 50%% fuel consumed = %0.2f ft/s = %0.2f Mach\n',v_max,M_max);
fprintf('Absolute flight ceiling at MTOW = %0.2f ft\n',ceiling_alt);
fprintf('Max rate of climb at SSL = %0.2f ft/s\n',v_vclimb);
fprintf('Time to climb to 10,000ft at MTOW = %0.2f s = %0.2f min\n',ttc,ttc/60);
fprintf('Max glide range = %0.2f ft = %0.2f mi\n',R_glide_max,R_glide_max/5280);

%% Plotting, Final Output

scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]) 
plot(vel_indep,PR/550,'DisplayName','Power Required Curve'); %PR in hp
hold on
grid on

%indicate the power required for the following at 50% fuel consumed
v_indicate = [v_stall v_R_max v_E_max v_PR_min]; %ft/s
PR_indicate = zeros(1,length(v_indicate));
for i = 1:length(v_indicate) %iterate velocity 
    q = 0.5*rho_alt*v_indicate(i)*v_indicate(i);
    cl = W0/(s*q);
    cd = cd0 + (cl*cl*k); 
    PR_indicate(i) = (W0/(cl/cd))*v_indicate(i); %ft*lb/s
end 
PR_indicate = PR_indicate/550; %hp
scatter(v_indicate(1),PR_indicate(1),75,'o','filled','DisplayName','Stall Velocity');
scatter(v_indicate(2),PR_indicate(2),75,'^','filled','DisplayName','Max Range');
scatter(v_indicate(3),PR_indicate(3),75,'d','DisplayName','Max Endurace');
scatter(v_indicate(4),PR_indicate(4),75,'s','filled','DisplayName','Max Rate of Climb');

legend('Location','bestoutside')
title("Figure 1: Power Required vs. Speed at 50% fuel consumed");
xlabel("Speed (ft/s)");
ylabel("Power Required (hp)");
hold off

fprintf('Plot of power required vs. speed at 50%% fuel consumed has been generated.\n'); 
fprintf('\n');
