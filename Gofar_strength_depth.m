close all
clear all
clc


%% Data
load('Gofar_T_interp_GSRM.mat');
load('Gofar_strain_rate.mat');
%%


%% Parameters
z = gofar_z_interp.*1e3; %m
T = gofar_T_interp + 273.15; %deg K
density_rock = 3300; %kg/m^3
density_water = 1000; %kg/m^3
g = 9.81; %m/s^2
P = density_rock.*g.*z; %Pressure Pa 
P_f_hydro = 0.4.*P; %Hydrostatic pore pressure
P_f_litho = 0.95.*P; %Lithostatic pore pressure
R = 8.314; %Ideal gas constant J/K/mol  

d_protolith_upper = 5000; %Grain size microns
d_protolith_lower = 1000; %Grain size microns

d_mylonite_upper = 10; %Grain size microns
d_mylonite_lower = 1; %Grain size microns
%%


%% Protolith Deformation
% Dry Dislocation Creep (Hirth & Kohlstedt, 2003)
A = 1.096e5;  %MPa^(1/n)*micron^(1/d)*s^-1   
n = 3.5;            
Q = 520e3;    %J/mol    
V = 17e-6;    %m^3*mol^-1     
p = 0; 
sigma_ol_dry_dis_upper = (strain_rate_proto_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_protolith_upper.^p).*A)).^(1./n);
sigma_ol_dry_dis_lower = (strain_rate_proto_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_protolith_lower.^p).*A)).^(1./n);

% Dry Grain Boundary Sliding (Hansen et al. 2011)
A = 6.310e4;  %MPa^(1/n)*micron^(1/d)*s^-1           
n = 2.9;            
Q = 445e3;    %J/mol        
V = 14e-6;    %m^3*mol^-1     
p = -0.7;        
sigma_ol_dry_gbs_upper = (strain_rate_proto_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_protolith_upper.^p).*A)).^(1./n);
sigma_ol_dry_gbs_lower = (strain_rate_proto_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_protolith_lower.^p).*A)).^(1./n);

% Dry Low Temperature Plasticity (Mei et al. 2010)
sigma_p = 5900e6; %Periels strss Pa
H = 5.4e5; %J/mol
B = 5.7e-11; %s^-1
q = 2;
sigma_ol_dry_ltp_upper = real((sigma_p.*(1-(((-R.*T)./H).*log(strain_rate_proto_lower./B)).^(1/q))));
sigma_ol_dry_ltp_lower = real((sigma_p.*(1-(((-R.*T)./H).*log(strain_rate_proto_upper./B)).^(1/q))));
%%


%% HT Mylonite
T_lower = 725 + 273.15; %deg K
T_upper = 775 + 273.15; %deg K
depth_lower = 4.5; %Corresponding depth from thermal model km
depth_upper = 4.9; %Corresponding depth from thermal model km
P_lower = (depth_lower.*1e3)*(density_rock-density_water)*g;
P_upper = (depth_upper.*1e3)*(density_rock-density_water)*g; 

a1=-1.194;
b1=2.263;
c1=0.128;
CH2O_lower=10.^(a1+b1.*(depth_lower.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_lower=round(CH2O_lower.*16.3);         %COH in H/10^6Si
CH2O_upper=10.^(a1+b1.*(depth_upper.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_upper=round(CH2O_upper.*16.3);         %COH in H/10^6Si

% Wet Diffsion Creep (Hirth & Kohlstedt, 2003)
A = 1e6;      %MPa^(1/n)*micron^(1/d)*H/10^6Si^(1/4)*s^-1
A_corrected = 4e5;
n = 1;
r = 1;
Q = 335e3;    %J/mol      
V = 21e-6;    %m^3*mol^-1      
p = -3;  
sigma_ol_wet_dif_ht_upper = (strain_rate_ht_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_upper.^p).*(COH_upper.^r).*A_corrected));
sigma_ol_wet_dif_ht_lower = (strain_rate_ht_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_lower.^p).*(COH_lower.^r).*A_corrected));
%%


%% MT Mylonite
T_lower = 575 + 273.15; %deg K
T_upper = 625 + 273.15; %deg K
depth_lower = 3.5; %Corresponding depth from thermal model km
depth_upper = 3.8; %Corresponding depth from thermal model km
P_lower = (depth_lower.*1e3)*(density_rock-density_water)*g;
P_upper = (depth_upper.*1e3)*(density_rock-density_water)*g; 

a1=-1.194;
b1=2.263;
c1=0.128;
CH2O_lower=10.^(a1+b1.*(depth_lower.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_lower=round(CH2O_lower.*16.3);         %COH in H/10^6Si
CH2O_upper=10.^(a1+b1.*(depth_upper.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_upper=round(CH2O_upper.*16.3);         %COH in H/10^6Si

% Wet Diffsion Creep (Hirth & Kohlstedt, 2003)
A = 1e6;      %MPa^(1/n)*micron^(1/d)*H/10^6Si^(1/4)*s^-1
A_corrected = 4e5;
n = 1;
r = 1;
Q = 335e3;    %J/mol      
V = 21e-6;    %m^3*mol^-1      
p = -3;  
sigma_ol_wet_dif_mt_upper = (strain_rate_mt_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_upper.^p).*(COH_upper.^r).*A_corrected));
sigma_ol_wet_dif_mt_lower = (strain_rate_mt_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_lower.^p).*(COH_lower.^r).*A_corrected));
%%


%% LT Mylonite
T_lower = 475 + 273.15; %deg K
T_upper = 525 + 273.15; %deg K
depth_lower = 3;   %Corresponding depth from thermal model km
depth_upper = 3.3; %Corresponding depth from thermal model km
P_lower = (depth_lower.*1e3)*(density_rock-density_water)*g; 
P_upper = (depth_upper.*1e3)*(density_rock-density_water)*g; 

a1=-1.194;
b1=2.263;
c1=0.128;
CH2O_lower=10.^(a1+b1.*(depth_lower.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_lower=round(CH2O_lower.*16.3);         %COH in H/10^6Si
CH2O_upper=10.^(a1+b1.*(depth_upper.^c1)); %Water content in ppm (Hirschmann et al., 2005)
COH_upper=round(CH2O_upper.*16.3);         %COH in H/10^6Si

% Wet Diffsion Creep (Hirth & Kohlstedt, 2003)
A = 1e6;      %MPa^(1/n)*micron^(1/d)*H/10^6Si^(1/4)*s^-1
A_corrected = 4e5;
n = 1;
r = 1;
Q = 335e3;    %J/mol      
V = 21e-6;    %m^3*mol^-1      
p = -3;  
sigma_ol_wet_dif_lt_upper = (strain_rate_lt_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_upper.^p).*(COH_upper.^r).*A_corrected));
sigma_ol_wet_dif_lt_lower = (strain_rate_lt_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_lower.^p).*(COH_lower.^r).*A_corrected));

% Antigorite Power Law (Hilariet et al. 2007)
A = 1.8e-17;    %Recalculated from Amiguet et al. (2012)
n = 3.8;        %Stress exponent
Q = 8.9e3;      %J/mol        
V = 3.2e-6;     %m^3*mol^-1     
p = 0;           
sigma_atg_pl_lt_upper = (strain_rate_lt_serp_lower./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_upper.^p).*A)).^(1./n);
sigma_atg_pl_lt_lower = (strain_rate_lt_serp_upper./(exp(-(Q+P.*V)./(R.*T)).*(d_mylonite_lower.^p).*A)).^(1./n);
%%


%% Brittle Deformation
sigma_n = P; %Normal stress
sigma_t = 21e6; %Tensile strength Pa (Boettcher et al. 2007) 
sigma_y = (1./sigma_ol_dry_ltp_lower + 1./(sigma_ol_dry_gbs_lower.*1e6)).^(-1); %Yield stress
sigma_eff = (sigma_n-P_f_hydro)./(1-(P_f_hydro./sigma_y)); %Effective stress hydrostatic
sigma_eff_litho = (sigma_n-P_f_litho)./(1-(P_f_litho./sigma_y)); %Effective stress lithostatic
alpha_hydro = 1 - sigma_eff./sigma_y; %Pore fluid factor hydrostatic
alpha_litho = 1 - sigma_eff_litho./sigma_y; %Pore fluid factor lithostatic

% Byerlee's law
mu_i1 = find(sigma_n < 200e6,1,'last');

mu = ones(length(P),1);
mu(1:mu_i1,1) = 0.85;
mu(mu_i1:end,1) = 0.6;

F = (2.*mu)./((mu.^2 + 1).^(1/2));
sigma_brittle_hydro = F.*(P) + sigma_t; %alpha = 0
sigma_brittle_hydro(mu_i1:end) = sigma_brittle_hydro(mu_i1:end) + sigma_brittle_hydro(mu_i1-1)-sigma_brittle_hydro(mu_i1) + sigma_brittle_hydro(mu_i1-1)-sigma_brittle_hydro(mu_i1-2); %Pa

sigma_brittle_eff_hydro = F.*(P-alpha_hydro.*P_f_hydro) + sigma_t;
sigma_brittle_eff_hydro(mu_i1:end) = sigma_brittle_eff_hydro(mu_i1:end) + sigma_brittle_eff_hydro(mu_i1-1)-sigma_brittle_eff_hydro(mu_i1) + sigma_brittle_eff_hydro(mu_i1-1)-sigma_brittle_eff_hydro(mu_i1-2); %Pa

sigma_brittle_eff_hydro_a = F.*(P-P_f_hydro) + sigma_t; %alpha = 1
sigma_brittle_eff_hydro_a(mu_i1:end) = sigma_brittle_eff_hydro_a(mu_i1:end) + sigma_brittle_eff_hydro_a(mu_i1-1)-sigma_brittle_eff_hydro_a(mu_i1) + sigma_brittle_eff_hydro_a(mu_i1-1)-sigma_brittle_eff_hydro_a(mu_i1-2); %Pa

sigma_brittle_eff_litho = F.*(P-alpha_litho.*P_f_litho) + sigma_t;
sigma_brittle_eff_litho(mu_i1:end) = sigma_brittle_eff_litho(mu_i1:end) + sigma_brittle_eff_litho(mu_i1-1)-sigma_brittle_eff_litho(mu_i1) + sigma_brittle_eff_litho(mu_i1-1)-sigma_brittle_eff_litho(mu_i1-2); %Pa

% Serpentine friction
mu_i1 = find(T - 273.15 < 100,1,'last');
mu_i2 = find(T - 273.15 < 300,1,'last');

mu_weak_lb = ones(length(P),1);
mu_weak_lb(1:mu_i1,1) = 0.1;
mu_weak_lb(mu_i1:mu_i2,1) = 0.3;
mu_weak_lb(mu_i2:end,1) = 0.6;

F_weak_lb = (2.*mu_weak_lb)./((mu_weak_lb.^2 + 1).^(1/2));
sigma_brittle_weak_lb = F_weak_lb.*(P) + sigma_t; %alpha = 0
sigma_brittle_weak_lb(mu_i1:end) = sigma_brittle_weak_lb(mu_i1:end) + sigma_brittle_weak_lb(mu_i1-1)-sigma_brittle_weak_lb(mu_i1) + sigma_brittle_weak_lb(mu_i1-1)-sigma_brittle_weak_lb(mu_i1-2); %Pa
sigma_brittle_weak_lb(mu_i2:end) = sigma_brittle_weak_lb(mu_i2:end) + sigma_brittle_weak_lb(mu_i2-1)-sigma_brittle_weak_lb(mu_i2) + sigma_brittle_weak_lb(mu_i2-1)-sigma_brittle_weak_lb(mu_i2-2); %Pa

sigma_brittle_eff_weak_lb = F_weak_lb.*(P-alpha_hydro.*P_f_hydro) + sigma_t; %alpha = 1
sigma_brittle_eff_weak_lb(mu_i1:end) = sigma_brittle_eff_weak_lb(mu_i1:end) + sigma_brittle_eff_weak_lb(mu_i1-1)-sigma_brittle_eff_weak_lb(mu_i1) + sigma_brittle_eff_weak_lb(mu_i1-1)-sigma_brittle_eff_weak_lb(mu_i1-2); %Pa
sigma_brittle_eff_weak_lb(mu_i2:end) = sigma_brittle_eff_weak_lb(mu_i2:end) + sigma_brittle_eff_weak_lb(mu_i2-1)-sigma_brittle_eff_weak_lb(mu_i2) + sigma_brittle_eff_weak_lb(mu_i2-1)-sigma_brittle_eff_weak_lb(mu_i2-2); %Pa

mu_weak_ub = ones(length(P),1);
mu_weak_ub(1:mu_i1,1) = 0.4;
mu_weak_ub(mu_i1:end,1) = 0.6;

F_weak_ub = (2.*mu_weak_ub)./((mu_weak_ub.^2 + 1).^(1/2));
sigma_brittle_weak_ub = F_weak_ub.*(P) + sigma_t; %alpha = 0
sigma_brittle_weak_ub(mu_i1:end) = sigma_brittle_weak_ub(mu_i1:end) + sigma_brittle_weak_ub(mu_i1-1)-sigma_brittle_weak_ub(mu_i1) + sigma_brittle_weak_ub(mu_i1-1)-sigma_brittle_weak_ub(mu_i1-2); %Pa

sigma_brittle_eff_weak_ub = F_weak_ub.*(P-alpha_hydro.*P_f_hydro) + sigma_t; %alpha = 1
sigma_brittle_eff_weak_ub(mu_i1:end) = sigma_brittle_eff_weak_ub(mu_i1:end) + sigma_brittle_eff_weak_ub(mu_i1-1)-sigma_brittle_eff_weak_ub(mu_i1) + sigma_brittle_eff_weak_ub(mu_i1-1)-sigma_brittle_eff_weak_ub(mu_i1-2); %Pa
%%


%% Strength Depth Plot
figure
box on
hold on

plot(z./1e3,real(sigma_brittle_eff_weak_lb)./1e6,'-','Color',[0 0 0]);
plot(z./1e3,real(sigma_brittle_eff_hydro)./1e6,'-','Color',[0 0 0]);

plot(z./1e3,sigma_ol_dry_gbs_upper,'-','Color',[153 153 153]/255);
plot(z./1e3,sigma_ol_dry_gbs_lower,'-','Color',[153 153 153]/255);

plot(z./1e3,sigma_ol_wet_dif_ht_upper,'-','Color',[227 26 28]/255);
plot(z./1e3,sigma_ol_wet_dif_ht_lower,'-','Color',[227 26 28]/255);

plot(z./1e3,sigma_ol_wet_dif_mt_upper,'-','Color',[253 141 60]/255);
plot(z./1e3,sigma_ol_wet_dif_mt_lower,'-','Color',[253 141 60]/255);

plot(z./1e3,sigma_ol_wet_dif_lt_upper,'-','Color',[254 204 92]/255);
plot(z./1e3,sigma_ol_wet_dif_lt_lower,'-','Color',[254 204 92]/255);

plot(z./1e3,sigma_atg_pl_lt_upper,'-','Color',[44 162 95]/255);
plot(z./1e3,sigma_atg_pl_lt_lower,'-','Color',[44 162 95]/255);

pbaspect([1 1 1])
xlabel('Depth (km)')
ylabel('Differential stress (MPa)')
set(gca,'FontSize',7.2)
set(gca,'YColor','k')
set(gca,'XColor','k')

axis([0 12 0 400])
set(gca,'YTick',[0:100:500])
set(gca,'XTick',[0:2:20])
set(gca, 'YAxisLocation', 'left')
set(gca,'TickLength',[0.008, 0.01])
set(gca,'TickDir','out')
view(90,90)
%%


%% Strength Alpha Depth Plot
figure
plot(sigma_brittle_eff_hydro./1e6,z./1e3,'-','Color',[0 0 0]);
hold on
plot(sigma_brittle_eff_litho./1e6,z./1e3,'--','Color',[0 0 0]);
plot(sigma_brittle_eff_hydro_a./1e6,z./1e3,'-','Color',[0.75 0.75 0.75]);
plot(sigma_brittle_hydro./1e6,z./1e3,'-','Color',[0.75 0.75 0.75]);
plot(sigma_y./1e6,z./1e3,'-','Color',[227 26 28]/255);

pbaspect([1 1 1])
xlabel('Differential stress (MPa)')
ylabel('Depth (km)')
set(gca,'FontSize',7.2)
set(gca,'YColor','k')
set(gca,'XColor','k')

axis([0 400 0 12])
set(gca,'YDir','reverse')
set(gca,'XTick',[0 100 200 300 400 500])
set(gca,'YTick',[0:2:20])
set(gca, 'YAxisLocation', 'left')
set(gca,'TickLength',[0.008, 0.01])
set(gca,'TickDir','out')

ax1 = gca;
ax1.XColor = 'k';
ax1.YColor = 'k';

ax1_pos = ax1.Position;
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');

ax2 = gca;
hold on
plot(ax2,alpha_hydro,z./1e3,'-','Color','b');
plot(ax2,alpha_litho,z./1e3,'--','Color','b');

pbaspect([1 1 1])
set(gca,'FontSize',7.2)
set(gca,'YDir','reverse')
set(gca,'TickLength',[0.008, 0.01])
set(gca,'TickDir','out')
ax2.XColor = 'k';
ax2.YColor = 'k';

axis([0 1 0 15])
yticks([])
xticks([0:0.2:1])
xlabel('Pore fluid factor, a')