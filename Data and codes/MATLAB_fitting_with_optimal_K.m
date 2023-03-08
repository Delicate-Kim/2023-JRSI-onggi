
clear all;
close all;
clc;
%% paramters

V_g = 1.9e-3; %volume [m^3]
V_o = 4.6e-3;
A_g = 4.71e-2; %side area [m^2]
A_o = 6.28e-2;
d_g = 5.5e-3; %thickness [m]
d_o = 5.5e-3; 

P_0 = 101.325; % ambient pressure [kPa]
T_0 = 298; % ambient temperature [K]
R = 8.31e-3; % gas constant [kPa m^3 mol^-1 K^-1]
mu_CO2 = 1.47e-8; % dynamic viscosity [kPa s]
mu_O2 = 1.95e-8;


%% read data

load Headspace_pure_and_salted_48hours_data.mat;

pure = table2array(pure_final_48hrs); % time(hrs) | glass co2(ppm) ave | std | onggi co2(ppm) ave | std | glass pressure(kPa) ave | std | onggi pressure(kPa) ave | std 
salted = table2array(salted_final_48hrs);

Np = length(pure(:,1));
Ns = length(salted(:,1));

pp_co2gp = pure(:,2); %partial pressure, glass
pp_co2gp_std = pure(:,3);
pp_co2op = pure(:,4); %onggi
pp_co2op_std = pure(:,5);

pp_co2gs = salted(:,2); %partial pressure, glass
pp_co2gs_std = salted(:,3);
pp_co2os = salted(:,4); %onggi
pp_co2os_std = salted(:,5);

time = pure(:,1);

%% Data fitting - find optimal k_bar for onggi and glass
% MATLAB fitting options: https://www.mathworks.com/help/curvefit/fitoptions.html
% current mehtod: nonlinearLeastSquares -> minimize rmse
%onggi
K_bar_oopt = 1e-8;
R_ave_oopt = 0; % random value which is big enough

for K_bar_temp=K_bar_oopt:1e-8:10001e-8
    tau_temp = (V_o*d_o)/(R*T_0*K_bar_temp*A_o);

    % co2op
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[pp_co2op(1,1),-inf,tau_temp],...
                   'Upper',[pp_co2op(1,1),inf,tau_temp],...
                   'Startpoint',[pp_co2op(1,1) 1 10]);
    f = fittype('P_0+a*(1-exp(-x/tau))','option', s);
    [c_co2op_temp, gof_co2op_temp] = fit(pure(:,1), pp_co2op(:,1),f);
    
    n_dot_co2op_temp = K_bar_temp*A_o*c_co2op_temp.a/d_o;
    R_avep_temp = gof_co2op_temp.rsquare;
    %rmse_avep_temp = gof_co2op_temp.rmse;

    % co2os
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[pp_co2os(1,1),-inf,tau_temp],...
                   'Upper',[pp_co2os(1,1),inf,tau_temp],...
                   'Startpoint',[pp_co2os(1,1) 1 10]);
    f = fittype('P_0+a*(1-exp(-x/tau))','option', s);
    [c_co2os_temp, gof_co2os_temp] = fit(pure(:,1), pp_co2os(:,1),f);
    
    n_dot_co2os_temp = K_bar_temp*A_o*c_co2os_temp.a/d_o;
    R_aves_temp = gof_co2os_temp.rsquare;

    R_ave_temp = (R_avep_temp+R_aves_temp)/2; 


    if R_ave_temp > R_ave_oopt && R_aves_temp > 0.85 % minimum bound for R-squared
        K_bar_oopt = K_bar_temp;
        R_ave_oopt = R_ave_temp;
        n_dot_co2op_opt = n_dot_co2op_temp;
        n_dot_co2os_opt = n_dot_co2os_temp;
        gof_co2op_opt = gof_co2op_temp;
        gof_co2os_opt = gof_co2os_temp;
    end
end

%glass
K_bar_gopt = 1e-8;
R_ave_gopt = 0; % random value which is big enough

for K_bar_temp=K_bar_gopt:1e-8:1001e-8
    tau_temp = (V_g*d_g)/(R*T_0*K_bar_temp*A_g);

    % co2gp
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[pp_co2gp(1,1),-inf,tau_temp],...
                   'Upper',[pp_co2gp(1,1),inf,tau_temp],...
                   'Startpoint',[pp_co2gp(1,1) 1 10]);
    f = fittype('P_0+a*(1-exp(-x/tau))','option', s);
    [c_co2gp_temp, gof_co2gp_temp] = fit(pure(:,1), pp_co2gp(:,1),f);
    
    n_dot_co2gp_temp = K_bar_temp*A_g*c_co2gp_temp.a/d_g;
    R_avep_temp = gof_co2gp_temp.rsquare;

    % co2gs
    s = fitoptions('Method','NonlinearLeastSquares',...
                   'Lower',[pp_co2gs(1,1),-inf,tau_temp],...
                   'Upper',[pp_co2gs(1,1),inf,tau_temp],...
                   'Startpoint',[pp_co2gs(1,1) 1 10]);
    f = fittype('P_0+a*(1-exp(-x/tau))','option', s);
    [c_co2gs_temp, gof_co2gs_temp] = fit(pure(:,1), pp_co2gs(:,1),f);
    
    n_dot_co2gs_temp = K_bar_temp*A_g*c_co2gs_temp.a/d_g;
    R_aves_temp = gof_co2gs_temp.rsquare;

    R_ave_temp = (R_avep_temp+R_aves_temp)/2; 

    if R_ave_temp > R_ave_gopt && R_aves_temp > 0.95 % minimum bound for R-squared
        K_bar_gopt = K_bar_temp;
        R_ave_gopt = R_ave_temp;
        n_dot_co2gp_opt = n_dot_co2gp_temp;
        n_dot_co2gs_opt = n_dot_co2gs_temp;
        gof_co2gp_opt = gof_co2gp_temp;
        gof_co2gs_opt = gof_co2gs_temp;
    end
end

%% plot

figure(1) %pure
%title('time vs CO2 partial pressure');
hold on;
%grid on;

coord_up = [time,pp_co2op(:,1)+pp_co2op_std(:,1)];
coord_low = [time,pp_co2op(:,1)-pp_co2op_std(:,1)];
coord_combine = [coord_up;flipud(coord_low)];
fill(coord_combine(:,1),coord_combine(:,2),'b', 'FaceColor', '#e8a9a5', 'EdgeColor','none','FaceAlpha',1)

p1 = plot(time,pp_co2op(:,1), 'LineWidth', 5, 'Color','#de382c');

coord_up = [time,pp_co2gp(:,1)+pp_co2gp_std(:,1)];
coord_low = [time,pp_co2gp(:,1)-pp_co2gp_std(:,1)];
coord_combine = [coord_up;flipud(coord_low)];
fill(coord_combine(:,1),coord_combine(:,2),'b', 'FaceColor', '#b6cdf2', 'EdgeColor','none','FaceAlpha',1)


p2 = plot(time,pp_co2gp(:,1), 'LineWidth', 5, 'Color','#2c62b8');


%x = linspace(0,28,280);
x = linspace(0, Np/120, 48*Np/120); % 2 data per minute, 120 data per hr

y_co2gp_opt = pp_co2gp(1,1) +(n_dot_co2gp_opt*d_g/K_bar_gopt/A_g)-(n_dot_co2gp_opt*d_g/K_bar_gopt/A_g)*exp(-x*R*T_0*K_bar_gopt*A_g/V_g/d_g);
y_co2op_opt = pp_co2op(1,1) +(n_dot_co2op_opt*d_o/K_bar_oopt/A_o)-(n_dot_co2op_opt*d_o/K_bar_oopt/A_o)*exp(-x*R*T_0*K_bar_oopt*A_o/V_o/d_o);

p3 = plot(x,y_co2gp_opt, '--k', 'LineWidth', 3);
p4 = plot(x,y_co2op_opt, '--k', 'LineWidth', 3);
p3.Color(4) = 0.5;
p4.Color(4) = 0.5;

set(gca, 'FontSize', 40)
xlabel('Time (hours)', 'FontSize',35);
ylabel('CO_2 Partial pressure (kPa)', 'FontSize',35);
% ay = ancestor(p1, 'axes');
% ay.YAxis.Exponent = 0;
% ytickformat('%,.0f');
%title('CO2 generation of Salted Cabbage')
%legend([p1 p2 p3],'Onggi', 'Glass', 'Model')


xlim([0 48]);
set(gca,'Xtick',0:4:48)
xticklabels({0, '', 8, '', 16, '', 24, '', 32, '', 40, '', 48})
set(gca,'Ytick',0:0.25:3)
yticklabels({0, '', 0.5, '', 1, '', 1.5, '', 2, '', 2.5, '', 3})

%ax = ancestor(p3, 'axes');
%ax.XAxis.TickVales = [];

%%
figure(2) %salted

hold on;
%grid on;

coord_up = [time,pp_co2os(:,1)+pp_co2os_std(:,1)];
coord_low = [time,pp_co2os(:,1)-pp_co2os_std(:,1)];
coord_combine = [coord_up;flipud(coord_low)];
fill(coord_combine(:,1),coord_combine(:,2),'b', 'FaceColor', '#e8a9a5', 'EdgeColor','none','FaceAlpha',1)

p1 = plot(time,pp_co2os(:,1), 'LineWidth', 5, 'Color','#de382c');

coord_up = [time,pp_co2gs(:,1)+pp_co2gs_std(:,1)];
coord_low = [time,pp_co2gs(:,1)-pp_co2gs_std(:,1)];
coord_combine = [coord_up;flipud(coord_low)];
fill(coord_combine(:,1),coord_combine(:,2),'b', 'FaceColor', '#b6cdf2', 'EdgeColor','none','FaceAlpha',1)

p2 = plot(time,pp_co2gs(:,1), 'LineWidth', 5, 'Color','#2c62b8');


x = linspace(0, Np/120, 48*Np/120); % 2 data per minute, 120 data per hr

y_co2gs_opt = pp_co2gs(1,1) +(n_dot_co2gs_opt*d_g/K_bar_gopt/A_g)-(n_dot_co2gs_opt*d_g/K_bar_gopt/A_g)*exp(-x*R*T_0*K_bar_gopt*A_g/V_g/d_g);
y_co2os_opt = pp_co2os(1,1) +(n_dot_co2os_opt*d_o/K_bar_oopt/A_o)-(n_dot_co2os_opt*d_o/K_bar_oopt/A_o)*exp(-x*R*T_0*K_bar_oopt*A_o/V_o/d_o);

p3 = plot(x,y_co2gs_opt, '--k', 'LineWidth', 3);
p4 = plot(x,y_co2os_opt, '--k', 'LineWidth', 3);
p3.Color(4) = 0.5;
p4.Color(4) = 0.5;

set(gca, 'FontSize', 40)
xlabel('Time (hours)', 'FontSize',35);
ylabel('CO_2 Partial pressure (kPa)', 'FontSize',35);
ay = ancestor(p1, 'axes');
ay.YAxis.Exponent = 0;
ytickformat('%,.0f');
%title('CO2 generation of Salted Cabbage')
%legend([p1 p2 p3],'Onggi', 'Glass', 'Model')

xlim([0 48]);
set(gca,'Xtick',0:4:48)
xticklabels({0, '', 8, '', 16, '', 24, '', 32, '', 40, '', 48})
set(gca,'Ytick',0:1:10)
yticklabels({0, '', 2, '', 4, '', 6, '', 8, '', 10})



