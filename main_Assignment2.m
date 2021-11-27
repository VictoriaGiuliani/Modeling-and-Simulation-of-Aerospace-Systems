%MSAS Assignment 2
% Author: Victoria Katia Giuliani

%%
%For plots

set(groot,'defaulttextinterpreter','Latex');
set(groot,'defaultAxesTickLabelInterpreter','Latex');
set(groot,'defaultLegendInterpreter','Latex');
set(0,'defaultAxesFontSize',20);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exercise 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; 
%Data:
J1 = 0.2; %kgm
J2 = 0.1; %kgm
T0 = 0.1; %Nm
t0 = 0; %s
tf = 10;% 10;  %s
f = 100; %Hz
k = 0.5;%0.1; %3.1
b = 1; %0.5;  %2.7
time = [t0 : 1/f : tf];

%Compute eigenvalues by considering the arbitrary values of k and b
A_r = [0 0 1 0; 0 0 0 1; -k/J1 k/J1 0 0; k/J2 -k/J2 0 -2*b/J2];
lambda_r = eig(A_r);
fprintf('Exercise 1\n\n');
fprintf('The eigenvalues with random k and b are:\n');
fprintf('lambda_1 = %f + %fi\n', real(lambda_r(1)), imag(lambda_r(1)));
fprintf('lambda_2 = %f + %fi\n', real(lambda_r(2)), imag(lambda_r(2)));
fprintf('lambda_3 = %f + %fi\n', real(lambda_r(3)), imag(lambda_r(3)));
fprintf('lambda_4 = %f + %fi\n\n', real(lambda_r(4)), imag(lambda_r(4)));

x0 = [0, 0, 0, 0]'; %initial condition at rest
load('samples')
%Integrate from t0 to tf --> system response with guess values for k and b
[~, x] = ode45(@(t, x)reactionwheel(t, x, J1, J2, k, b, T0), time, x0);

%Plot the behaviour of the system
figure('Name', 'Exercise 1', 'WindowState', 'maximized')
plot(time, x(:, 1)*180/pi, 'b');
hold on
plot(time, x(:, 2)*180/pi, 'm');
grid on
title('Rotation of the disks (with arbitrary value for k and b)');
legend('$\theta_1$', '$\theta_2$', 'Location', 'best');
xlabel('Time [s]'); ylabel('Angle of rotation [deg]');
figure('Name','Exercise 1', 'WindowState', 'maximized')
plot(time, x(:, 3), 'b');
hold on
plot(time, x(:, 4), 'm');
grid on
title('Angular velocity of the disks (with arbitrary value for k and b)');
legend('$\dot{\theta_1}$', '$\dot{\theta_2}$', 'Location', 'Best');
xlabel('Time [s]'); ylabel('Angular velocity [rad/s]');

%Point 3
dd_theta1 = k/J1*(x(:, 2) - x(:, 1));
dd_theta2 = -k/J2*(x(:, 2) - x(:, 1)) - b/J2*sign(x(:, 4)).*x(:, 4).^2 + T0/J2*ones(size(x(:, 4)));
%Plot experimental data against numerical results for first guess of k and b
figure('Name', 'Exercise 1', 'WindowState', 'maximized')
scatter(time, samples(:, 2), 8, 'g');
hold on
plot(time,dd_theta1, 'b', 'LineWidth', 1.5);
title('Angular acceleration of disk 1');
xlabel('Time [s]'); ylabel('Acceleration [rad/s$^2$]');
legend('Sampled acceleration', 'Acceleration with random k and b');
grid on

figure('Name', 'Exercise 1',  'WindowState', 'maximized')
scatter(time, samples(:, 3), 8, 'g');
hold on
plot(time, dd_theta2, 'm', 'LineWidth', 1.5);
title('Angular acceleration of disk 2');
xlabel('Time [s]'); ylabel('Acceleration [rad/s$^2$]');
legend('Sampled acceleration', 'Acceleration with random k and b');
grid on

%Point 3: determine value of k and b that allows retracing the experimental
%data
k0 = k;
b0 = b; 
kb0 = [k0; b0];
acc_sample1 = samples(:, 2);
acc_sample2 = samples(:, 3);
fun = @(var) costfunction(var, J1, J2, T0, time, x0, 1, acc_sample1);
options = optimset('Algorithm', 'levenberg-marquardt');
[kb, J] = fsolve(fun, kb0, options); %find the right values of k and b
%Plot the solution to see if experimental data and acceleration with the
%found values of k and b actually 
K = kb(1);
B = kb(2);
A = [0 0 1 0; 0 0 0 1; -K/J1 K/J1 0 0; K/J2 -K/J2 0 -2*B/J2];%-2*B/J2
lambda = eig(A);
fprintf('\nValue of the parameters retracing the experimental data:\n');
fprintf('k = %f\n', K);
fprintf('b = %f\n', B);
fprintf('Accuracy = %f\n', J);
fprintf('The eigenvalues with exact k and b are:\n');
fprintf('lambda_1 = %f + %fi\n', real(lambda(1)), imag(lambda(1)));
fprintf('lambda_2 = %f + %fi\n', real(lambda(2)), imag(lambda(2)));
fprintf('lambda_3 = %f + %fi\n', real(lambda(3)), imag(lambda(3)));
fprintf('lambda_4 = %f + %fi\n\n', real(lambda(4)), imag(lambda(4)));
[~, x_kb] = ode45(@(t, x)reactionwheel(t, x, J1, J2, K, B, T0), time, x0);
dd_theta1_kb = K/J1*(x_kb(:, 2) - x_kb(:, 1));
dd_theta2_kb = -K/J2*(x_kb(:, 2) - x_kb(:, 1)) - B/J2*sign(x_kb(:, 4)).*x_kb(:, 4).^2 + T0/J2*ones(size(x_kb(:, 4)));
%Plots
figure('Name','Exercise 1', 'WindowState', 'maximized');
plot(time, acc_sample1, 'og');
hold on
plot(time, dd_theta1_kb, 'r', 'LineWidth', 1.5);
grid on
title('Angular acceleration of disk 1');
kbstr = sprintf('Acceleration with k = %.3f, b = %.3f', K, B);
legend('Sampled Acceleration', kbstr);
xlabel('Time [s]'); ylabel('Acceleration [rad/s$^2$]');

figure('Name','Exercise 1', 'WindowState', 'maximized');
scatter(time, acc_sample2, 'g');
hold on
plot(time, dd_theta2_kb, 'r', 'LineWidth', 1.5);
grid on
title('Angular acceleration of disk 2');
xlabel('Time [s]'); ylabel('Acceleration [rad/s$^2$]');
kbstr = sprintf('Acceleration with k = %.3f, b = %.3f', K, B);
legend('Sampled Acceleration', kbstr)
ax=axes;
set(ax,'units','normalized','position',[0.2,0.6,0.3,0.3])
box(ax,'on')
plot(time, acc_sample2, 'og', 'MarkerFaceColor', 'g', 'parent', ax);
hold on
plot(time, dd_theta2_kb, 'r','LineWidth', 1.5, 'parent', ax)
grid on
set(ax,'xlim',[0,0.1],'ylim',[0.4, 1.1], 'FontSize', 10)
%Plot angular rotation
figure('Name', 'Exercise 1', 'Windowstate', 'maximized')
plot(time, x_kb(:, 1)*180/pi, time, x_kb(:, 2)*180/pi, 'linewidth', 1.5);
grid on
kbstr1 = sprintf('Rotation of the disks with k = %.3f, b = %.3f', K, B);
title(kbstr1);
legend('$\theta_1$', '$\theta_2$', 'Location', 'best');
xlabel('Time [s]'); ylabel('Angle of rotation [deg]');
%Plot the difference between the two rotations
figure('Name', 'Exercise 1', 'Windowstate', 'maximized')
plot(time, x_kb(:, 1)*180/pi -x_kb(:, 2)*180/pi, 'linewidth', 1.5);
grid on
title('Difference in rotation ($\theta_1 - \theta_2$)');
xlabel('Time [s]'); ylabel('[deg]');
%Plot angular velocity
figure('Name', 'Exercise 1', 'Windowstate', 'maximized')
plot(time, x_kb(:, 3), time, x_kb(:, 4), 'linewidth', 1.5);
grid on
kbstr2 = sprintf('Angular velocity of the disks with k = %.3f, b = %.3f', K, B);
title(kbstr2);
legend('$\dot{\theta_1}$', '$\dot{\theta_2}$', 'Location', 'Best');
xlabel('Time [s]'); ylabel('Angular velocity [rad/s]');

%LINEAR MODEL --> fsolve does not converge to a minimum

%Try to find the right k and b with the linear model of the damping
% fun = @(var) costfunction_linear(var, J1, J2, T0, time, x0, 1, acc_sample1);
% options = optimset('Algorithm', 'levenberg-marquardt');
% [kb_linear, J_linear] = fsolve(fun, kb0, options); %find the right values of k and b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Exercise 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; 
%%% --- Data --- %%%
%Fluid
data.skydrol.rho = 890;  %kg/m^3
%Accumulator
data.accumulator.V_inf = 10e-3; %m^3
data.accumulator.p_inf = 2.5e6; %Pa
data.accumulator.p0 = 21e6; %Pa
data.accumulator.gamma = 1.2; 
%For isothermal transformation
data.accumulator.V0 = data.accumulator.p_inf*data.accumulator.V_inf/data.accumulator.p0;
%Delivery
data.delivery.kA = 1.12; 
data.delivery.kcv = 2;
data.delivery.D23 = 18e-3; %m
data.delivery.L23 = 2; %m 
data.delivery.f23 = 0.032;
%Area of the delivery line
data.delivery.A23 = (pi*(data.delivery.D23/2)^2); %m^2  
%Distributor
data.distributor.kd = 12; 
data.distributor.d0 = 5e-3; %m
%Actuator
data.actuator.Dc = 50e-3; %m
data.actuator.Dr = 22e-3; %m
data.actuator.m = 2; %kg
data.actuator.xmax = 200e-3; %m
data.actuator.F0 = 1e3; %N 
data.actuator.k = 120e3; %N/m 
data.actuator.Ac = pi*(data.actuator.Dc/2)^2; %m^2 area of the cylinder
data.actuator.As = pi/4*(data.actuator.Dc^2 - data.actuator.Dr^2); %m^2 area of the section
%Return
data.return.Dr = 18e-3; %m
data.return.L67 = 15; %m 
data.return.f67 = 0.035;
data.return.A67 = pi*(data.return.Dr/2)^2; %m^2
%Tank
data.tank.Pt = 0.1e6; %Pa 
data.tank.Vt0 = 1e-3; %m^3
data.tank.kt = 1.12; 
%Time instants and initial conditions
data.t0 = 0;
data.t1 = 1;
data.t2 = 1.5;
data.tf = 3;
data.x0 = 0; 
data.v0 = 0;

%Initial conditions:
q0 = [data.accumulator.V0; data.tank.Vt0; data.x0; data.v0]';
tspan = [data.t0, data.tf];
%First integration to find the pilot event which is when the piston strokes
options = odeset('event', @(t, q)stroke_event(t, q, data));
[tprova, qprova, te, pe, ie] = ode23s(@(t, q) hydraulic(t, q, data), tspan, q0, options);
fprintf('The piston reaches the maximum stroke at t = %f s\n', te);

%%% --- This commented section perform the integration for each interval of time ---%%%
% %Integration from t0 to t1 when command z is still switched off
% tspan_i = [data.t0, data.t1];
% [Ti, Qi] = ode23s(@(t, q) hydraulic(t, q, data), tspan_i, q0);
% %Integration during the ramp of z
% tspan_ramp = [data.t1, data.t2];
% [Tr, Qr] = ode23s(@(t, q) hydraulic(t, q, data), tspan_ramp, q0);
% %Integration from the end of the ramp to the event
% tspan_e = [data.t2, te];
% qi = Qr(end, :);
% [Ts, Qs] = ode23s(@(t, q) hydraulic(t, q, data), tspan_e, qi);
% %Integration from the event to tf
% tspan_f = [te, data.tf];
% qii = [Qs(end, 1:2), data.actuator.xmax, 0]; %impose the velocity equal to zero because it reached the end of the stroke
% [Tf, Qf] = ode23s(@(t, q) hydraulic(t, q, data), tspan_f, qii);
% T = [Ti; Tr; Ts; Tf];
% Q = [Qi; Qr; Qs; Qf];

%Integration from t0 to the time of the event
tspan_tot = [data.t0, te];
[Tinit, Qinit] = ode23s(@(t, q) hydraulic(t, q, data), tspan_tot, q0);
%Integration from the event to tf
tspan_f = [te, data.tf];
qii = [Qinit(end, 1:2), data.actuator.xmax, 0]; %impose the velocity equal to zero because it reached the end of the stroke
[Tf, Qf] = ode23s(@(t, q) hydraulic(t, q, data), tspan_f, qii);

T = [Tinit; Tf];
Q = [Qinit; Qf];

%Plot of the position and velocity of the piston
figure('Name', 'Exercise 2')
colororder({'b','r'})
yyaxis left
p1 = plot(T, Q(:, 3)*10^3, 'LineWidth', 1.5);
ylabel('Position [mm]'); ylim([-5, 390]);
hold on
yyaxis right
p2 = plot(T, Q(:, 4)*10^3, 'LineWidth', 1.5);
ylabel('Velocity [mm/s]');
hold on
plot(te*ones(10,1), linspace(-5, 390, 10), '--k')
grid on
legend([p1, p2], {'Position of the piston', 'Velocity of the piston'}, 'Location', 'northwest');
xlabel('Time [s]'); ylim([-5, 390]);

%Recover the pressures and the flow rates in the system
paraout = zeros(length(T), 10);
for ii = 1: length(T)
    [~, paraout(ii, :)] = hydraulic(T(ii), Q(ii, :), data);
end

Q4 = paraout(:, 1);
Q7 = paraout(:, 2);
PA = paraout(:, 3);
P1 = paraout(:, 4);
P2 = paraout(:, 5);
P3 = paraout(:, 6);
P4 = paraout(:, 7);
P5 = paraout(:, 8);
P6 = paraout(:, 9);
P7 = paraout(:, 10);
%Flow rates
figure('Name', 'Exercise 2')
plot(T, Q4*1e3, 'c', T, Q7*1e3, 'm', 'LineWidth', 1.5); 
xlabel('Time [s]'); ylabel('Flow rate [dm$^3$/s]');
title('Volumetric flow rate in the lines');
grid on
legend('Flow rate in the delivery line: $Q_4$', 'Flow rate in the return line: $Q_5$');
%Plot of all the pressures
figure('Name', 'Exercise 2');
plot(T, paraout(:, 3:end)*10^(-6), 'linewidth', 1.5); %plot of the pressures
grid on
legend('$p_A$', '$p_1$', '$p_2$', '$p_3$', '$p_4$', '$p_5$', '$p_6$', '$p_7$');
xlabel('Time [s]'); ylabel('Pressure [MPa]');
title('Pressures');
%Pressures in the delivery line
figure('Name', 'Exercise 2');
plot(T, PA*1e-6, T, P1*1e-6, T, P2*1e-6, T, P3*1e-6, 'Linewidth', 1.5);
grid on
title('Pressures in the delivery line');
legend('$p_A$', '$p_1$', '$p_2$', '$p_3$');
xlabel('Time [s]'); ylabel('Pressure [MPa]');
ax=axes;
set(ax,'units','normalized','position',[0.2,0.2,0.3,0.23])
box(ax,'on')
plot(T, PA*1e-6, T, P1*1e-6, T, P2*1e-6, T, P3*1e-6, 'LineWidth', 1.5, 'parent', ax)
set(ax,'xlim',[1.44, 1.465],'ylim',[17.76, 17.85], 'FontSize', 10)
%Pressures in the piston
figure('Name', 'Exercise 2');
plot(T, P4*10^(-6), T, P5*10^(-6), 'LineWidth', 1.5)
grid on
legend('$p_4$', '$p_5$');
xlabel('Time [s]'); ylabel('Pressure [MPa]');
title('Pressures in the piston');
%Pressures in the return line
figure('Name', 'Exercise 2');
plot(T, P6*10^(-6), T, P7*10^(-6), 'LineWidth', 1.5)
grid on
legend('$p_6$', '$p_7$');
xlabel('Time [s]'); ylabel('Pressure [MPa]');
title('Pressures in the return line');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----- Exercise 3 -------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; 

R1 = 1000; %om
R2 = 100; %om
L = 1e-3; %H
C = 1e-3; %F
% Point 1 : current generator
x0 = [0; -1/(L*(R2/R1 + 1))]; %initial conditions [i_L0, di_L0]  
% Check on the eigenvalues of the state matrix to assess if the system is
% stiff
A = [ 0                  ,                               1;
     -1/(L*C*(R2/R1 + 1)), -(R2*C + L/R1)/(L*C*(R2/R1 + 1))];
lambda = eig(A);
fprintf('\nExercise 3\n');
fprintf(['The eigenvalues of the state matrix are: \nlambda_1 = %.1f, lambda_2 = %.1f\n'], lambda(1), lambda(2));
fprintf('The system is stiff\n');

[t_c, x_c] = ode23s(@(t, x)electric_circuit(t, x, L, C, R1, R2), [0, 2], x0);

figure('Name', 'Exercise 3', 'Windowstate', 'maximized')
plot(t_c, abs(x_c(:, 1)), 'b', 'LineWidth', 2); %plot of the absolute value of the 
%current passing through the capacitor because the sign only indicates the
%direction of flow. But this depends just on the convention adopted.
title('Time history of the current passing in the inductance $i_L$');
xlabel('Time [s]'); ylabel('$i_L$ [A]'); ylim([-0.002, 0.012]);
grid on

%Time history of Vc
Vc_c = -L*(R2/R1 + 1)*x_c(:, 2) - R2*x_c(:, 1);

figure('Name','Exercise 3', 'WindowState', 'maximized')
plot(t_c, Vc_c, 'r', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('$V_c$ [V]'); 
grid on
title('Time history of the voltage $V_c$')
ylim([-0.2, 1.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Point 3 : Voltage switch
syms t
f = 5; %Hz
v = sin(2*pi*f*t)*atan(t);
dv = matlabFunction(diff(v));
%State matrix --> is the same of point 1 so stiff system
%The initial conditions are the same of point 1 because v(t) is null at t = 0
options =  odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[t_v, x_v] = ode15s(@(t, x)circuit_voltage(t, x, L, C, R1, R2, dv), [0, 12], x0, options);

%Time history of Vc
v = matlabFunction(v);
figure('Name', 'Exercise 3', 'WindowState', 'maximized')
fplot(v, [0, 12], 'LineWidth', 2); 
title('Voltage source');
xlabel('Time [s]'); ylabel('$v(t)$ [V]');
grid on

Vc_v = v(t_v) - (1 + R2/R1)*L*x_v(:, 2) - R2*x_v(:, 1);

figure('Name', 'Exercise 3', 'WindowState', 'maximized')
plot(t_v, x_v(:, 1), 'r', 'LineWidth', 2);
title('Time history of the current passing through the inductance $i_L$');
xlabel('Time [s]'); ylabel('$i_L$ [A]');
grid on

figure('Name', 'Exercise 3', 'WindowState', 'maximized')
plot(t_v, Vc_v, 'b', 'LineWidth', 2)
title('Time history of the voltage on the conductor $V_c$');
grid on
xlabel('Time [s]'); ylabel('$V_c$ [V]');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ Exercise 4 -------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all; clc; 

%Thickness of the layers: 
data.layer.l1 = 2e-3; %[m]
data.layer.l2 = 22e-3; %[m]
data.layer.l4 = 13e-3; %[m]
data.layer.l5 = 3e-3; %[m]
%Dimensional characteristic of the area through which heat transfer occur
data.nozzle.Dext = 2; %external diameter [m]
data.nozzle.thick = data.layer.l1 + data.layer.l2 + data.layer.l4 + data.layer.l5; %9e-3; %total thickness [m]
data.nozzle.L = 3.4; %length of the nozzle [m]
data.nozzle.A = pi*data.nozzle.Dext*data.nozzle.L; %area of heat transfer
%Thermal conductivity
data.cond.k1 = 12; %[W/(mK)] Silicon carbide fiber 
data.cond.k2 = 470; %[W/(mK)] Graphite
data.cond.k4 = 0.5598; %[W/(mK)] Carbon Cloth Phenolic
data.cond.k5 = 40;% Niobium alloy
%Specific heat capacity
data.cap.cp2 = 710; %[J/(kgK)] of the conductor layer
data.cap.cp4 = 1100;%1.42e3; %[J/(kgK)] of the insulator layer
%Density 
data.nozzle.rho2 = 1760; %[kg/m^3]
data.nozzle.rho4 = 1400;%1380; %[kg/m^3]
%Thermal capacitance
data.cap.C2 = data.nozzle.rho2*data.nozzle.A*data.layer.l2*data.cap.cp2; %[J/K]
data.cap.C4 = data.nozzle.rho4*data.nozzle.A*data.layer.l4*data.cap.cp4; %[J/K]
%Data
data.temp.Ti_0 = 20 + 273.15; %[K]
data.temp.To = 20 + 273.15; %[K]
data.temp.Ti_f = 1000 + 273.15; %[K]

data.time.t1 = 1; %[s]
data.time.tf = 60; %[s]

%Resistances --> radial conduction in a plate
data.res.R1 = data.layer.l1/(data.cond.k1*data.nozzle.A);
data.res.R2 = data.layer.l2/(data.cond.k2*data.nozzle.A);
data.res.R3 = 3.5e-5; %data.layer.l3/(data.cond.k3*data.nozzle.A);
data.res.R4 = data.layer.l4/(data.cond.k4*data.nozzle.A);
data.res.R5 = data.layer.l5/(data.cond.k5*data.nozzle.A);
data.res.Req23 = data.res.R2/2 + data.res.R3/2;
data.res.Req34 = data.res.R3/2 + data.res.R4/2;
data.res.Req2_b3 = data.res.R2/3 + data.res.R3/2;
data.res.Req34_a = data.res.R3/2 + data.res.R4/3;

%Compute the response: 
%Eingevalues:
a11 = 2/(data.cap.C2*data.res.R2)*(2/data.res.R2*(2/data.res.R2 + 1/data.res.R1)^(-1) - 1) - 1/(data.cap.C2*data.res.Req23) + 1/(data.cap.C2*data.res.Req23^2)*(1/data.res.Req23 + 1/data.res.Req34)^(-1);
a12 = 1/(data.cap.C2*data.res.Req23*data.res.Req34)*(1/data.res.Req23 + 1/data.res.Req34)^(-1);
a21 = 1/(data.cap.C4*data.res.Req23*data.res.Req34)*(1/data.res.Req23 + 1/data.res.Req34)^(-1);
a22 = 1/(data.cap.C4*data.res.Req34^2)*(1/data.res.Req23 + 1/data.res.Req34)^(-1) - 1/(data.cap.C4*data.res.Req34) - 2/(data.cap.C4*data.res.R4)*(1 - 2/data.res.R4*(2/data.res.R4 + 1/data.res.R5)^(-1));
A = [a11, a12; a21, a22];
lambda = eig(A);
fprintf('Exercise 4\n')
fprintf('The eigenvalues of the system are:\n');
fprintf('lambda_1 = %f\n', lambda(1));
fprintf('lambda_2 = %f\n', lambda(2));
tspan = [0, data.time.tf];

%Since at the initial time both inner and outer layer are at 20°C, every
%point of every section is at 20°C
T0 = [data.temp.Ti_0, data.temp.Ti_0];
[t, T] = ode15s(@(t, x) thermal(t, x, data), tspan, T0);
T2 = T(:, 1);
T4 = T(:, 2);
%Plot over time
%T2 and T4
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(t, T(:, 1)-273.15, 'r', t, T(:, 2)-273.15, 'b', 'Linewidth', 2);
grid on
legend('$T_2$(t)', '$T_4$(t)');
title('Case 1');
xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]');

paraout = zeros(length(t), 3);
for ii = 1: length(t)
    [~, paraout(ii, :)] = thermal(t(ii), T(ii, :), data);
end
T1 = paraout(:, 1);
T3 = paraout(:, 2);
T5 = paraout(:, 3);
%T1, T3, T5
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(t, T1-273.15, 'g', t, T3-273.15, 'm', t, T5-273.15, 'c', 'Linewidth', 2);
grid on
xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]');
title('Case 1');
legend('$T_1$(t)', '$T_3$(t)', '$T_5$(t)');
%All temperatures
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(t, T1-273.15, t, T(:, 1)-273.15, t, T3-273.15, t, T(:, 2)-273.15, t, T5-273.15, 'Linewidth', 1.5);
title('Temperatures'); xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]'); grid on
legend('$T_1$(t)', '$T_2$(t)', '$T_3$(t)', '$T_4$(t)', '$T_5$(t)');
%Plot over space at final time
x = [0, data.layer.l1*1e3, (data.layer.l1 + data.layer.l2/2)*1e3,...
    (data.layer.l1 + data.layer.l2)*1e3, (data.layer.l1 + data.layer.l2 + data.layer.l4/2)*1e3, (data.nozzle.thick - data.layer.l5)*1e3, (data.nozzle.thick)*1e3];
T_space_final = [data.temp.Ti_f, T1(end), T2(end), T3(end), T4(end), T5(end), data.temp.To];
figure('Name', 'Exercise 4', 'Windowstate', 'maximized')
plot(x, T_space_final-273.15, '-or', 'LineWidth', 1.5);
grid on
xlabel('Thickness [mm]'); ylabel('Temperature [$^{\circ}$C]');
title('Evolution of the temperature across the profile at final time');
%Plot vertical dotted lines to delineate thicknesses
X = linspace(0, 1050, 10);

for ii = 1 : length(x)
    hold on
    plot(x(ii)*ones(1, 10), X, '--k');
end
ylim([0, 1050]);
%magari con una colorbar per vedere come è l'andamento nello spazio 

%Plot over space but at different times which are the following:
%[2 sec, 10 sec, 20 sec, 30 sec, 40 sec, 50 sec, 60 sec]
T0 = [data.temp.Ti_0, data.temp.Ti_0];
tspan_sec = [0, 2, 10, 20, 30, 40, 50, 60];
[t_sec, T_sec] = ode15s(@(t, x) thermal(t, x, data), tspan_sec, T0);
T2_sec = T_sec(:, 1)-273.15;
T4_sec = T_sec(:, 2)-273.15;
paraout_sec = zeros(length(t_sec), 3);
for ii = 1: length(t_sec)
    [~, paraout_sec(ii, :)] = thermal(t_sec(ii), T_sec(ii, :), data);
end
T1_sec = paraout_sec(:, 1)-273.15;
T3_sec = paraout_sec(:, 2)-273.15;
T5_sec = paraout_sec(:, 3)-273.15;
Ti_sec = [data.temp.Ti_0-273.15; data.temp.Ti_f*ones(7, 1)-273.15];
T_mat = [Ti_sec, T1_sec, T2_sec, T3_sec, T4_sec, T5_sec, data.temp.To*ones(8,1)-273.15];
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(x, T_mat, '-o', 'linewidth', 1.5);
grid on
xlabel('Thickness [mm]'); ylabel('Temperature [$^{\circ}$C]');
title('Evolution of the temperature across the profile at different times');
%Plot vertical dotted lines to delineate thicknesses
X = linspace(0, 1050, 10);

for ii = 1 : length(x)
    hold on
    plot(x(ii)*ones(1, 10), X, '--k');
end
legend('t = 0 s','t = 2 s', 't = 10 s','t = 20 s','t = 30 s','t = 40 s', ...
    't = 50 s','t = 60 s', '', '', '', '', '', '', '', 'location', 'northeast')
ylim([0, 1050]);

%%%%% ---- Point 3 ----- %%%%%
%Multiple nodes in layer 2 and 4
T0_mn = data.temp.Ti_0*ones(1, 4);
[t_mn, T_mn] = ode15s(@(t, x) thermal_mn(t, x, data), tspan, T0_mn);
T2_a = T_mn(:, 1);
T2_b = T_mn(:, 2);
T4_a = T_mn(:, 3);
T4_b = T_mn(:, 4);
%Plot over time
%T2_a, T2_b, T4_a and T4_b: state variables
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(t_mn, T2_a-273.15, t_mn, T2_b-273.15, t_mn, T4_a-273.15, t_mn, T4_b-273.15, 'Linewidth', 2);
grid on
title('Case 2');
legend('$T_{2a}$(t)', '$T_{2b}$(t)', '$T_{4a}$(t)', '$T_{4b}$(t)', 'location', 'best');
xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]'); ylim([0, 1100]);
paraout_mn = zeros(length(t_mn), 3);
for ii = 1: length(t_mn)
    [~, paraout_mn(ii, :)] = thermal_mn(t_mn(ii), T_mn(ii, :), data);
end
T1_mn = paraout_mn(:, 1);
T3_mn = paraout_mn(:, 2);
T5_mn = paraout_mn(:, 3);
%T1, T3, T5
figure('Name', 'Exercise 4', 'WindowState', 'maximized')
plot(t_mn, T1_mn-273.15, 'g', t_mn, T3_mn-273.15, 'm', t_mn, T5_mn-273.15, 'c', 'Linewidth', 2);
grid on
xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]'); ylim([0, 1100]);
title('Case 2');
legend('$T_1$(t)', '$T_3$(t)', '$T_5$(t)', 'location', 'best');
%All temperatures
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(t_mn, T1_mn-273.15, t_mn, T2_a-273.15, t_mn, T2_b - 273.15, t_mn, T3_mn-273.15, t_mn, T4_a-273.15, t_mn, T4_b -273.15, t_mn, T5_mn-273.15, 'Linewidth', 1.5);
title('Temperatures'); xlabel('Time [s]'); ylabel('Temperature [$^{\circ}$C]'); grid on
legend('$T_1$(t)', '$T_{2a}$(t)', '$T_{2b}$(t)', '$T_3$(t)', '$T_{4a}$(t)', '$T_{4b}$(t)', '$T_5$(t)');
%Plot over space at final time
x_mn = [0, data.layer.l1*1e3, (data.layer.l1 + data.layer.l2/3)*1e3,(data.layer.l1 + 2*data.layer.l2/3)*1e3, ...
    (data.layer.l1 + data.layer.l2)*1e3, (data.layer.l1 + data.layer.l2 + data.layer.l4/3)*1e3, (data.layer.l1 + data.layer.l2 + 2*data.layer.l4/3)*1e3, (data.nozzle.thick - data.layer.l5)*1e3, (data.nozzle.thick)*1e3];
T_space_final_mn = [data.temp.Ti_f, T1_mn(end), T2_a(end),  T2_b(end), T3(end), T4_a(end), T4_b(end), T5(end), data.temp.To];
figure('Name', 'Exercise 4', 'Windowstate', 'maximized')
plot(x_mn, T_space_final_mn-273.15, '-or', 'LineWidth', 1.5);
grid on
xlabel('Thickness [mm]'); ylabel('Temperature [$^{\circ}$C]'); ylim([0, 1050]);
title('Evolution of the temperature across the profile at final times with multiple nodes');
%Plot vertical dotted lines to delineate thicknesses
X = linspace(0, 1050, 10);

for ii = 1 : length(x_mn)
    hold on
    plot(x_mn(ii)*ones(1, 10), X, '--k');
end

%Spatial trend at different time instants
%[2 sec, 10 sec, 20 sec, 30 sec, 40 sec, 50 sec, 60 sec]
T0_mn = [data.temp.Ti_0, data.temp.Ti_0, data.temp.Ti_0, data.temp.Ti_0];
tspan_sec = [0, 2, 10, 20, 30, 40, 50, 60];
[tmn_sec, Tmn_sec] = ode15s(@(t, x) thermal_mn(t, x, data), tspan_sec, T0_mn);
T2amn_sec = Tmn_sec(:, 1)-273.15;
T2bmn_sec = Tmn_sec(:, 1)-273.15;
T4amn_sec = Tmn_sec(:, 3)-273.15;
T4bmn_sec = Tmn_sec(:, 4)-273.15;
paraout_secmn = zeros(length(tmn_sec), 3);
for ii = 1: length(tmn_sec)
    [~, paraout_secmn(ii, :)] = thermal_mn(tmn_sec(ii), Tmn_sec(ii, :), data);
end
T1mn_sec = paraout_secmn(:, 1)-273.15;
T3mn_sec = paraout_secmn(:, 2)-273.15;
T5mn_sec = paraout_secmn(:, 3)-273.15;
Timn_sec = [data.temp.Ti_0-273.15; data.temp.Ti_f*ones(7, 1)-273.15];
Tmn_mat = [Timn_sec, T1mn_sec, T2amn_sec, T2bmn_sec, T3mn_sec, T4amn_sec, T4bmn_sec, T5mn_sec, data.temp.To*ones(8,1)-273.15];
figure('Name', 'Exercise 4', 'windowstate', 'maximized')
plot(x_mn, Tmn_mat, '-o', 'linewidth', 1.5);
grid on
xlabel('Thickness [mm]'); ylabel('Temperature [$^{\circ}$C]');
title('Evolution of the temperature across the profile at different times for multiple nodes');
%Plot vertical dotted lines to delineate thicknesses
X = linspace(0, 1500, 10);

for ii = 1 : length(x_mn)
    hold on
    plot(x_mn(ii)*ones(1, 10), X, '--k');
end

legend('t = 0 s','t = 2 s', 't = 10 s','t = 20 s','t = 30 s','t = 40 s','t = 50 s','t = 60 s', '', '', '', '', '', '', '')
ylim([0, 1050]);
%%
function dx = reactionwheel(~, x, J1, J2, k, b, T)
%This functions aims at providing the dynamics of the reaction wheel
%INPUT:
% x = [th_1 (rad), th_2 (rad), th_dot1 (rad/s), th_dot2 (rad/s)]
% J1 = rotational inertia for first disk
% J2 = rotational inertia for second disk
% k = shaft stiffness
% b = torsional damping coefficient
% T = external torque
%OUTPUT: 
% dx = derivative of x

dx = [x(3);
      x(4);
      -k/J1*x(1) + k/J1*x(2);
      k/J2*x(1) - k/J2*x(2) - sign(x(4))*b/J2*x(4)^2 + T/J2];
end

function dx = reactionwheel_linear(~, x, J1, J2, k, b, T)
%This functions aims at providing the dynamics of the reaction wheel
%considering a linear damper model
%INPUT:
% x = [th_1 (rad), th_2 (rad), th_dot1 (rad/s), th_dot2 (rad/s)]
% J1 = rotational inertia for first disk
% J2 = rotational inertia for second disk
% k = shaft stiffness
% b = torsional damping coefficient
% T = external torque
%OUTPUT: 
% dx = derivative of x

dx = [x(3);
      x(4);
      -k/J1*x(1) + k/J1*x(2);
      k/J2*x(1) - k/J2*x(2) - sign(x(4))*b/J2*x(4) + T/J2];
end

function J = costfunction(var, J1, J2, T0, time, x0, n, acc_sample)
%This function returns the value of the cost function which has to be
%minimized. The cost function is the sum of all the differences between the
%experimental data and the acceleration obtained through numerical
%integration. Depending on which data you want to consider. J changes its
%definition. There are three cases: 
% CASE 1) Only experimental data coming from the first disk
% CASE 2) Only experimental data coming from the second disk
% CASE 3) Experimental data coming from both the disks

%INPUTS: 
% var        : variables of the problem. In this case: [k, b]
% J1         : rotational inertia of the first disk
% J2         : rotational inertia of the second disk
% T0         : torque
% time       : timespan vector [t0 : step : tf]
% x0         : initial condition for the state of the dynamics 
%              x0 = [th0_1 (rad), th0_2 (rad), th0_dot1 (rad/s), th0_dot2 (rad/s)]
% n          : case. It indicates what samples to take into account
% acc_sample : experimental data to provide to compare with the numerical
%              results. In case 1, 2 : [Number of samples x 1]. 
%              In case 3 : [Number of samples x 2].

%OUTPUT:
% J          : cost function

k = var(1);
b = var(2);
[~, x] = ode45(@(t, x)reactionwheel(t, x, J1, J2, k, b, T0), time, x0);

switch n
    case 1
        acc_sample1 = acc_sample; 
        dd_theta1 = k/J1*(x(:, 2) - x(:, 1));
        J = sum((acc_sample1 - dd_theta1).^2);
    case 2
        acc_sample2 = acc_sample;
        dd_theta2 = -k/J2*(x(:, 2) - x(:, 1)) - b/J2*sign(x(:, 4)).*x(:, 4).^2 + T0/J2*ones(size(x(:, 4)));
        J = sum((acc_sample2 - dd_theta2).^2);
    case 3
        acc_sample1 = acc_sample(:, 1);
        acc_sample2 = acc_sample(:, 2);
        dd_theta1 = k/J1*(x(:, 2) - x(:, 1));
        dd_theta2 = -k/J2*(x(:, 2) - x(:, 1)) - b/J2*sign(x(:, 4)).*x(:, 4).^2 + T0/J2*ones(size(x(:, 4)));
        J = sum((acc_sample1 - dd_theta1).^2) + sum((acc_sample2 - dd_theta2).^2);
end

end

function J = costfunction_linear(var, J1, J2, T0, time, x0, n, acc_sample)
%This function returns the value of the cost function which has to be
%minimized. The cost function is the sum of all the differences between the
%experimental data and the acceleration obtained through numerical
%integration. Depending on which data you want to consider. J changes its
%definition. There are three cases: 
% CASE 1) Only experimental data coming from the first disk
% CASE 2) Only experimental data coming from the second disk
% CASE 3) Experimental data coming from both the disks

%INPUTS: 
% var        : variables of the problem. In this case: [k, b]
% J1         : rotational inertia of the first disk
% J2         : rotational inertia of the second disk
% T0         : torque
% time       : timespan vector [t0 : step : tf]
% x0         : initial condition for the state of the dynamics 
%              x0 = [th0_1 (rad), th0_2 (rad), th0_dot1 (rad/s), th0_dot2 (rad/s)]
% n          : case. It indicates what samples to take into account
% acc_sample : experimental data to provide to compare with the numerical
%              results. In case 1, 2 : [Number of samples x 1]. 
%              In case 3 : [Number of samples x 2].

%OUTPUT:
% J          : cost function

k = var(1);
b = var(2);
[~, x] = ode45(@(t, x)reactionwheel_linear(t, x, J1, J2, k, b, T0), time, x0);

switch n
    case 1
        acc_sample1 = acc_sample; 
        dd_theta1 = k/J1*(x(:, 2) - x(:, 1));
        J = sum((acc_sample1 - dd_theta1).^2);
    case 2
        acc_sample2 = acc_sample;
        dd_theta2 = -k/J2*(x(:, 2) - x(:, 1)) - b/J2*sign(x(:, 4)).*x(:, 4).^2 + T0/J2*ones(size(x(:, 4)));
        J = sum((acc_sample2 - dd_theta2).^2);
    case 3
        acc_sample1 = acc_sample(:, 1);
        acc_sample2 = acc_sample(:, 2);
        dd_theta1 = k/J1*(x(:, 2) - x(:, 1));
        dd_theta2 = -k/J2*(x(:, 2) - x(:, 1)) - b/J2*sign(x(:, 4)).*x(:, 4).^2 + T0/J2*ones(size(x(:, 4)));
        J = sum((acc_sample1 - dd_theta1).^2) + sum((acc_sample2 - dd_theta2).^2);
end

end

function dx = electric_circuit(~, x, L, C, R1, R2)
%This function represents the dynamics of the electric circuit

% INPUT:
% x     : state vector [i_L; di_L]
% L     : inductance
% C     : capacitance
% R1    : resistance 1
% R2    : resistance 2

% OUTPUT: 
% dx    : derivative of the state vector [di_L; ddi_L]

dx = [ x(2);
      -(R2*C + L/R1)/(L*C*(R2/R1 + 1))*x(2) - x(1)/(L*C*(R2/R1 + 1))];

end

function dx = circuit_voltage(t, x, L, C, R1, R2, dv)
%This function represents the dynamics of the electric circuit when it is
%present a sinusoidal voltage source

% INPUTS: 
% t   : time
% x   : [i_L; di_L]
% L   : inductance
% C   : capacitance
% R1  : first resistance
% R2  : second resistance
% dv  : derivative of the voltage

% OUTPUT:
% dx  : derivative of the state vector

dx = [x(2);
     -(C*R2 + L/R1)/(C*L*(1 + R2/R1))*x(2) - x(1)/(C*L*(1 + R2/R1)) + C/(C*L*(1 + R2/R1))*dv(t)];
end

function [dq, parout] = hydraulic(t, q, data)
%This function represents the dynamics of the hydraulic system
% INPUTS: 
% t   : time
% q   : state vector [Q_4, Q_5, dot_x, ddot_x]
% data : struct containing all the data of the hydraulic system
%OUTPUTS:
%dq   : derivative of state vector 
%paraout : flowrates and pressure along the lines

%Data definition
Ac = data.actuator.Ac;
As = data.actuator.As;
A23 = data.delivery.A23;
m_act = data.actuator.m;
P0 = data.accumulator.p0;
V0 = data.accumulator.V0;
gamma = data.accumulator.gamma;
rho = data.skydrol.rho;
kA = data.delivery.kA;
kcv = data.delivery.kcv;
kd = data.distributor.kd;
f23 = data.delivery.f23;
L23 = data.delivery.L23;
D23 = data.delivery.D23;
A67 = data.return.A67;
D67 = data.return.Dr; 
L67 = data.return.L67;
f67 = data.return.f67;
F0 = data.actuator.F0;
k = data.actuator.k;
kt = data.tank.kt;
t1 = data.t1;
t2 = data.t2;
PT = data.tank.Pt;
xmax = data.actuator.xmax;
r0 = data.distributor.d0/2;

%Flow rate
Q4 = Ac*q(4);
Q7 = As*q(4);

%Pressures
Pacc = P0*(V0/(V0 + Ac*q(3)))^gamma;
PA = Pacc; %- 0.5*kA*rho/A23^2*Q4*abs(Q4);
P1 = PA - 0.5*kA*rho/A23^2*Q4*abs(Q4);
P2 = P1 - 0.5*kcv*rho/A23^2*Q4*abs(Q4);
P3 = P2 - 0.5*f23*L23/D23*rho/A23^2*Q4*abs(Q4);
P7 = PT + 0.5*kt*rho/A67^2*Q7*abs(Q7);
P6 = P7 + 0.5*f67*L67/D67*rho/A67^2*Q7*abs(Q7);

if q(3)<0  %outside the operative range
    q(3) = 0;
elseif q(3)<=0 && q(4)<0
    q(3) = 0;
    q(4) = 0;
elseif q(3) >= xmax
    q(3) = xmax; 
    q(4) = 0;
elseif q(3) >= xmax && q(4) >0
    q(3) = xmax; 
    q(4) = 0;
end

if t<=t1
    P4 = (F0 + PT*As)/Ac;
    P5 = PT;
    dq = [0; 0; 0; 0];
else
    if t>t1 && t<t2
        alpha = 2*acos(1 - 2*(1 + (t - t2)/(t2 - t1)));
        Ad = 0.5*r0^2*(alpha - sin(alpha));
        if Ad<sqrt(eps)
            P4 = (F0 + PT*As)/Ac;
            P5 = PT;
        else
            P4 = P3 - 0.5*kd*rho/Ad^2*Q4*abs(Q4);
            P5 = P6 + 0.5*kd*rho/Ad^2*Q7*abs(Q7);
        end
    else
        Ad = pi*r0^2;
        P4 = P3 - 0.5*kd*rho/Ad^2*Q4*abs(Q4);
        P5 = P6 + 0.5*kd*rho/Ad^2*Q7*abs(Q7);
    end
dq(1,1) = -Ac*q(4);
dq(2,1) = As*q(4);
dq(3,1) = q(4);
dq(4,1) =  1/m_act*(P4*Ac - P5*As - (F0 + k*q(3)));                   

    if q(3) == xmax && dq(4) >0
        dq(3,1) = 0;
        dq(4,1) = 0;
    end
end

parout = [Q4, Q7, PA, P1, P2, P3, P4, P5, P6, P7];
end

function [value, isterminal, direction] = stroke_event(~,q, data)
%This function is needed to find the time instant when the piston reaches
%the end of the stroke

xmax = data.actuator.xmax;
value = xmax - q(3) ;
isterminal = 1 ;
direction = -1 ;

end

function [dx, paraout] = thermal(t, x, data)
%This function represents the dynamics of the temperatures across the
%surface of the nozzle
% INPUTS: 
% t    :  time
% x    :  state variable [T2, T4]
% data : struct containing all the data about the thermal system
% OUTPUTS:
% dx   : derivative of the state vector
%paraout : temperatures along the profile

T2 = x(1);
T4 = x(2);
%Check on the time to take into account the ramp of Ti
if t>= 0 && t<= data.time.t1
    Ti = (data.temp.Ti_f - data.temp.Ti_0)*(t-data.time.t1)/data.time.t1 + data.temp.Ti_f;
else 
    Ti = data.temp.Ti_f;
end

T1 = (Ti/data.res.R1 + 2*T2/data.res.R2)/(2/data.res.R2 + 1/data.res.R1);
T3 = (T2/data.res.Req23 + T4/data.res.Req34)/(1/data.res.Req23 + 1/data.res.Req34);
T5 = (2*T4/data.res.R4 + data.temp.To/data.res.R5)/(2/data.res.R4 + 1/data.res.R5);

T2_dot = 1/data.cap.C2*(2/data.res.R2*(T1 - T2) - 1/data.res.Req23*(T2 - T3));
T4_dot = 1/data.cap.C4*(1/data.res.Req34*(T3 - T4) - 2/data.res.R4*(T4 - T5));

dx = [T2_dot;
      T4_dot];

paraout = [T1, T3, T5];
end

function [dx, paraout] = thermal_mn(t, x, data)
%This function represents the dynamics of the temperatures across the
%surface of the nozzle
%Used for the point 3 of the exercise 4: conductor and insulator have two
%nodes
% INPUTS: 
% t    :  time
% x    :  state variable [T2, T4]
% data : struct containing all the data about the thermal system
% OUTPUTS:
% dx   : derivative of the state vector
%paraout : temperatures along the profile

T2_a = x(1);
T2_b = x(2);
T4_a = x(3);
T4_b = x(4);

if t>= 0 && t<= data.time.t1
    Ti = (data.temp.Ti_f - data.temp.Ti_0)*(t-data.time.t1)/data.time.t1 + data.temp.Ti_f;
else 
    Ti = data.temp.Ti_f;
end

T1 = (Ti/data.res.R1 + 3*T2_a/data.res.R2)/(3/data.res.R2 + 1/data.res.R1);
T3 = (T2_b/data.res.Req2_b3 + T4_a/data.res.Req34_a)/(1/data.res.Req2_b3 + 1/data.res.Req34_a);
T5 = (3*T4_b/data.res.R4 + data.temp.To/data.res.R5)/(3/data.res.R4 + 1/data.res.R5);

T2_adot = 2/data.cap.C2*(3/data.res.R2*(T1 - T2_a) - 3/data.res.R2*(T2_a - T2_b));
T2_bdot = 2/data.cap.C2*(3/data.res.R2*(T2_a - T2_b) - 1/data.res.Req2_b3*(T2_b - T3));
T4_adot = 2/data.cap.C4*(1/data.res.Req34_a*(T3 - T4_a) - 3/data.res.R4*(T4_a - T4_b));
T4_bdot = 2/data.cap.C4*(3/data.res.R4*(T4_a - T4_b) - 3/data.res.R4*(T4_b - T5));

dx = [T2_adot;
      T2_bdot;
      T4_adot;
      T4_bdot];

paraout = [T1, T3, T5];
end

