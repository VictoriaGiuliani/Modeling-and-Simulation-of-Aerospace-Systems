%%%% MODELING AND SIMULATION OF AEROSPACE SYSTEMS (2021/2022) %%%%
% Assignment 1 
% Author: Victoria Katia Giuliani
%%
%For plots

set(groot,'defaulttextinterpreter','Latex');
set(groot,'defaultAxesTickLabelInterpreter','Latex');
set(groot,'defaultLegendInterpreter','Latex');
set(0,'defaultAxesFontSize',20);

%%
%%% EXERCISE 1 %%%

clearvars; close all; clc; 
format long 
%Data
f = @(x) cos(x) - x; 
x = linspace(-10, 10, 50);
a = -2; 
b = 8;
%Plot of the function
figure('Name', 'Exercise 1');
plot(a, f(a), 'or', 'MarkerFaceColor', 'r', 'Markersize', 8);
hold on
plot(b, f(b), 'ob', 'MarkerFaceColor', 'b', 'Markersize', 8);
hold on
plot(x, f(x), 'k'); grid on; 
title('$f = cos(x) - x$'); xlabel('$x$'); ylabel('$f(x)$');
legend('$f(a)$', '$f(b)$');

nmax = 35;
tol = 1e-8;

%Bisection method
tic
[x_b, n_iter_b, fevals_b] = bisection(a, b, f, nmax, tol);
time_b = toc;

%Secant method
tic
[x_s, n_iter_s, fevals_s] = secant(a, b, f, tol, nmax);
time_s = toc; 

%Regula falsi method
tic
[x_rf, n_iter_rf, fevals_rf] = regula_falsi(a, b, f, tol, nmax);
time_rf = toc;

% Print results;
fprintf('Exercise 1:\n');
methods = {'Bisection'; 'Secant'; 'Regula Falsi'};
Solution = [x_b; x_s; x_rf];
Iterations = [n_iter_b; n_iter_s; n_iter_rf];
FunctionEvaluations = [fevals_b; fevals_s; fevals_rf];
Time = [time_b; time_s; time_rf];
Results = table(Solution, Iterations, FunctionEvaluations, Time, 'RowNames', methods)

%Take a fixed and let b vary 
B = 1 : 0.5: 10;
N_b = zeros(length(B), 1);   T_b = N_b;
N_s = zeros(length(B), 1);   T_s = N_s; 
N_rf = zeros(length(B), 1);  T_rf = N_rf;

for ii = 1 : length(B)
    tic
    [~, ~, N_b(ii)] = bisection(a, B(ii), f, nmax, tol);
    T_b(ii) = toc;
    tic
    [~, ~, N_s(ii)] = secant(a, B(ii), f, tol, nmax);
    T_s(ii) = toc; 
    tic
    [~, ~, N_rf(ii)] = regula_falsi(a, B(ii), f, tol, nmax);
    T_rf(ii) = toc;
end
%Plot of function evaluations as function of the distance between 
%the extremes a and b
figure('Name', 'Exercise 2')
scatter(abs(a*ones(length(B),1) - B'), N_b, 60, 'r', 'filled');
hold on
scatter(abs(a*ones(length(B),1) - B'), N_s, 60, 'b', 'filled');
hold on 
scatter(abs(a*ones(length(B),1) - B'), N_rf, 60, 'g', 'filled');
grid on
title('Function evaluations as function of I = $|b - a|$ with fixed a = $-2$')
legend({'Bisection', 'Secant', 'Regula  Falsi'}, 'FontSize', 20, 'Location', 'best')
xlabel('$|b - a|$'); ylabel('Function evaluations');
%Plot of the computational time required by the three methods as function 
%of the number of function evaluations
figure('Name', 'Exercise 1')
scatter(N_b, T_b, 60, 'r', 'filled');
hold on
scatter(N_s, T_s, 60, 'b', 'filled');
hold on
scatter(N_rf, T_rf, 60, 'g', 'filled');
grid on
title('Computational time');
xlabel('Number of function evaluations'); ylabel('Time');
legend('Bisection', 'Secant', 'Regula Falsi', 'Location', 'best')
%%
%%% EXERCISE 2 %%%
clearvars; close all; clc; 
format long
syms 'x1' 'x2'

f = [x1^2 - x1 - x2; 
    x1^2/16 + x2^2 - 1];
f2 = 0*x1 + 0*x2;
%Plot 3D of the vector valued function
figure('Name', 'Exercise 1');
fs1 = fsurf(f(1));
hold on
fs2 = fsurf(f(2));
hold on
fs3 = fsurf(f2);
hold on
scatter3(1.58101, 0.918573, 0, 200, 'r', 'filled'); %first root
hold on
scatter3(-0.612743, 0.988197, 0, 200, 'r', 'filled'); %second root
fs1.FaceColor = 'c'; fs2.FaceColor = 'y'; fs3.FaceColor = 'b';
xlabel('$x_1$'); ylabel('$x_2$'); zlabel('$f(x_1, x_2)$');
%From the plot it can be noticed that there exist two zeros

x_var = [x1; x2];
%Initial conditions for the two zeros
x0 = [18 -10;    % x0 = [x1_initial_z1, x1_initial_z2;
      15 14];    %       x2_initial_z1, x2_initial_z2]
nmax = 30;
tol = 1e-8;

X_an = zeros(2, 2); X_fd = X_an; X_bd = X_an; X_cd = X_an;
N = zeros(2, 4);
ERR_an = zeros(2, 2); ERR_fd = ERR_an; ERR_bd = ERR_an; ERR_cd = ERR_an;
for ii = 1 : 2
    %Zeros found by using Newton's method with the derivatives computed
    %analytically
    [x_an, n_iter_an, err_an] = newton_analytical(f, x_var, x0(:, ii), nmax, tol);

    %Zeros found by using Netwon's method with the derivatives computed with
    %forward differences
    [x_fd, n_iter_fd, err_fd] = newton_fd(f, x_var, x0(:, ii), nmax, tol);

    %Zeros found by using Netwon's method with the derivatives computed with
    %backward differences
    [x_bd, n_iter_bd, err_bd] = newton_bd(f, x_var, x0(:, ii), nmax, tol);

    %Zeros found by using Netwon's method with the derivatives computed with
    %centered differences
    [x_cd, n_iter_cd, err_cd] = newton_cd(f, x_var, x0(:, ii), nmax, tol);
    
    X_an(:, ii) = x_an; X_fd(:, ii) = x_fd; X_bd(:, ii) = x_bd; X_cd(:, ii) = x_cd;
    N(ii, :) = [n_iter_an, n_iter_fd, n_iter_bd, n_iter_cd];
    ERR_an(:, ii) = err_an; ERR_fd(:, ii) = err_fd; ERR_bd(:, ii) = err_bd; ERR_cd(:, ii) = err_cd;
end

%Print results
fprintf('Exercise 2:\n');
fprintf('There exist two zeros of the function f\n\n');
fprintf('FIRST ZERO:\n');
fprintf('   Analytical Solution: x1 = %.14f, x2 = %.14f\n', X_an(1, 1), X_an(2, 1));
fprintf('   Error = %.8e\n\n', norm(ERR_an(:, 1)));
fprintf('   Solution with Forward Difference: x1 = %.14f, x2 = %.14f\n', X_fd(1, 1), X_fd(2, 1));
fprintf('   Error = %.8e\n\n', norm(ERR_fd(:, 1)));
fprintf('   Solution with Backward Difference: x1 = %.14f, x2 = %.14f\n', X_bd(1, 1), X_bd(2, 1)); 
fprintf('   Error = %.8e\n\n', norm(ERR_bd(:, 1)));
fprintf('   Solution with Centered Difference: x1 = %.14f, x2 = %.14f\n', X_cd(1, 1), X_cd(2, 1)); 
fprintf('   Error = %.8e\n\n', norm(ERR_cd(:, 1)));
fprintf('SECOND ZERO:\n');
fprintf('   Analytical Solution: x1 = %.8f, x2 = %.8f\n', X_an(1, 2), X_an(2, 2));
fprintf('   Error = %.8e\n\n', norm(ERR_an(:, 2)));
fprintf('   Solution with Forward Difference: x1 = %.8f, x2 = %.8f\n', X_fd(1, 2), X_fd(2, 2)); 
fprintf('   Error = %.8e\n\n', norm(ERR_fd(:, 2)));
fprintf('   Solution with Backward Difference: x1 = %.8f, x2 = %.8f\n', X_bd(1, 2), X_bd(2, 2)); 
fprintf('   Error = %.8e\n\n', norm(ERR_bd(:, 2)));
fprintf('   Solution with Centered Difference: x1 = %.8f, x2 = %.8f\n', X_cd(1, 2), X_cd(2, 2)); 
fprintf('   Error = %.8e\n\n', norm(ERR_cd(:, 2)));


%%
%%% EXERCISE 3 %%%
clearvars; close all; clc; 

f = @(x, t) x - t^2 + 1;
x0 = 0.5;

t_i = 0; t_f = 2; 
h = [0.5, 0.2, 0.05, 0.01];

%Cell structure where the first row is the solution x(t)_Heun 
%with a specific step size h, the second row reports the time instants
%of the integration ([t_i: h : t_f]). Every column correspond to a different
%step size
X_H = cell([2, length(h)]);


%Solution of the IVP with Heun's method

for ii = 1: length(h)
    [x, t] = heun(f, h(ii), x0, t_i, t_f);
    X_H(1, ii) = {x};
    X_H(2, ii) = {t}; %save the relative time vector
end

%Computation of the CPU time
time_h = zeros(length(h), 1);
for ii = 1 : 100
    [x, ~] = heun(f, h(1), x0, t_i, t_f);

    for jj = 1 : length(h)
        tic
        [x, ~] = heun(f, h(jj), x0, t_i, t_f);
        time_h(jj) = toc;
    end
end
 
%Analytical solution
tt = linspace(t_i, t_f, 100);
f_an = @(t) t.^2 + 2*t + 1 - 0.5*exp(t); 

%Comparison between Heun's solution and analytical solution
figure('Name', 'Exercise 3', 'WindowState', 'maximized');
subplot(2, 2, 1)
plot(tt, f_an(tt), 'k', X_H{2, 1}, X_H{1, 1}, 'r', 'LineWidth', 3);
title(['h = ', num2str(h(1))]);
grid on
xlabel('t'); ylabel('x(t)');
subplot(2, 2, 2)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_H{2, 2}, X_H{1, 2}, 'r', 'LineWidth', 3);
title(['h = ', num2str(h(2))]);
grid on
xlabel('t'); ylabel('x(t)');
subplot(2, 2, 3)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_H{2, 3}, X_H{1, 3}, 'r', 'LineWidth', 3);
grid on
title(['h = ', num2str(h(3))]);
xlabel('t'); ylabel('x(t)');
subplot(2, 2, 4)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_H{2, 4}, X_H{1, 4}, 'r', 'LineWidth', 3);
grid on
legend('Analytical Solution', 'RK2 Solution', 'Location', 'NorthWest');
title(['h = ', num2str(h(4))]);
xlabel('t'); ylabel('x(t)');
sgtitle('Heun Method', 'FontSize', 25);

%Integration error at final time (which corresponds to the maximum)
ERR_H = zeros(length(h), 1);
for ii = 1 : length(h)
    ERR_H(ii) = max(abs(X_H{1, ii} - f_an(X_H{2, ii})));
end

fprintf('Exercise 3\n');
fprintf('HEUN METHOD');
H = {'h = 0.5', 'h = 0.2', 'h = 0.05', 'h = 0.01'};
CPU_time = time_h;
IntegrationError = ERR_H;
results = table(CPU_time, IntegrationError, 'RowNames', H)

%Solution of the IVP with RK4

%Cell structure where the first row is the solution x(t)_RK4 
%with a specific step size h, the second row reports the time instants
%of the integration ([t_i: h : t_f]). Every column correspond to a different
%step size
X_RK4 = cell([2, length(h)]);

for ii = 1: length(h)
    [x, t] = RK4(f, h(ii), x0, t_i, t_f);
    X_RK4(1, ii) = {x};
    X_RK4(2, ii) = {t};
end


%Computation of the CPU time
time_RK4 = zeros(length(h), 1);
for ii = 1 : 100
    [x, ~] = RK4(f, h(1), x0, t_i, t_f);

    for jj = 1 : length(h)
        tic
        [x, ~] = RK4(f, h(jj), x0, t_i, t_f);
        time_RK4(jj) = toc;
    end
end  

%Comparison between RK4 solution and analytical solution
figure('Name', 'Exercise 3', 'WindowState', 'maximized');
subplot(2, 2, 1)
plot(tt, f_an(tt), 'k', X_RK4{2, 1}, X_RK4{1, 1}, 'g', 'LineWidth', 3);
title(['h = ', num2str(h(1))]);
xlabel('t'); ylabel('x(t)');
grid on
subplot(2, 2, 2)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_RK4{2, 2}, X_RK4{1, 2}, 'g', 'LineWidth', 3);
title(['h = ', num2str(h(2))]);
xlabel('t'); ylabel('x(t)');
grid on
subplot(2, 2, 3)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_RK4{2, 3}, X_RK4{1, 3}, 'g', 'LineWidth', 3);
title(['h = ', num2str(h(3))]);
xlabel('t'); ylabel('x(t)');
grid on
subplot(2, 2, 4)
%figure('Name', 'Exercise 3');
plot(tt, f_an(tt), 'k', X_RK4{2, 4}, X_RK4{1, 4}, 'g', 'LineWidth', 3);
title(['h = ', num2str(h(4))]);
xlabel('t'); ylabel('x(t)');
grid on
legend('Analytical Solution', 'RK2 Solution', 'Location', 'NorthWest');
sgtitle('Runge Kutta 4', 'FontSize', 25);
%Integration error at final time (which corresponds to the maximum)
ERR_RK4 = zeros(length(h), 1);
for ii = 1 : length(h)
    ERR_RK4(ii) = max(abs(X_RK4{1, ii} - f_an(X_RK4{2, ii})));
end

fprintf('RK4 METHOD')
CPU_time = time_RK4;
IntegrationError = ERR_RK4;
results = table(CPU_time, IntegrationError, 'RowNames', H)
%%
%%% EXERCISE 4 %%%
clearvars; close all; clc;
syms alpha 

%Matrix A
A = [0,  1;
     -1,  2*cos(alpha)];
L = eig(A);

%Definition of the vector containing all the alphas
alpha_v = pi: -pi/180 : 0;

%Computation of the eiganvalues of A
L_v = eval(subs(L, alpha, alpha_v));
L_v_real = [real(L_v(1, :)), fliplr(real(L_v(2, :)))];
L_v_imag = [imag(L_v(1, :)), fliplr(imag(L_v(2, :)))];

%Plot of the locus of the eigenvalues of A
figure('Name', 'Exercise 4')
plot(L_v_real, L_v_imag, 'b', 'LineWidth', 2);
xlabel('$\Re(\lambda)$');
ylabel('$\Im(\lambda)$');
title('Eigenvalues of A as function of $\alpha$');
axis equal
grid on

% Point 2: with alpha = pi, find the h such that max(|eig(F_RK2(h, alpha)|) = 1
A_pi = eval(subs(A, alpha, pi));               %Evaluate A in alpha = pi
F_RK2 = @(h) eye(2) + h*A_pi + 0.5*h^2*A_pi^2; %Forward operator F_RK2(alpha = pi)
fcn = @(h) max(abs(eig(F_RK2(h)))) - 1;  %Zero finding problem
%Plot fcn as function of h to choose a valid initial guess for h
h = 0: 0.1: pi;
figure('Name', 'Exercise 4')
FCN = zeros(length(h), 1);
for ii = 1 : length(h)
    FCN(ii) = fcn(h(ii));
end
plot(h, FCN, 'k', 'LineWidth', 2);
hold on
plot(h, 0*h, '--g');
hold on
annotation('ellipse', [.56 .24 .15 .1], 'Color', 'r', 'LineWidth', 2);
annotation('textarrow', [.5 .6], [.45 .35], 'String', ...
   'Region of interest for initial guess', 'Color', 'r', 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
grid on
title('$f(x) = max(|eig(F_{RK2}(h)|) - 1$');
xlabel('$h$'); ylabel('$f(x)$');
clear FCN
%From the plot a proper initial guess is:
h0 = 2;
%h solution of max(|eig(F_RK2(h, pi)|) = 1
h_pi = fzero(fcn, h0);  

%Plot of the mechanism of mapping done by RK2
figure('Name', 'Exercise 4')
x_lambda = real(eig(A_pi)); y_lambda = imag(eig(A_pi));
scatter(real(h_pi*eig(A_pi)), imag(h_pi*eig(A_pi)), 80, 'r', 'filled');
hold on
scatter(real(eig(A_pi)), imag(eig(A_pi)), 80, 'b', 'filled');
x_arr = x_lambda(1); y_arr = y_lambda(1);
dx_arr = (h_pi*x_lambda(1)) - x_lambda(1);
dy_arr = h_pi*y_lambda(1) - y_lambda(1);
quiver(x_arr, y_arr, dx_arr, dy_arr, 0, 'k', 'MaxHeadSize', 0.1);
text(x_lambda(1) - 0.15, 0.02, '$\lambda = -1$', 'Color', 'b', 'FontSize', 20);
text(-1.60, 0.02, 'h = 2', 'FontSize', 20);
text(-2.09, 0.02, '$h\lambda = -2$', 'Color', 'r', 'FontSize', 20);
grid on
xlim([-2.1, -0.9]); ylim([-0.25, 0.25]);
title('Mapping of RK2 with $\alpha = \pi$');

% Point 3 : Stability region of Runge-Kutta 2
A = @(alpha) [0, 1; -1, 2*cos(alpha)];
F_RK2 = @(h, alpha) eye(2) + h*A(alpha) + 0.5*h^2*A(alpha)^2; %F_RK2

%Find the points to plot in the {h lambda}-plane
[~, LAMBDA_RK2] = stability(alpha_v, F_RK2, h0, L_v);

figure('Name', 'Exercise 4')
plot(real(LAMBDA_RK2), imag(LAMBDA_RK2), 'k', real(LAMBDA_RK2), -imag(LAMBDA_RK2), 'k', 'LineWidth', 3);
hold on
fill([real(LAMBDA_RK2), real(LAMBDA_RK2)], [imag(LAMBDA_RK2), -imag(LAMBDA_RK2)], 'red', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
hold on
plot(zeros(10, 1), linspace(-3, 3, 10), '--k', 'LineWidth', 2);
set(gca, 'Color', rgb('PowderBlue'));
text(-1.5, 0.75, 'Stable', 'FontSize', 20);
text(2, -2.3, 'Unstable', 'FontSize', 20);
title('Stability region for RK2');
xlim([-3, 3]); ylim([-3, 3]);
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
grid on
axis equal

% Point 4: with alpha = pi, find the h such that max(|eig(F_RK4(h, alpha)|) = 1

%Forward operator F_RK4(alpha = pi)
F_RK4 = @(h) eye(2) + h*A_pi + 0.5*h^2*A_pi^2 + 1/6*h^3*A_pi^3 + 1/24*h^4*A_pi^4;
fcn = @(h) max(abs(eig(F_RK4(h)))) - 1; %Zero finding problem
%plot fcn as function of h to choose a valid initial guess 
h = 0 : 0.1: 4;
figure('Name', 'Exercise 4')
FCN = zeros(length(h), 1);
for ii = 1 : length(h)
    FCN(ii) = fcn(h(ii));
end
plot(h, FCN, 'k', 'LineWidth', 2);
hold on
plot(h, 0*h, '--g');
hold on
annotation('ellipse', [.6 .22 .15 .1], 'Color', 'r', 'LineWidth', 2);
annotation('textarrow', [.53 .63], [.45 .32], 'String', ...
   'Region of interest for initial guess', 'Color', 'r', 'Interpreter', 'latex', 'HorizontalAlignment', 'left');
grid on
title('$f(x) = max(|eig(F_{RK4}(h)|) - 1$');
xlabel('$h$'); ylabel('$f(x)$');

%From the plot, a proper initial guess for h is:
h0 = 3;
%h solution of max(|eig(F_RK2(h, pi)|) = 1
h_RK4 = fzero(fcn, h0);

%Plot of the mechanism of mapping done by RK4
figure('Name', 'Exercise 4')
scatter(real(h_RK4*eig(A_pi)), imag(h_RK4*eig(A_pi)), 80, 'r', 'filled');
hold on
scatter(real(eig(A_pi)), imag(eig(A_pi)), 80, 'b', 'filled');
dx_arr = h_RK4*x_lambda(1) - x_lambda(1);
dy_arr = h_RK4*y_lambda(1) - y_lambda(1);
quiver(x_arr, y_arr, dx_arr, dy_arr, 0, 'k', 'MaxHeadSize', 0.1);
text(x_lambda(1) - 0.2, 0.02, '$\lambda = -1$', 'Color', 'b', 'FontSize', 14);
text(-2.2, 0.02, '$h = 2.78$ ', 'FontSize', 14);
text(h_RK4*x_lambda(1) - 0.2, 0.02, '$h\lambda = -2.78$', 'Color', 'r', 'FontSize', 14);
grid on
xlim([-3, -0.9]); ylim([-0.25, 0.25]);
title('Mapping of RK4 with $\alpha = \pi$');

%Plot stability region of Runge Kutta 4
A = @(alpha) [0, 1; -1, 2*cos(alpha)];
F_RK4 = @(h, alpha) eye(2) + h*A(alpha) + 0.5*h^2*A(alpha)^2 + ...
        1/6*h^3*A(alpha)^3 + 1/24*h^4*A(alpha)^4;

[~, LAMBDA_RK4] = stability_fzero(alpha_v, F_RK4, h0, L_v);
save('Lambda_RK4', 'LAMBDA_RK4');
figure('Name', 'Exercise 4')
plot(real(LAMBDA_RK4), imag(LAMBDA_RK4), 'k', real(LAMBDA_RK4), -imag(LAMBDA_RK4), 'k', 'LineWidth', 3);
hold on
fill([real(LAMBDA_RK4), real(LAMBDA_RK4)], [imag(LAMBDA_RK4), -imag(LAMBDA_RK4)], 'r', ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
set(gca, 'Color', rgb('PowderBlue'));
text(-1.9, 1, 'Stable', 'FontSize', 20);
text(2, -1.5, 'Unstable', 'FontSize', 20);
title('Stability region for RK4');
xlim([-3, 3]); ylim([-3, 3]);
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
grid on
axis equal

%Points hi*lambda_ex3 with t = 0
lambda_ex3 = 1; %eigenvalue of Exercise 3
h_ex3 = [0.5, 0.2, 0.05, 0.01]; %values of h from Exercise 3
hL_ex3 = lambda_ex3*h_ex3;
%Plot of the points in RK2 stability region
figure(4)
scatter(real(hL_ex3), imag(hL_ex3), 'kx', 'LineWidth', 5);
xlim([-3, 5]);
ax=axes;
set(ax,'units','normalized','position',[0.5,0.65,0.25,0.2])
box(ax,'on')
scatter(real(hL_ex3), imag(hL_ex3), 'kx',  'LineWidth', 5, 'parent', ax)
hold on
plot(real(LAMBDA_RK2), imag(LAMBDA_RK2), 'k', real(LAMBDA_RK2), -imag(LAMBDA_RK2), 'k', 'LineWidth', 3, 'parent', ax);
hold on
fill([real(LAMBDA_RK2), real(LAMBDA_RK2)], [imag(LAMBDA_RK2), -imag(LAMBDA_RK2)], 'red', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
set(ax,'xlim',[-0.05,0.55],'ylim',[-0.1, 0.1], 'FontSize', 10)
set(gca, 'Color', rgb('PowderBlue'));
grid on
%Plot of the points in RK4 stability region
figure(7)
scatter(real(hL_ex3), imag(hL_ex3), 'kx', 'LineWidth', 5);
xlim([-3.5, 5]);
ax=axes;
set(ax,'units','normalized','position',[0.52,0.55,0.25,0.2])
box(ax,'on')
scatter(real(hL_ex3), imag(hL_ex3), 'kx',  'LineWidth', 5, 'parent', ax)
hold on
plot(real(LAMBDA_RK4), imag(LAMBDA_RK4), 'k', real(LAMBDA_RK4), -imag(LAMBDA_RK4), 'k', 'LineWidth', 3, 'parent', ax);
hold on
fill([real(LAMBDA_RK4), real(LAMBDA_RK4)], [imag(LAMBDA_RK4), -imag(LAMBDA_RK4)], 'red', 'FaceAlpha', 0.8, 'EdgeColor', 'none'); 
set(ax,'xlim',[-0.05,0.55],'ylim',[-0.1, 0.1], 'FontSize', 10)
set(gca, 'Color', rgb('PowderBlue'));
grid on

%%
%%% EXERCISE 5 %%%
clear; close all; clc;
% % Definition of the matrix A
syms alpha
A = [0,            1; 
    -1, 2*cos(alpha)];
% Eigenvalues of A
L = eig(A);

% Vector of alpha values
alpha_v = pi: -pi/180 : 0; 

% Eigenvalues of A evaluated at each alpha
L_v = eval(subs(L, alpha, alpha_v));  
L_vreal = real(L_v(1, :))';
L_vimag = imag(L_v(1, :))';

% figure()
% plot(real(L_v(1, :)), imag(L_v(1, :)));
% axis equal

% Data
x0 = [1, 1]';
t0 = 0; tf = 1; 
tol = [1e-3, 1e-4, 1e-5, 1e-6];

% ------------------------------------------------------------------------ % 
% RK1

% Plot of the function ||x_an(1) - x_RK1(1)|| = tol for alpha = pi
% to choose a proper initial value for h

% Global view of the function for tol = 1e-3 and for alpha = pi
H_g = linspace(eps, 1, 100);
FCN_g = zeros(length(H_g), 1);
A_pi = eval(subs(A, alpha, pi));
x_an = expm(A_pi*tf)*x0;
% Operator linking the solution at tf = 1 with the initial condition x0
x_RK1 = @(h) (eye(2) + (1 - floor(1/h)*h)*A_pi)*(eye(2) + h*A_pi)^(floor(1/h))*x0;
fcn = @(h) max(abs(x_an - x_RK1(h))) - tol(1); 
for ii = 1 : length(H_g)
    FCN_g(ii) = fcn(H_g(ii));
end
figure('Name', 'Exercise 5')
plot(H_g, FCN_g, 'LineWidth', 2);
hold on
plot(H_g, 0*ones(length(H_g)), '--r', 'LineWidth', 2); 
grid on
legend('f(h)', 'Zero line', 'Location', 'best');
xlabel('$h$'); ylabel('$max(|x_{an(1)} - x_{RK1(1)}|) - tol$');
title('$tol = 10^{-3}, \alpha = \pi$');

% Local view of the same function to focus better on the zero 
H = linspace(0.5*1e-3, 5*1e-3, 5);
FCN = zeros(length(H), 1);

for ii = 1 : length(H)
    FCN(ii) = fcn(H(ii));
end
  
figure('Name', 'Exercise 5')
plot(H, FCN, 'LineWidth', 2);
hold on
plot(H, 0*ones(length(H)), 'r', 'LineWidth', 2); 
grid on
legend('f(h)', 'Zero line', 'Location', 'best');
xlabel('$h$'); xlim([0.5*1e-3, 5*1e-3]); ylabel('$max(|x_{an(1)} - x_{RK1(1)}|) - tol$');
title('Closer look at the zero for $tol = 10^{-3}, \alpha = \pi$');

% Search for h solution of the problem

%matrix containing values of h found for each tolerance and each value of alpha
%where every column corresponds to a different tolerance, each row corresponds
%to a different value of alpha
H = zeros(length(alpha_v), length(tol)); 
%RK1
for ii = 1 : length(tol)
    h0 = 10^(-2-ii); %initial condition
    for jj = 1 : length(alpha_v)
        A_v = eval(subs(A, alpha, alpha_v(jj)));
        x_an = expm(A_v*tf)*x0;
        x_RK1 = @(h) (eye(2) + (1 - floor(1/h)*h)*A_v)*(eye(2) + h*A_v)^(floor(1/h))*x0;
        fcn = @(h) max(abs(x_an - x_RK1(h))) - tol(ii);                    
        H(jj, ii) = fzero(fcn, [h0/10, h0*10]);       
        h0 = H(jj, ii);
    end
end

% Plot of the locus of solution in the {hlambda} plane
hlambda_r = H.*L_vreal;
hlambda_i = H.*L_vimag;

figure('Name', 'Exercise 5')
p1 = plot(hlambda_r(:, 1), hlambda_i(:, 1), 'k'); hold on; plot(hlambda_r(:, 1), -hlambda_i(:, 1), 'k');
hold on
p2 = plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b'); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b');
hold on
p3 = plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r'); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r');
hold on
p4 = plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g'); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g');
grid on
axis equal
title('$||x_{an}(1) - x_{RK1}(1)||_\infty = tol$');
legend([p1 p2 p3 p4], '$tol = 10^{-3}$', '$tol = 10^{-4}$', '$tol = 10^{-5}$', '$tol = 10^{-6}$', 'Location', 'best');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% ----- SINGLE PLOTS ----- %
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b', 'LineWidth', 2); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK1}(1)||_\infty = 10^{-4}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r', 'LineWidth', 2); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK1}(1)||_\infty = 10^{-4}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g', 'LineWidth', 2); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK1}(1)||_\infty = 10^{-4}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% Function evaluations vs tol for alpha = pi
fcn_eval_RK1 = ceil((tf - t0)./H(1, :)); %1 fcn evaluation for step
figure('Name', 'Exercise 5')
loglog(tol, fcn_eval_RK1, '-or', 'LineWidth', 2);
grid on
title('Function evaluations vs tol for $\alpha = \pi$');
axis fill
xlabel('Tolerance'); ylabel('Number of fcn evaluations')

% Plot of the analytical solution for some value of alpha as function of time 
% to understand how its behaviour changes
A_31 = eval(subs(A, alpha, 31*pi/180));  %alpha = 31°
A_0 = eval(subs(A, alpha, 0*pi/180));    %alpha = 0°
time = 0 : 0.1 : 1;
x_an31 =@(T) expm(A_31*T)*x0; %analytical solution for alpha = 31°
x_an0 =@(T) expm(A_0*T)*x0;   %analytical solution for alpha = 0°
X31 = zeros(2, length(time));
X0 = zeros(2, length(time));
for t = 1 : length(time)
    X31(:, t) = x_an31(time(t));
    X0(:, t) = x_an0(time(t));
end
figure('Name', 'Exercise 5')
p1 = plot(time, X31(1, :), 'r', 'LineWidth', 2); hold on; plot(time, X31(2, :), 'r', 'LineWidth', 2);
p2 = plot(time, X0(1, :), 'b', 'LineWidth', 2); hold on; plot(time, X0(2, :), 'b', 'LineWidth', 2);
title('Evolution of the analytical solution as function of time');
legend([p1, p2], {'x$_{an}$ for $\alpha = 31\deg$', 'x$_{an}$ for $\alpha = 0\deg$'});
xlabel('t'); ylabel('x$_{an}$');
grid on

% ----------------------------------------------------------------------- %
% RK2
% Plot of the function ||x_an(1) - x_RK2(1)|| = tol for alpha = pi
% to choose a proper initial value for h

% Global view of the function for tol = 1e-3 and for alpha = pi
H_g = linspace(eps, 1, 100);
FCN_g = zeros(length(H_g), 1);
A_pi = eval(subs(A, alpha, pi));
x_an = expm(A_pi*tf)*x0;
% Operator linking the solution at tf = 1 with the initial condition x0
x_RK2 = @(h) (eye(2) + (1 - floor((tf - t0)/h)*h)*A_pi + 0.5*(1 - floor((tf - t0)/h)*h)^2*A_pi^2)*...
    (eye(2) + h*A_pi + 0.5*h^2*A_pi^2)^(floor((tf - t0)/h))*x0;
fcn = @(h) max(abs(x_an - x_RK2(h))) - tol(1); %It can be represented for all
%the other values of tolerances (just change the index)
for ii = 1 : length(H_g)
    FCN_g(ii) = fcn(H_g(ii));
end
figure('Name', 'Exercise 5')
plot(H_g, FCN_g, 'LineWidth', 2);
hold on
plot(H_g, 0*ones(length(H_g)), '--r', 'LineWidth', 2); 
grid on
legend('f(h)', 'Zero line', 'Location', 'best');
xlabel('h');  ylabel('$max(|x_{an(1)} - x_{RK2(1)}|) - tol$');
title('$tol = 10^{-3}, \alpha = \pi$');

% Local view of the same function to focus better on the zero 

%Depending on the value of tolerance chosen these are the intervals of h where
%the zero is located: 
%for 1e-3 : [3*1e-2, 8*1e-2]; for 1e-4 : [1e-2, 4*1e-2]
%for 1e-5 : [3*1e-3, 8*1e-3]; for 1e-6 : [1e-3, 3*1e-3]
H = linspace(3*1e-2, 8*1e-2, 5); 
FCN = zeros(length(H), 1);

for ii = 1 : length(H)
    FCN(ii) = fcn(H(ii));
end
  
figure('Name', 'Exercise 5')
plot(H, FCN, 'LineWidth', 2);
hold on
plot(H, 0*ones(length(H)), '--r', 'LineWidth', 2); 
grid on
legend('f(h)', 'Zero line', 'Location', 'best');
xlabel('h'); ylabel('$max(|x_{an(1)} - x_{RK2(1)}|) - tol$');
title('Closer look for $tol = 10^{-3}, \alpha = \pi$');

% Search for h solution of the problem

%matrix containing values of h found for each tolerance and each value of alpha
%where every column corresponds to a different tolerance, each row corresponds
%to a different value of alpha
H = zeros(length(alpha_v), length(tol)); 
%RK2

for ii = 1 : length(tol)
    for jj = 1 : length(alpha_v)
        A_v = eval(subs(A, alpha, alpha_v(jj)));
        x_an = expm(A_v*tf)*x0;
        x_RK2 = @(h) (eye(2) + (1 - floor((tf - t0)/h)*h)*A_v + 0.5*(1 - floor((tf - t0)/h)*h)^2*A_v^2)*...
                     (eye(2) + h*A_v + 0.5*h^2*A_v^2)^(floor((tf - t0)/h))*x0;
        fcn = @(h) max(abs(x_an - x_RK2(h))) - tol(ii);                    
        h_sol = fzero(fcn, [0.001, 0.3]);  
        H(jj, ii) = h_sol;
        h0 = H(jj, ii);
    end
end

% Plot of the locus of solution in the {hlambda} plane
hlambda_r = H.*L_vreal;
hlambda_i = H.*L_vimag;
figure('Name', 'Exercise 5')
p1 = plot(hlambda_r(:, 1), hlambda_i(:, 1), 'k'); hold on; plot(hlambda_r(:, 1), -hlambda_i(:, 1), 'k');
hold on
p2 = plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b'); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b');
hold on
p3 = plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r'); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r');
hold on
p4 = plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g'); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g');
grid on
axis equal
title('$||x_{an}(1) - x_{RK2}(1)||_\infty = tol$');
legend([p1 p2 p3 p4], '$tol = 10^{-3}$', '$tol = 10^{-4}$', '$tol = 10^{-5}$', '$tol = 10^{-6}$', 'Location', 'best');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% ----- SINGLE PLOTS ----- %
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b', 'LineWidth', 2); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK2}(1)||_\infty = 10^{-4}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r', 'LineWidth', 2); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK2}(1)||_\infty = 10^{-5}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g', 'LineWidth', 2); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK2}(1)||_\infty = 10^{-6}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% Function evaluations vs tol for alpha = pi
fcn_eval_RK2 = 2*ceil((tf - t0)./H(1, :)); %2 fcn evaluation for step
figure('Name', 'Exercise 5')
loglog(tol, fcn_eval_RK2, '-or', 'LineWidth', 2);
grid on
title('Function evaluations vs tol for $\alpha = \pi$');
axis fill
xlabel('Tolerance'); ylabel('Number of fcn evaluations')


% ----------------------------------------------------------------------- %
% RK4
clear H hlambda_r hlambda_i
H = zeros(length(alpha_v), length(tol)); 

for ii = 1 : length(tol)
    for jj = 1 : length(alpha_v)
        A_v = eval(subs(A, alpha, alpha_v(jj)));
        x_an = expm(A_v*tf)*x0;
        x_RK4 = @(h) (eye(2) + (1 - floor((tf - t0)/h)*h)*A_v + 0.5*(1 - floor((tf - t0)/h)*h)^2*A_v^2 ...
        + 1/6*(1 - floor((tf - t0)/h)*h)^3*A_v^3 + 1/24*(1 - floor((tf - t0)/h)*h)^4*A_v^4) ...
        *(eye(2) + h*A_v + 0.5*h^2*A_v^2 + 1/6*h^3*A_v^3 + 1/24*h^4*A_v^4)^(floor((tf - t0)/h))*x0;
        fcn = @(h) max(abs(x_an - x_RK4(h))) - tol(ii);                    
        h_sol = fzero(fcn, [0.001, 2]);  
        H(jj, ii) = h_sol;
        h0 = H(jj, ii);
    end
end

% Plot of the locus of solution in the {hlambda} plane
hlambda_r = H.*L_vreal;
hlambda_i = H.*L_vimag;
figure('Name', 'Exercise 5')
p1 = plot(hlambda_r(:, 1), hlambda_i(:, 1), 'k'); hold on; plot(hlambda_r(:, 1), -hlambda_i(:, 1), 'k');
hold on
p2 = plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b'); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b');
hold on
p3 = plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r'); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r');
hold on
p4 = plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g'); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g');
grid on
axis equal
title('$||x_{an}(1) - x_{RK4}(1)||_\infty = tol$');
legend([p1 p2 p3 p4], '$tol = 10^{-3}$', '$tol = 10^{-4}$', '$tol = 10^{-5}$', '$tol = 10^{-6}$', 'Location', 'best');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% ----- SINGLE PLOTS ----- %
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 2), hlambda_i(:, 2), 'b', 'LineWidth', 2); hold on; plot(hlambda_r(:, 2), -hlambda_i(:, 2), 'b', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK4}(1)||_\infty = 10^{-4}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 3), hlambda_i(:, 3), 'r', 'LineWidth', 2); hold on; plot(hlambda_r(:, 3), -hlambda_i(:, 3), 'r', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK4}(1)||_\infty = 10^{-5}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
figure('Name', 'Exercise 5')
plot(hlambda_r(:, 4), hlambda_i(:, 4), 'g', 'LineWidth', 2); hold on; plot(hlambda_r(:, 4), -hlambda_i(:, 4), 'g', 'LineWidth', 2);
grid on
axis equal
title('$||x_{an}(1) - x_{RK4}(1)||_\infty = 10^{-6}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

% Function evaluations vs tol for alpha = pi
fcn_eval_RK4 = 4*ceil((tf - t0)./H(1, :)); %4 fcn evaluation for step
figure('Name', 'Exercise 5')
loglog(tol, fcn_eval_RK4, '-or', 'LineWidth', 2);
grid on
title('Function evaluations vs tol for $\alpha = \pi$');
axis fill
xlabel('Tolerance'); ylabel('Number of fcn evaluations')

%Global plot of the function evaluations 
figure('Name', 'Exercise 5')
loglog(tol, fcn_eval_RK1, '-or', 'LineWidth', 2); hold on
loglog(tol, fcn_eval_RK2, '-og', 'LineWidth', 2); hold on
loglog(tol, fcn_eval_RK4, '-ob', 'LineWidth', 2); hold on
grid on
xlabel('Tolerance'); ylabel('Number of fcn evaluations')
legend('RK1', 'RK2', 'RK4');
   
%%
%%% EXERCISE 6 %%%
clear; close all; clc; 
% All values of theta
theta = [0.1, 0.3, 0.4, 0.7, 0.9];

syms alpha 
% Matrix A in symbolic form
A = [0,  1;
     -1,  2*cos(alpha)];
% Computation of the eigenvalues (symbolic form)
L = eig(A);

% Definition of vector of alphas
alpha_v = 0: pi/90 : pi;

% Evaluation of the eigenvalues of A for every alpha
L_v = eval(subs(L, alpha, alpha_v));
L_v_real = [real(L_v(1, :)), fliplr(real(L_v(2, :)))];
L_v_imag = [imag(L_v(1, :)), fliplr(imag(L_v(2, :)))];

A = @(alpha) [ 0,            1;
              -1, 2*cos(alpha)];
          
% plot the function to understand a proper initial guess for h
h = linspace(-8, 11);
FCN = zeros(length(h), length(theta));
for jj = 1 : length(theta)
    B_BI2 = @(h, alpha) (eye(2) - (1 - theta(jj))*h*A(alpha) + 0.5*(1 - theta(jj))^2*h^2*A(alpha)^2)\...
             (eye(2) + theta(jj)*h*A(alpha) + 0.5*theta(jj)^2*h^2*A(alpha)^2);
    fcn = @(h) max(abs(eig(B_BI2(h, 0)))) - 1; 
    for hh = 1 : length(h)
    FCN(hh, jj) = fcn(h(hh));
    end
end
figure('Name', 'Exercise 6')
plot(h, FCN(:, 1), 'r', h, FCN(:, 2), 'b', h, FCN(:, 3), 'g', h, FCN(:, 4), 'm', h, FCN(:, 5), 'c', 'LineWidth', 2);
hold on
plot(h, zeros(length(h), 1), '--k', 'LineWidth', 2);
grid on
legend('$\theta = 0.1$', '$\theta = 0.3$', '$\theta = 0.4$', '$\theta = 0.7$', '$\theta = 0.9$');
xlabel('$h$'); ylabel('$f(h)$'); xlim([-7, 11]); ylim([-1, 5])
title('$f(h) = max(|eig(B_{BI2}(h)|) - 1$ for $\alpha$ = 0');

% From the analysis of the previous plot, the initial guess of h is
% deduced for each value of theta
h_guess = [2, 5, 10, -5, -2.5];
LAMBDA = zeros(length(alpha_v), length(theta));
% Plot of the stability region for BI2 for each value of theta
for jj = 1 : length(theta)
    B_BI2 = @(h, alpha) (eye(2) - (1 - theta(jj))*h*A(alpha) + 0.5*(1 - theta(jj))^2*h^2*A(alpha)^2)\...
             (eye(2) + theta(jj)*h*A(alpha) + 0.5*theta(jj)^2*h^2*A(alpha)^2);      
    [~, LAMBDA(:, jj)] = stability(alpha_v, B_BI2, h_guess(jj), L_v); 
end

figure('Name', 'Exercise 6', 'Windowstate', 'maximized')
p1 = plot([real(LAMBDA(:, 1)), real(LAMBDA(:, 1))], [imag(LAMBDA(:, 1)), -imag(LAMBDA(:, 1))], 'r', 'LineWidth', 2);
hold on
p2 = plot([real(LAMBDA(:, 2)), real(LAMBDA(:, 2))], [imag(LAMBDA(:, 2)), -imag(LAMBDA(:, 2))], 'b', 'LineWidth', 2); 
hold on
p3 = plot([real(LAMBDA(:, 3)), real(LAMBDA(:, 3))], [imag(LAMBDA(:, 3)), -imag(LAMBDA(:, 3))], 'g', 'LineWidth', 2);
hold on
p4 = plot([real(LAMBDA(:, 4)), real(LAMBDA(:, 4))], [imag(LAMBDA(:, 4)), -imag(LAMBDA(:, 4))], 'm', 'LineWidth', 2);
hold on
p5 = plot([real(LAMBDA(:, 5)), real(LAMBDA(:, 5))], [imag(LAMBDA(:, 5)), -imag(LAMBDA(:, 5))], 'c', 'LineWidth', 2);
axis equal
grid on  
title('Stability region of BI2 with different values of $\theta$');
text(1.2, 0.3, '$BI2_{0.1}$', 'Color', 'r', 'FontSize', 15);
text(3.6, 0.3, '$BI2_{0.3}$', 'Color', 'b', 'FontSize', 15);
text(8.7, 0.3, '$BI2_{0.4}$', 'Color', 'g', 'FontSize', 15);
text(-4.8, 0.3, '$BI2_{0.7}$', 'Color', 'm', 'FontSize', 15);
text(-2.2, 0.3, '$BI2_{0.9}$', 'Color', 'c', 'FontSize', 15);
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
%Plot of the stability region of BI2_0.4
figure('Name', 'Exercise 6')
fill([real(LAMBDA(:, 3)), real(LAMBDA(:, 3))], [imag(LAMBDA(:, 3)), -imag(LAMBDA(:, 3))], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
hold on
plot([real(LAMBDA(:, 3)), real(LAMBDA(:, 3))], [imag(LAMBDA(:, 3)), -imag(LAMBDA(:, 3))], 'b', 'LineWidth', 2);
hold on
plot(zeros(1, length(LAMBDA)), linspace(-6, 6, length(LAMBDA)), '--k', 'LineWidth', 2)
text(1.2, 2.5, 'Unstable', 'Color', 'b', 'FontSize', 20); 
text(8.5, 4.5, 'Stable', 'Color', rgb('MediumVioletRed'), 'FontSize', 20);
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');
title('Stability region of BI2$_{0.4}$');
grid on
axis equal
ylim([-6, 6]);
set(gca, 'Color', rgb('AntiqueWhite'));
LAMBDA_BI2 = LAMBDA;
save('Lambda_BI2', 'LAMBDA_BI2');
%%
%%% EXERCISE 7 %%%
%%% ----- ATTENTION ------ %%%
%IN ORDER FOR THIS EXERCISE TO WORK PROPERLY, IT NEEDS EXERCISE 4 AND 
%EXERCISE 6 TO BE RUNNED FIRST BECAUSE IT LOADS VARIABLES SAVED FROM THERE
%OTHERWISE IT IS NOT ABLE TO PLOT THE STABILITY REGIONS OF RK4 AND BI2
clear; close all; clc; 

%Definition of the IVP
B = [-180.5, 219.5; 179.5, -220.5];
x0 = [1, 1]';   %initial condition
t0 = 0; tf = 5; %initial and final time 
H = 0.1;        %step size
f = @(x, t) B*x;

%Analytical solution
t_vec = t0: H/100 : tf; 
x_an = @(t) expm(B*t)*x0;

X_an = zeros(length(x0), length(t_vec));
for ii = 1 : length(t_vec)
    X_an(:, ii) = x_an(t_vec(ii));  %[2 x 1] vector
end

%Solution of the IVP with RK4
[X_RK4, t] = RK4(f, H, x0, t0, tf);

% --- Plot of the analytical solution --- %
%first component
figure('Name', 'Exercise 7');
plot(t_vec, X_an(1, :), 'b', 'LineWidth', 2);
xlabel('t'); ylabel('$x_1$'); grid on
title('Analytical solution for x1');
ax=axes;
set(ax,'units','normalized','position',[0.3,0.6,0.2,0.2])
box(ax,'on')
plot(t_vec, X_an(1, :), 'b', 'LineWidth', 2, 'parent', ax)
set(ax,'xlim',[0,0.15],'ylim',[1, 1.1], 'FontSize', 10)
%second component
figure('Name', 'Exercise 7');
plot(t_vec, X_an(2, :), 'b', 'LineWidth', 2); 
xlabel('t'); ylabel('$x_2$'); grid on
title('Analytical solution for x2');
ax=axes;
set(ax,'units','normalized','position',[0.3,0.45,0.1,0.4])
box(ax,'on')
plot(t_vec, X_an(2, :), 'b', 'LineWidth', 2, 'parent', ax)
set(ax,'xlim',[0,0.01],'ylim',[0.88, 1], 'FontSize', 10)

% --- Plot of RK4 solution against analytical solution --- %
%first component
figure('Name', 'Exercise 7', 'WindowState', 'maximized')
plot(t_vec, X_an(1, :), 'b', 'LineWidth', 2);
hold on
plot(t, X_RK4(1, :), '-or', 'LineWidth', 2);
grid on
title('Comparison between RK4 and analytical solutions');
legend('Analytical Solution', 'RK4 Solution', 'Location', 'best');
xlabel('t'); ylabel('$x_1$');
%second component
figure('Name', 'Exercise 7', 'WindowState', 'maximized');
plot(t_vec, X_an(2, :), 'b', 'LineWidth', 2);
hold on
plot(t, X_RK4(2, :), '-or', 'LineWidth', 2);
grid on
title('Comparison between RK4 and analytical solutions');
legend('Analytical Solution', 'RK4 Solution', 'Location', 'best');
xlabel('t'); ylabel('$x_2$');

%Solution of the IVP with BI2_0.1
theta = 0.1;
[X_BI2, t] = BI2(theta, H, x0, B, t0, tf);
% --- Plot of BI2_0.1 solution against analytical solution --- %
%first component
figure('Name', 'Exercise 7', 'Windowstate', 'maximized')
plot(t_vec, X_an(1, :), 'b', 'LineWidth', 2);
hold on
plot(t, X_BI2(1, :), '-og', 'LineWidth', 2);
xlabel('t'); ylabel('$x_1$'); grid on
legend('Analytical Solution', 'BI2 Solution', 'Location', 'best');
title('Comparison between BI2 and analytical solutions');
ax=axes;
set(ax,'units','normalized','position',[0.3,0.5,0.3,0.3])
box(ax,'on')
plot(t_vec, X_an(1, :), 'b', t, X_BI2(1, :), '-og', 'LineWidth', 2, 'parent', ax)
set(ax,'xlim',[0,0.15],'ylim',[0.9, 1.1], 'FontSize', 10)
%second component
figure('Name', 'Exercise 7', 'Windowstate', 'maximized')
plot(t_vec, X_an(2, :), 'b', 'LineWidth', 2);
hold on
plot(t, X_BI2(2, :), '-og', 'LineWidth', 2);
xlabel('t'); ylabel('$x_2$'); grid on
legend('Analytical Solution', 'BI2 Solution', 'Location', 'best');
title('Comparison between BI2 and analytical solutions');
ax=axes;
set(ax,'units','normalized','position',[0.3,0.45,0.1,0.4])
box(ax,'on')
plot(t_vec, X_an(2, :), 'b', t, X_BI2(2, :), '-og', 'LineWidth', 2, 'parent', ax)
set(ax,'xlim',[0,0.1],'ylim',[0.8, 1], 'FontSize', 10)

%Plot of the eigevanlues of B in the (h lambda)-plane for RK4 and for BI2
b = eig(B);

%RK4
load('Lambda_RK4');
figure('Name', 'Exercise 7', 'WindowState', 'maximized')
plot([real(LAMBDA_RK4), real(LAMBDA_RK4)], [imag(LAMBDA_RK4), -imag(LAMBDA_RK4)], 'k', 'LineWidth', 2)
hold on
fill([real(LAMBDA_RK4), real(LAMBDA_RK4)], [imag(LAMBDA_RK4), -imag(LAMBDA_RK4)], 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
scatter(real(H*b), imag(H*b), 60, 'filled', 'b'); hold on
text(-23, 12, 'Unstable', 'Color', 'r', 'FontSize', 25); hold on
annotation('textarrow', [.7 .75], [.35 .5], 'String', '\lambda_1', 'FontSize', 25, 'Color', 'b'); hold on
annotation('textarrow', [.33 .28], [.35 .5], 'String', '\lambda_2', 'FontSize', 25, 'Color', 'b'); hold on
axis equal
grid on
set(gca, 'Color', rgb('AntiqueWhite'));
title('RK4 stability region and mapping of $\lambda_{Bi}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');

%BI2
load('Lambda_BI2');
figure('Name', 'Exercise 7', 'Windowstate', 'maximized')
plot([real(LAMBDA_BI2(:, 1)), real(LAMBDA_BI2(:, 1))], [imag(LAMBDA_BI2(:, 1)), -imag(LAMBDA_BI2(:, 1))], 'k', 'LineWidth', 2);
hold on
fill([real(LAMBDA_BI2(:, 1)), real(LAMBDA_BI2(:, 1))], [imag(LAMBDA_BI2(:, 1)), -imag(LAMBDA_BI2(:, 1))], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.4);
hold on
scatter(real(H*b), imag(H*b), 60, 'filled', 'b');
text(-23, 12, 'Stable', 'Color', 'b', 'FontSize', 25); hold on
annotation('textarrow', [.68 .73], [.35 .5], 'String', '\lambda_1', 'FontSize', 25, 'Color', 'b'); hold on
annotation('textarrow', [.33 .28], [.35 .5], 'String', '\lambda_2', 'FontSize', 25, 'Color', 'b'); hold on
axis equal 
grid on
set(gca, 'Color', rgb('PaleTurquoise'));
title('$BI2_{0.1}$ stability region and mapping of $\lambda_{Bi}$');
xlabel('$\Re(h\lambda)$'); ylabel('$\Im(h\lambda)$');


%%
%%% --- FUNCTIONS --- %%%

function [x, n_iter, n_fevals] = bisection(a, b, fcn, nmax, tol)
%This function applies the bisection method in order to find the root
%of the function fcn

%INPUT: 
% a   : left extreme of the considered interval
% b   : right extreme of the considered interval 
% fcn : function for which the zero shall be found
% nmax: maximum number of iterations allowed
% tol : tolerance value to assure 8-digit accuracy

%OUTPUT:
% x      : final solution found by bisection method
% n_iter : number of iterations required by the method to arrive at x
% n_fevals : number of function evaluations required by the method


%check validity of input values a and b for the application of the
%bisection method
if fcn(a)*fcn(b) > 0
    error('Error in values of a and b. Bisection method not applicable, change values\n');
end

n_iter = -1;
n_fevals = 2; %f(a) and f(b)
%nmax = ceil( log2((b-a)/tol)-1);   %maximum number of iterations allowed
err = tol + 1;

while n_iter <= nmax && err > tol      
    x = (a + b)/2;    
    n_iter = n_iter + 1;
    err = abs(fcn(x));       

            if  fcn(x) == 0
                n_fevals = n_fevals + 1;
               return
            elseif fcn(a)*fcn(x) < 0
                b = x;
            else 
                a = x; 
            end
    n_fevals = n_fevals + 1; 
 end

end

function [x, n_iter, n_fevals] = secant(x0, x1, fcn, tol, nmax)
%This function aims at finding the root of the function fcn by applying the
%secant method
%INPUT: 
% x0   : initial guess near the root
% x1   : second initial guess near the root
% fcn  : function for which the zero shall be found
% tol  : tolerance value to assure 8-digit accuracy
% nmax : maximum number of iterations allowed

%OUTPUT:
% x        : final solution found by secant method
% n_iter   : number of iterations required by the method to arrive at x
% n_fevals : number of function evaluations required by the method

n_iter = -1;
err = tol + 1;
x_k = x1;
x_prev = x0;
n_fevals = 1;  %f(x_prev)

while n_iter <= nmax && err > tol
    n_iter = n_iter + 1;
    x = x_k - (x_k - x_prev)/(fcn(x_k) - fcn(x_prev))*fcn(x_k);
    err = abs(fcn(x));
    x_prev = x_k;
    x_k = x;
    n_fevals = n_fevals + 1; %f(x_k)
end
end

function [x, n_iter, n_fevals] = regula_falsi(x0, x1, fcn, tol, nmax)
%This function aims at finding the root of the function fcn by applying the
%regula falsi method

%INPUT: 
% x0   : initial guess near the root
% x1   : second initial guess near the root
% fcn  : function for which the zero shall be found
% tol  : tolerance value to assure 8-digit accuracy
% nmax : maximum number of iterations allowed

%OUTPUT:
% x        : final solution found by regula falsi method
% n_iter   : number of iterations required by the method to arrive at x
% n_fevals : number of function evaluations required by the method

n_iter = -1;
err = tol + 1;
x_k = x1;
x_prev = x0;

if fcn(x0)*fcn(x1) > 0
    error('Error in values of x0 and x1. Regula falsi method not applicable, change values\n');
end
n_fevals = 2; %f(x0) and f(x1)

while n_iter < nmax && err > tol
    n_iter = n_iter + 1;
    x = x_k - (x_k - x_prev)/(fcn(x_k) - fcn(x_prev))*fcn(x_k);
    err = abs(fcn(x));
    if fcn(x) == 0
        n_fevals = n_fevals + 1;
        return
    elseif fcn(x_prev)*fcn(x) < 0
        x_k = x;
    else
        x_prev = x; 
    end
    n_fevals = n_fevals + 1;
end
end

function [x, n_iter, err] = newton_analytical(f, x_var, x0, nmax, tol)
%The aim of this function is to compute the root of the function f through 
%Newton's method where the partial derivatives are computed analytically.
%To use this function for any input function f substitute to J the
%following: J = jacobian(f).

%INPUT: 
% f      : vectorial function in symbolic form
% x_var  : vector containing the symbolic variables [Nx1]
% x0     : vector of initial guesses [Nx1]
% nmax   : maximum number of iterations
% tol    : digit accuracy

%OUTPUT: 
% x      : root of the function [Nx1]
% n_iter : number of iterations done by Newton's method
% err    : accuracy defined as |x_end - x_previous|

%Jacobian
J = [2*x_var(1)-1, -1; x_var(1)/8, 2*x_var(2)]; %J = jacobian(f)
n_iter = -1;
err = [tol +1, tol + 1]; %in this way it enters for sure in the cycle

while n_iter < nmax && any(err > tol)
    n_iter = n_iter + 1;
    x = x0 - eval(subs(J\f, x_var, x0)); %Newton's iterate
    err = abs(x - x0);
    x0 = x;    
end
end

    
function [x_sol, n_iter, err] = newton_fd(F, x_var, x0, nmax, tol)
%The aim of this function is to compute the root of the function f through 
%Newton's method where the partial derivatives are approximated
%with the forward finite difference.

%INPUT: 
% F     : vectorial function in symbolic form
% x_var : vector containing the symbolic variables [Nx1]
% x0    : guess vector [Nx1]
% nmax  : maximum number of iterations
% tol   : digit accuracy

%OUTPUT:
% x_sol : root vector [Nx1]
% n_iter: number of iterations performed by the method
% err   : accuracy defined as |x_end - x_previous|

n_iter = -1;
n = length(F);
J = zeros(n, n);
xk = x0;  
xk_pert = x0; 
err = [tol + 1, tol + 1]; %in this way it enters in the cycle for sure
while n_iter < nmax && any(err > tol)
    n_iter = n_iter + 1;
    % Construction of the Jacobian matrix computing the derivatives with
    % forward difference
    for ii = 1 : n
        f = F(ii);
        for jj = 1 : n
            x = xk(jj);
            epsilon = 1e-6;%max(sqrt(eps), sqrt(eps)*abs(x)); %perturbation
            g = eval(subs(f, x_var, xk));  
            xk_pert(jj) = x + epsilon;
            h = eval(subs(f, x_var, xk_pert));
            J(ii, jj) = (h - g)/epsilon; %forward difference
            xk_pert = xk;
        end
    end
    xkk = xk - J\eval(subs(F, x_var, xk)); %Newton's iterate
    err = abs(xkk - xk);
    xk = xkk;
    xk_pert = xk;
end
x_sol = xkk;
end

function [x_sol, n_iter, err] = newton_bd(F, x_var, x0, nmax, tol)
%The aim of this function is to compute the root of the function f through 
%Newton's method where the partial derivatives are approximated
%with the backward finite difference.

%INPUT: 
% F     : vectorial function in symbolic form
% x_var : vector containing the symbolic variables [Nx1]
% x0    : guess vector [Nx1]
% nmax  : maximum number of iterations
% tol   : digit accuracy

%OUTPUT:
% x_sol : root vector [Nx1]
% n_iter: number of iterations performed by the method
% err   : accuracy defined as |x_end - x_previous|

n_iter = -1;
n = length(F);
J = zeros(n, n);
xk = x0;  
xk_pert = x0;
err = [tol + 1, tol+ 1]; %in this way it enters in the cycle for sure
while n_iter < nmax && any(err > tol)
    n_iter = n_iter + 1;
    % Construction of the Jacobian matrix computing the derivatives with
    % backward difference
    for ii = 1 : n
        f = F(ii);
        for jj = 1 : n
            x = xk(jj);
            epsilon = max(sqrt(eps), sqrt(eps)*abs(x));
            g = eval(subs(f, x_var, xk));
            xk_pert(jj) = x - epsilon;
            h = eval(subs(f, x_var, xk_pert));
            J(ii, jj) = (g - h)/epsilon; %backward difference
            xk_pert = xk;
        end
    end
    xkk = xk - J\eval(subs(F, x_var, xk));
    err = abs(xkk - xk);
    xk = xkk;
    xk_pert = xk;
end
x_sol = xkk;
end

function [x_sol, n_iter, err] = newton_cd(F, x_var, x0, nmax, tol)
%The aim of this function is to compute the root of the function f through 
%Newton's method where the partial derivatives are approximated
%with the centered finite difference.

%INPUT: 
% F     : vectorial function in symbolic form
% x_var : vector containing the symbolic variables [Nx1]
% x0    : guess vector [Nx1]
% nmax  : maximum number of iterations
% tol   : digit accuracy

%OUTPUT:
% x_sol : root vector [Nx1]
% n_iter: number of iterations performed by the method
% err   : accuracy defined as |x_end - x_previous|
n_iter = -1;
n = length(F);
J = zeros(n, n);
xk = x0;
xk_pertf = x0;
xk_pertb = x0;
err = [1, 1];
while n_iter < nmax && any(err > tol)
    n_iter = n_iter + 1;
    % Construction of the Jacobian matrix computing the derivatives with
    % centered differences
    for ii = 1 : n
        f = F(ii);
        for jj = 1 : n
            x = xk(jj);
            epsilon = max(sqrt(eps), sqrt(eps)*abs(x));
            xk_pertf(jj) = x + epsilon;
            xk_pertb(jj) = x - epsilon;
            g = eval(subs(f, x_var, xk_pertf));
            h = eval(subs(f, x_var, xk_pertb));
            J(ii, jj) = (g - h)/(2*epsilon);  %centered difference approx
            xk_pertf = xk; 
            xk_pertb = xk;
        end
    end
    xkk = xk - J\eval(subs(F, x_var, xk));
    err = abs(xkk - xk);
    xk = xkk;
    xk_pertf = xk;
    xk_pertb = xk;
end
x_sol = xkk;
end


function [x, t] = heun(f, h, x0, t_i, t_f)
%This function aims at solving an ODE with a given initial value (IVP)
%adopting the Heun's method. 
%INPUT: 
%f   : function f(x, t). It can be also a vector-valued function
%h   : time step
%x0  : initial condition 
%t_i : initial time 
%t_f : final time

%OUTPUT:
%x : solution vector at each time instant [1 x n]
%t : [t_i : h : t_f]

t = t_i : h : t_f;
xk = x0;
x = zeros(length(x0), length(t));
x(:, 1) = x0;

for ii = 1 : length(t) - 1
    x_P = xk + h*f(xk, t(ii));
    x(:, ii+1) = xk + 0.5*h*(f(xk, t(ii)) + f(x_P, t(ii + 1)));
    xk = x(:, ii+1);
end
end      

function [x, t] = RK4(f, h, x0, t_i, t_f)
%This function aims at solving an ODE with a given initial value (IVP)
%adopting the Heun's method. 
%INPUT: 
%f   : function f(x, t). It can be also a vector-valued function
%h   : time step
%x0  : initial condition 
%t_i : initial time 
%t_f : final time

%OUTPUT:
%x : solution vector at each time instant [1 x n]
%t : [t_i : h : t_f]

t = t_i : h : t_f;
xk = x0;
x = zeros(length(x0), length(t));
x(:, 1) = x0;

for ii = 1 : length(t) - 1
    k1 = f(xk, t(ii));
    k2 = f(xk + 0.5*h*k1, t(ii) + 0.5*h);
    k3 = f(xk + 0.5*h*k2, t(ii) + 0.5*h);
    k4 = f(xk + h*k3, t(ii) + h);
    x(:, ii + 1) = xk + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    xk = x(:, ii + 1);
end
end

function [H_sol, LAMBDA] = stability(alpha_v, F, h0, lambda_A)
%This function aims at finding the points that trace the contour of the
%stability region of the method described by the operator F.
%To do this, the function uses fsolve
%INPUT: 
% alpha_v   : vector of values of alpha [1 x N]
% F         : operator of the method
% h0        : initial guess
% lambda_A  : vector of the complex conjugate eigenvalues of A [2 x N]

%OUTPUT:
% H_sol     : solution of the function max(|eig(F)|) for each value of
%             alpha
% LAMBDA    : Points to plot in the h-lambda plane
 
LAMBDA = zeros(length(alpha_v),1); %vector containing h*lambda
H_sol = zeros(length(alpha_v),1); %vector containing values of h for each alpha
H_sol(1) = h0;
options = optimoptions('fsolve', 'Display', 'off', 'StepTolerance', 1e-10);
 for ii = 1 : length(alpha_v)
     fcn = @(h) max(abs(eig(F(h, alpha_v(ii))))) - 1; 
     H_sol(ii) = fsolve(fcn, h0, options);  
     LAMBDA(ii) = H_sol(ii)*lambda_A(1, ii);
     h0 = H_sol(ii);
 end
end

function [H_sol, LAMBDA] = stability_fzero(alpha_v, F, h0, lambda_A)
%This function aims at finding the points that trace the contour of the
%stability region of the method described by the operator F.
%To do this, the function uses fzero

%INPUT: 
% alpha_v   : vector of values of alpha [1 x N]
% F         : operator of the method, function of h and alpha
% h0        : initial guess
% lambda_A  : vector of the complex conjugate eigenvalues of A [2 x N]

%OUTPUT:
% H_sol     : solution of the function max(|eig(F)|) for each value of
%             alpha
% LAMBDA    : Points to plot in the h-lambda plane
 
LAMBDA = zeros(length(alpha_v),1); %vector containing h*lambda
H_sol = zeros(length(alpha_v),1); %vector containing values of h for each alpha
H_sol(1) = h0;
 for ii = 1 : length(alpha_v)
     fcn = @(h) max(abs(eig(F(h, alpha_v(ii))))) - 1; 
     H_sol(ii) = fzero(fcn, h0);
     LAMBDA(ii) = H_sol(ii)*lambda_A(1, ii);
     h0 = H_sol(ii);
 end
end

function [X_BI2, t_vec] = BI2(theta, H, x0, B, t0, tf)

%This function allows to find the approximate solution of the IVP
%x_dot = B*x by applying backinterpolation method of the second order

%INPUT: 
%theta : portion of the step size for the theta-method
%H     : step size
%x0    : initial condition
%B     : matrix 
%t0    : initial time
%tf    : final time

%OUTPUT: 
%X_BI2 : solution of BI2_theta
%t_vec : vector of time instants

%Time vector
t_vec = t0 : H : tf;

%Operator
B_BI2 = @(h) (eye(2) - (1 - theta)*h*B + 0.5*(1 - theta)^2*h^2*B^2)\...
             (eye(2) + theta*h*B + 0.5*theta^2*h^2*B^2);
         
X_BI2 = zeros(length(x0), length(t_vec));
xk = x0;
X_BI2(:, 1) = x0;
for ii = 1 : length(t_vec)-1
    X_BI2(:, ii+1) = B_BI2(H)*xk;
    xk = X_BI2(:, ii+1);
end
end


%Function downloaded to have more colors for the plot
function rgb = rgb(s)
  persistent num name
  if isempty(num) % First time rgb is called
    [num,name] = getcolors();
    name = lower(name);
    num = reshape(hex2dec(num), [], 3);
    % Divide most numbers by 256 for "aesthetic" reasons (green=[0 0.5 0])
    I = num < 240;  % (interpolate F0--FF linearly from 240/256 to 1.0)
    num(I) = num(I)/256;
    num(~I) = ((num(~I) - 240)/15 + 15)/16; + 240;
  end
  if strcmpi(s,'chart')
    showcolors()
  else
    k = find(strcmpi(s, name));
    if isempty(k)
      error(['Unknown color: ' s]);
    else
      rgb = num(k(1), :);
    end
  end
end

function showcolors()
  [num,name] = getcolors();
  grp = {'White', 'Gray', 'Red', 'Pink', 'Orange', 'Yellow', 'Brown'...
    , 'Green', 'Blue', 'Purple', 'Grey'};
  J = [1,3,6,8,9,10,11];
  fl = lower(grp);
  nl = lower(name);
  for i=1:length(grp)
    n(i) = strmatch(fl{i}, nl, 'exact'); 
  end
  clf
  p = get(0,'screensize');
  wh = 0.6*p(3:4);
  xy0 = p(1:2)+0.5*p(3:4) - wh/2;
  set(gcf,'position', [xy0 wh]);
  axes('position', [0 0 1 1], 'visible', 'off');
  hold on
  x = 0;
  N = 0;
  for i=1:length(J)-1
    N = max(N, n(J(i+1)) - n(J(i)) + (J(i+1) - J(i))*1.3); 
  end
  h = 1/N;
  w = 1/(length(J)-1);
  d = w/30;
  for col = 1:length(J)-1;
    y = 1 - h;
    for i=J(col):J(col+1)-1
      t = text(x+w/2, y+h/10 , [grp{i} ' colors']);
      set(t, 'fontw', 'bold', 'vert','bot', 'horiz','cent', 'fontsize',10);
      y = y - h;
      for k = n(i):n(i+1)-1
        c = rgb(name{k});
        bright = (c(1)+2*c(2)+c(3))/4;
        if bright < 0.5, txtcolor = 'w'; else txtcolor = 'k'; end
        rectangle('position',[x+d,y,w-2*d,h],'facecolor',c);
        t = text(x+w/2, y+h/2, name{k}, 'color', txtcolor);
        set(t, 'vert', 'mid', 'horiz', 'cent', 'fontsize', 9);
        y = y - h;
      end
      y = y - 0.3*h;
    end
    x = x + w;
  end
end

function [hex,name] = getcolors()
  css = {
    %White colors
    'FF','FF','FF', 'White'
    'FF','FA','FA', 'Snow'
    'F0','FF','F0', 'Honeydew'
    'F5','FF','FA', 'MintCream'
    'F0','FF','FF', 'Azure'
    'F0','F8','FF', 'AliceBlue'
    'F8','F8','FF', 'GhostWhite'
    'F5','F5','F5', 'WhiteSmoke'
    'FF','F5','EE', 'Seashell'
    'F5','F5','DC', 'Beige'
    'FD','F5','E6', 'OldLace'
    'FF','FA','F0', 'FloralWhite'
    'FF','FF','F0', 'Ivory'
    'FA','EB','D7', 'AntiqueWhite'
    'FA','F0','E6', 'Linen'
    'FF','F0','F5', 'LavenderBlush'
    'FF','E4','E1', 'MistyRose'
    %Grey colors'
    '80','80','80', 'Gray'
    'DC','DC','DC', 'Gainsboro'
    'D3','D3','D3', 'LightGray'
    'C0','C0','C0', 'Silver'
    'A9','A9','A9', 'DarkGray'
    '69','69','69', 'DimGray'
    '77','88','99', 'LightSlateGray'
    '70','80','90', 'SlateGray'
    '2F','4F','4F', 'DarkSlateGray'
    '00','00','00', 'Black'
    %Red colors
    'FF','00','00', 'Red'
    'FF','A0','7A', 'LightSalmon'
    'FA','80','72', 'Salmon'
    'E9','96','7A', 'DarkSalmon'
    'F0','80','80', 'LightCoral'
    'CD','5C','5C', 'IndianRed'
    'DC','14','3C', 'Crimson'
    'B2','22','22', 'FireBrick'
    '8B','00','00', 'DarkRed'
    %Pink colors
    'FF','C0','CB', 'Pink'
    'FF','B6','C1', 'LightPink'
    'FF','69','B4', 'HotPink'
    'FF','14','93', 'DeepPink'
    'DB','70','93', 'PaleVioletRed'
    'C7','15','85', 'MediumVioletRed'
    %Orange colors
    'FF','A5','00', 'Orange'
    'FF','8C','00', 'DarkOrange'
    'FF','7F','50', 'Coral'
    'FF','63','47', 'Tomato'
    'FF','45','00', 'OrangeRed'
    %Yellow colors
    'FF','FF','00', 'Yellow'
    'FF','FF','E0', 'LightYellow'
    'FF','FA','CD', 'LemonChiffon'
    'FA','FA','D2', 'LightGoldenrodYellow'
    'FF','EF','D5', 'PapayaWhip'
    'FF','E4','B5', 'Moccasin'
    'FF','DA','B9', 'PeachPuff'
    'EE','E8','AA', 'PaleGoldenrod'
    'F0','E6','8C', 'Khaki'
    'BD','B7','6B', 'DarkKhaki'
    'FF','D7','00', 'Gold'
    %Brown colors
    'A5','2A','2A', 'Brown'
    'FF','F8','DC', 'Cornsilk'
    'FF','EB','CD', 'BlanchedAlmond'
    'FF','E4','C4', 'Bisque'
    'FF','DE','AD', 'NavajoWhite'
    'F5','DE','B3', 'Wheat'
    'DE','B8','87', 'BurlyWood'
    'D2','B4','8C', 'Tan'
    'BC','8F','8F', 'RosyBrown'
    'F4','A4','60', 'SandyBrown'
    'DA','A5','20', 'Goldenrod'
    'B8','86','0B', 'DarkGoldenrod'
    'CD','85','3F', 'Peru'
    'D2','69','1E', 'Chocolate'
    '8B','45','13', 'SaddleBrown'
    'A0','52','2D', 'Sienna'
    '80','00','00', 'Maroon'
    %Green colors
    '00','80','00', 'Green'
    '98','FB','98', 'PaleGreen'
    '90','EE','90', 'LightGreen'
    '9A','CD','32', 'YellowGreen'
    'AD','FF','2F', 'GreenYellow'
    '7F','FF','00', 'Chartreuse'
    '7C','FC','00', 'LawnGreen'
    '00','FF','00', 'Lime'
    '32','CD','32', 'LimeGreen'
    '00','FA','9A', 'MediumSpringGreen'
    '00','FF','7F', 'SpringGreen'
    '66','CD','AA', 'MediumAquamarine'
    '7F','FF','D4', 'Aquamarine'
    '20','B2','AA', 'LightSeaGreen'
    '3C','B3','71', 'MediumSeaGreen'
    '2E','8B','57', 'SeaGreen'
    '8F','BC','8F', 'DarkSeaGreen'
    '22','8B','22', 'ForestGreen'
    '00','64','00', 'DarkGreen'
    '6B','8E','23', 'OliveDrab'
    '80','80','00', 'Olive'
    '55','6B','2F', 'DarkOliveGreen'
    '00','80','80', 'Teal'
    %Blue colors
    '00','00','FF', 'Blue'
    'AD','D8','E6', 'LightBlue'
    'B0','E0','E6', 'PowderBlue'
    'AF','EE','EE', 'PaleTurquoise'
    '40','E0','D0', 'Turquoise'
    '48','D1','CC', 'MediumTurquoise'
    '00','CE','D1', 'DarkTurquoise'
    'E0','FF','FF', 'LightCyan'
    '00','FF','FF', 'Cyan'
    '00','FF','FF', 'Aqua'
    '00','8B','8B', 'DarkCyan'
    '5F','9E','A0', 'CadetBlue'
    'B0','C4','DE', 'LightSteelBlue'
    '46','82','B4', 'SteelBlue'
    '87','CE','FA', 'LightSkyBlue'
    '87','CE','EB', 'SkyBlue'
    '00','BF','FF', 'DeepSkyBlue'
    '1E','90','FF', 'DodgerBlue'
    '64','95','ED', 'CornflowerBlue'
    '41','69','E1', 'RoyalBlue'
    '00','00','CD', 'MediumBlue'
    '00','00','8B', 'DarkBlue'
    '00','00','80', 'Navy'
    '19','19','70', 'MidnightBlue'
    %Purple colors
    '80','00','80', 'Purple'
    'E6','E6','FA', 'Lavender'
    'D8','BF','D8', 'Thistle'
    'DD','A0','DD', 'Plum'
    'EE','82','EE', 'Violet'
    'DA','70','D6', 'Orchid'
    'FF','00','FF', 'Fuchsia'
    'FF','00','FF', 'Magenta'
    'BA','55','D3', 'MediumOrchid'
    '93','70','DB', 'MediumPurple'
    '99','66','CC', 'Amethyst'
    '8A','2B','E2', 'BlueViolet'
    '94','00','D3', 'DarkViolet'
    '99','32','CC', 'DarkOrchid'
    '8B','00','8B', 'DarkMagenta'
    '6A','5A','CD', 'SlateBlue'
    '48','3D','8B', 'DarkSlateBlue'
    '7B','68','EE', 'MediumSlateBlue'
    '4B','00','82', 'Indigo'
    %Gray repeated with spelling grey
    '80','80','80', 'Grey'
    'D3','D3','D3', 'LightGrey'
    'A9','A9','A9', 'DarkGrey'
    '69','69','69', 'DimGrey'
    '77','88','99', 'LightSlateGrey'
    '70','80','90', 'SlateGrey'
    '2F','4F','4F', 'DarkSlateGrey'
    };
  hex = css(:,1:3);
  name = css(:,4);
end