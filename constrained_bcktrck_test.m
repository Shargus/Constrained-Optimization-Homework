clear 
close all
clc

%% DEFINE THE PROBLEM

n = 10^3;   % n. of variables
f = @(x) [1:n] * (x.^2);     % f. handle for the objective function
gradf = @(x) [1:n]' .* (2*x);    % f. handle for the gradient of the objective function
box_mins = -5.12 * ones(n,1);   % upper boundaries of the feasible box
box_maxs = 5.12 * ones(n,1);    % lower boundaries of the feasible box



%% DEFINE PARAMETERS FOR THE PROJECTED GRADIENT METHOD

rng(0);    % set seed (for replicability)
x0 = unifrnd(-10,10,n,1);

% for the projection step
gamma = 0.9; % 0.8, 0.9, 1
Pi_X = @(x) max(min(x, box_maxs), box_mins);  % f. handle for the box projection function

% For the purpose of comparing the projected gradient method and the
% unconstrained gradient method, here is a function handle for an identity
% projection function
%Pi_X = @(x) x;

% for the backtracking line search
c1 = 1e-4;
rho = 0.8;
btmax = 100;

% for the stopping criterion
tolgrad = 1e-5;     % iterate until ||gradf(x(k))|| < tolgrad
tolx = 1e-5;        % iterate until ||x(k+1)-x(k)|| < tolx
kmax = 3000;

% Choose how to compute the gradient
% - 'forw'  --> approximate grad(f) with the forward difference approximation
% - 'backw' --> approximate grad(f) with the backward difference approximation
% - 'centr' --> approximate grad(f) with the centered difference approximation
% - 'exact' (or just every other string) --> use the exact grad(f)
FDgrad = 'exact';
k_FDgrad = 2; % the exponent which appears in the expression of the increment h=10^(-k)*||x||

verbose = true;     % if true, print the current iteration every 10



%% RUN THE PGM

disp('**** CONSTR. STEEPEST DESCENT: START *****')
tic     % start timer

[xk_final, fk_final, gradfk_norm_final, deltaxk_norm_final, k_final, ...
    xseq_final, btseq_final, projection_count] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, kmax, tolgrad, c1, rho, ...
    btmax, gamma, tolx, Pi_X, verbose, FDgrad, k_FDgrad);

disp('**** CONSTR. STEEPEST DESCENT: FINISHED *****')
disp('')

disp('**** CONSTR. STEEPEST DESCENT: RESULTS *****')
disp('************************************')
toc     % stop timer
if n==2
    disp(['xk: ', mat2str(xk_final), ' (actual minimizer: vector of 0s);'])
end
disp(['norm(xk): ', mat2str(norm(xk_final)), ' (actual norm of the minimizer: 0; NB: this is also the norm of the error, ||xk-x*||);'])
disp(['f(xk): ', num2str(fk_final), ' (actual min. value: 0);'])
disp(['N. of iterations: ', num2str(k_final),'/',num2str(kmax), ';'])
disp(['N. of projections: ', num2str(projection_count), ';'])
disp(['gradient norm: ', num2str(gradfk_norm_final), ';'])
disp(['length of last step: ', num2str(deltaxk_norm_final), ';'])
disp('************************************')



%% PLOTS

% Barplot of btseq
fig_bt_iters = figure();
bar(btseq_final)
title('N. of backtracking shrinkings of alpha per iteration')

if n==2     % in this very simple case, we can visualize the objective function in 3D

    % Creation of the data to plot the domain boundaries
    t = linspace(0, 1, 25);
    dom_xy_1 = box_mins + t .* ([box_mins(1); box_maxs(2)] - box_mins);
    dom_xy_2 = [box_mins(1); box_maxs(2)] + t .* (box_maxs - [box_mins(1); box_maxs(2)]);
    dom_xy_3 = box_maxs + t .* ([box_maxs(2); box_mins(1)] - box_maxs);
    dom_xy_4 = [box_maxs(2); box_mins(1)] + t .* (box_mins - [box_maxs(2); box_mins(1)]);
    
    dom_xy = [dom_xy_1, dom_xy_2, dom_xy_3, dom_xy_4];
    f_z = f(dom_xy);
    
    % Projection of the starting point
    Pi_X_x0 = Pi_X(x0);
    
    % Creation of the meshgrid for the contour-plot
    [X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 6, 500));
    
    % Computation of the values of f for each point of the mesh
    Z = X.^2 + 2*Y.^2;

    % Simple Plot
    fig_contour = figure();
    % Contour plot with curve levels for each point in xseq
    [C1, ~] = contour(X, Y, Z);
    hold on
    % plot of the points in xseq
    plot([x0(1), Pi_X_x0(1)], [x0(2), Pi_X_x0(2)], 'r--*')
    plot([Pi_X_x0(1) xseq_final(1, :)], [Pi_X_x0(2) xseq_final(2, :)], 'b--*')
    plot(dom_xy(1, :), dom_xy(2, :), 'k')
    hold off
    title('PGM - x^2 + 2y^2')
    
    % Much more interesting plot
    fig_surface = figure();
    surf(X, Y, Z,'EdgeColor','none')
    hold on
    plot3([x0(1) Pi_X_x0(1)], [x0(2) Pi_X_x0(2)], [f(x0), f(Pi_X_x0)], 'y--*')
    plot3([Pi_X_x0(1) xseq_final(1, :)], [Pi_X_x0(2) xseq_final(2, :)], [f(Pi_X_x0), f(xseq_final)], 'r--*')
    plot3(dom_xy(1, :), dom_xy(2, :), f_z, 'k')
    hold off
    title('PGM - x^2 + 2y^2')

end
