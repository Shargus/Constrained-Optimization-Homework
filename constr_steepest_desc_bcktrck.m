function [xk, fk, gradfk_norm, deltaxk_norm, k, xseq, btseq, projection_count] = ...
    constr_steepest_desc_bcktrck(x0, f, gradf, ...
    kmax, tolgrad, c1, rho, btmax, gamma, tolx, Pi_X, verbose, FDgrad, k_FDgrad)
%
% Projected gradient method (steepest descent) for constrained optimization.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
%           gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, less than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
%         backtracking strategy.
% gamma = the initial factor that multiplies the descent direction at each
%         iteration;
% tolx = value used as stopping criterion w.r.t. the norm of the
%        steps (xnew - xk); infact, in constrained optimization the global
%        minimizer could also be a non-stationary point, so the stopping
%        criterion based on a tolerance on the norm of the gradient is not
%        sufficient attained in a non-stationary point
% Pi_X = projection function
% verbose = if true, print the current iteration every 10
% FDgrad = choose how to compute the gradient (exactly or approximately)
% k_FDgrad = the exponent which appears in the expression of the increment
%            h=10^(-k_FDgrad)*||x||
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% deltaxk_norm = length of the last step of the sequence
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.

switch FDgrad
    case 'forw'
        % overwrite gradf with a f. handle that uses the forward difference
        % approximation
        h = @(x) (10^(-k_FDgrad))*norm(x);
        n = size(x0,1);
        
        if n < 10^4
            % general-purpose forw. diff. approx.
            gradf = @(x) (f(x + h(x)*eye(n)) - f(x))' / h(x);
        else
            % problem-dependent forw. diff. approx.
            gradf = @(x) [1:n]' .* (2*x + h(x));
        end

    case 'backw'
        % overwrite gradf with a f. handle that uses the backward difference
        % approximation
        h = @(x) (10^(-k_FDgrad))*norm(x);
        n = size(x0,1);
        
        if n < 10^4
            % general-purpose backw. diff. approx.
            gradf = @(x) (f(x) - f(x - h(x)*eye(n)))' / h(x);
        else
            % problem-dependent backw. diff. approx.
            gradf = @(x) [1:n]' .* (2*x - h(x));
        end

    case 'centr'
        % overwrite gradf with a f. handle that uses the centered difference
        % approximation
        h = @(x) (10^(-k_FDgrad))*norm(x);
        n = size(x0,1);

        if n < 10^4
            % general-purpose centr. diff. approx.
            gradf = @(x) (f(x + h(x)*eye(n)) - f(x - h(x)*eye(n)))' / (2*h(x));
        else
            % problem-dependent centr. diff. approx.
            gradf = @(x) [1:n]' .* (2*x);   % exact gradient!
        end

    otherwise
        % we use the input function handle gradf
end

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
projection_count = 0;

xk = Pi_X(x0); % Project the starting point if outside the constraints
if ~isequal(xk, x0)
    projection_count = projection_count + 1;
end
fk = f(xk);

% compute the gradient (exactly or approximately, depending on FDgrad)
gradfk = gradf(xk);

k = 0;
gradfk_norm = norm(gradfk);
deltaxk_norm = tolx + 1;    % this is to ensure that at least a first iteration is performed

if verbose
    disp(['iteration: ', num2str(k)])
    disp(['norm(xk): ', num2str(norm(xk))])
    disp(['gradient norm: ', num2str(gradfk_norm)])
    if strcmp(FDgrad,'forw') || strcmp(FDgrad,'backw') || strcmp(FDgrad,'centr')
        disp(['h(xk): ', num2str(h(xk))])
    end
end

while k < kmax && gradfk_norm >= tolgrad && deltaxk_norm >= tolx
    % Compute the descent direction (exactly or approximately, depending on FDgrad)
    pk = -gradf(xk);
    
    % Take a step in the descent direction and project the resulting vector
    % onto the feasible set X
    xhatk = xk + gamma * pk
    xbark = Pi_X(xhatk);
    if ~isequal(xbark, xhatk)
        projection_count = projection_count + 1;
    end
    

    % Backtracking line search

    % Reset the value of alpha
    alpha = 1;
    
    % Compute the candidate new xk
    pik = xbark - xk;
    xnew = xk + alpha * pik;
    
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    
    bt = 0;     % counter of backtracking iterations
    
    % 2nd condition is the Armijo (w.r.t. pik) condition not satisfied
    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pik)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pik;
        fnew = f(xnew);
        
        % Increase the counter by one
        bt = bt + 1;
        
    end
    
    % Update xk, fk, gradfk_norm, deltaxk_norm
    deltaxk_norm = norm(xnew - xk);
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);
    
    % Increase the step by one
    k = k + 1;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;

    if verbose && mod(k,10) == 0
        disp(['iteration: ', num2str(k)])
        disp(['norm(xk): ', num2str(norm(xk))])
        disp(['gradient norm: ', num2str(gradfk_norm)])
        if strcmp(FDgrad,'forw') || strcmp(FDgrad,'backw') || strcmp(FDgrad,'centr')
            disp(['h(xk): ', num2str(h(xk))])
        end
    end
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);

end
