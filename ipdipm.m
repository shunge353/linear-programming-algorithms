function [xsol, fval, existflag, iters] = ipdipm(A, c, b, Eqin, MinMaxLP, c0, maxIters, tol, etaMin)
% Description: the function is an implementation of
% Mehrotra's Predictor-Corrector infeasible primal-dual interior-point method
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the onjective function (size n x 1)
% -- b: vector of the right-hand side of the constraints (size m x 1)
% -- Eqin: vector of the type of the constraints (size m x 1)
% -- MinMaxLP: type of the optimization (optional: default value -1 - minimization)
% -- c0: constant term of the objective function (optional: default value 0)
% -- maxIters: maximum number of iterations to perform if an optimal solution if not found
%    (optional: default value 100)
% -- tol: tolerance for the termination criterion (optional: defalut value 1e-08)
% etaMin: parameter rta (optiona: default value 0.995)
%
% Output:
% -- xsol: the solution found (size n x 1)
% -- fval: the value of the objective function at the solution xsol
% -- existflag: the reason that the algorithm terminated
%    (0: optimal solution found, 1: the LP problem is infeasible)
% -- iters: the number of iterations

% initialize output variable
xsol = [];
fval = [];
existflag = 0;
iters = 0;

[m, n] = size(A);

% set default values to missing inputs
if ~exist('MinMaxLP')
    MinMaxLP = -1;
end
if ~exist('c0')
    c0 = 0;
end
if ~exist('maxIters');
    maxIters = 100;
end
if ~exist('tol')
    tol = 1e-08;
end
if ~exist('etaMin')
    etaMin = 0.995;
end

% if the type of optimization is maximization, then multiply vector c and constant c0 by -1
if MinMaxLP == 1
    c = -c;
    c0 = -c0;
end

% transform a general LP to its standard form
if ~isempty(find(Eqin ~= 0))
    slack = nnz(Eqin);
    c(n + (1:slack)) = 0;
    A(:, n + (1:slack)) = zeros(m, slack);
    curcol = 1;
    for i = 1:m
        % 'greater than or equal to' inequality constraint
        if Eqin(i) == 1
            A(i, n + curcol) = -1;
            curcol = curcol + 1;
        % 'less than or equal to' inequality constraint
        elseif Eqin(i) == -1
            A(i, n + curcol) = 1;
            curcol = curcol + 1;
        end
    end
end

[m, n] = size(A); % new size of the matrix A

% check the consistency of linear systems Ax = b
if rank(A) ~= rank([A b])
    existflag = 1;
    return;
end

% tolerance for detecting infeasible LP
infTol = 1e+15 * max([norm(A) norm(b) norm(c)]);

% calculate the starting point using Mehrotra's heuristic
x = A' * ((A * A') \ b);
w = (A * A') \ (A * c);
s = c - A' * w;
delta_x = max(-1.5*min(x), 0);
delta_s = max(-1.5*min(s), 0);
e = ones(n, 1);
delta_x_s = 0.5 * (x + delta_x * e)' * (s + delta_s * e);
delta_x_c = delta_x + delta_x_s / (sum(s) + n * delta_s);
delta_s_c = delta_s + delta_x_s / (sum(x) + n * delta_x);
x = x + delta_x_c * e;
s = s + delta_s_c * e;

% iterate until maxIters or the LP problem is optimal or infeasible
for i = 1:maxIters
    iters = iters + 1;
    % calculate the residuals and the duality measures
    rp = A * x - b; % primal residual
    rd = A' * w + s - c; % dual residual
    rc = x .* s; % complementarity residual
    mu = (x' * s) / n; % duality measure
    residual = norm(rc, 1) / (1 + abs(b' * w));
    % heuristic to detect infeasible LP
    if ( isnan(residual) ) || ( norm(x) + norm(s) >= infTol )
        existflag = 1;
        disp('The LP problem is infeasible');
        return;
    end
    % termination criterion
    if max([mu, norm(rp), norm(rd)]) <= tol
        xsol = x;
        if MinMaxLP == 1
            fval = -(c' * x + c0);
        else
            fval = c' * x + c0;
        end
        existflag = 0;
        disp('The LP problem is optimal');
        return;
    end
    % formulate the coeffiient matrix of the linear systems
    M = A * diag(x ./ s) * A';
    [R, p] = chol(M);
    % if the matrix is not positive definite, we assume that the LP problem is optimal
    if p > 0
        xsol = x;
        if MinMaxLP == 1
            fval = -(c' * x + c0);
        else
            fval = c' * x + c0;
        end
        existflag = 0;
        diap('The LP problem is optimal');
        return;
    end
    % predictor step
    rhs = rp - A * ((rc - x .* rd) ./ s);
    dw = R \ (R' \ rhs);
    ds = rd - A' * dw;
    dx = (rc - x .* ds) ./ s;
    ind_x = find(dx > 0);
    ind_s = find(ds > 0);
    alpha_p = min(1, min(x(ind_x) ./ dx(ind_x)));
    alpha_d = min(1, min(s(ind_s) ./ ds(ind_s)));
    % centering parameter step
    sigma = ((x - alpha_p * dx)' * (s - alpha_d * ds) / (n * mu)) ^3;
    % corrector step
    rc = rc - sigma * mu * e + dx .* ds;
    rhs = rp - A * ((rc - x .* rd) ./ s);
    dw = R \ (R' \ rhs);
    ds = rd - A' * dw;
    dx = (rc - x .* ds) ./ s;
    eta = max(etaMin, 1 - mu);
    ind_x = find(dx > 0);
    ind_s = find(ds > 0);
    alpha_p = min(1, eta * min(x(ind_x) ./ dx(ind_x)));
    alpha_d = min(1, eta * min(s(ind_s) ./ ds(ind_s)));
    % update step
    x = x - alpha_p * dx;
    w = w - alpha_d * dw;
    s = s - alpha_d * ds;
end
end