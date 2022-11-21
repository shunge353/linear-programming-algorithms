function [xsol, fval, existflag, iters] = DantzigWolfeDecomp(extr, sub, n)
% Description: this function is an implementation of Dantzig-Wolfe decomposition
%              reference: 1961-George B.Dantzig and Phlip Wolfe-The decomposition algorithm for linear programs
%              url: https://www.jstor.org/stable/1911818
%
% Input:
% -- extr: struct that contains data of the extremal program, i.e. Aj an m by nj matrix, j = 1, ... ,n
%                                                                  cj an nj-vector, j = 1, ... ,n
%                                                                  b an m-vector
% -- sub: struct that contains data of the subproblems, i.e. Bj an mj by nj matrix, j = 1, ... ,n
%                                                            bj an mj-vector, j = 1, ... ,n
% -- n: the number of subproblems
%
% Output:
% -- xsol: optimal basic feasible solution of the decomposed program (size sum_j{nj} x 1)
% -- fval: optimal objective value of the decomposed program
% -- existflag: the reason that the algorithm terminate (0: optimal basic feasible solution found,
%               1: the decomposed program is infeasible, 2: the decomposed program is unbounded)
% -- iters: the number of iterations
%
% The decomposed program
%    minimize   c1 * x1 + c2 * x2 + ... + cn * xn
%    subject to A1 * x1 + A2 * x2 + ... + An * xn = b
%               B1 * x1                           = b1
%                         B2 * x2                 = b2
%                                   ...
%                                         Bn * xn = bn
%               x1, x2, ... ,xn >= 0
% where xj is a variable nj-vector, 
%       cj is an nj-vector
%       Aj is an m by nj matrix, b is an m-vector
%       Bj is an mj by nj matrix, bj is an mj-vector
% The number of variables sum_j{nj}, equality constraints m + sum_j{mj}.
%
% The extremal program
%    minimize   c11 * s11 + ... + c1K1 * s1K1 + c21 * s21 + ... + c2K2 * s2K2 + ...... + cn1 * sn1 + ... + CnKn * snKn
%    subject to P11 * s11 + ... + P1K1 * s1K1 + P21 * s21 + ... + P2K2 * s2K2 + ...... + Pn1 * sn1 + ... + PnKn * snKn = b  (pi)
%                 1 * s11 + ... +   1  * s1K1                                                                          = 1  (pi_bar_1)
%                                                 1 * s21 + ... +   1  * s2K2                                          = 1  (pi_bar_2)
%                                                                               ......
%                                                                                          1 * sn1 + ... +   1  * snKn = 1  (pi_bar_n)
%               s11, ... ,s1K1, s21, ... ,s2K2, ...... ,sn1, ... ,snKn >= 0
% where sjk is a variable Kj-vector for each j = 1, ... ,n or sjk is a variable scalar for each j = 1, ... ,n, k = 1, ... ,Kj,
%       Wj = {xj1, ... ,xjKj} is the set of all extreme points of the convex polyhedron Sj = { xj | Bj * xj = bj, xj >= 0 }
%       cjk = cj * xjk is a scalar for each j = 1, ... ,n, k = 1, ... ,Kj
%       Pjk = Aj * xjk is an m-vector for each j = 1, ... ,n, k = 1, ... ,Kj
%       b is an m-vector
%       pi is an m-vector simplex multipliers or prices
%       pi_bar is an n-vector simplex multipliers or prices
% The number of variables sum_j{Kj}, equality constraints m + n.
%
% The modified costs subproblem j
%    minimize   (cj - pi * Aj) * xj
%    subject to Bj * xj = bj
%               xj >= 0
% where xj is a variable nj-vector
%       cj is an nj-vector
%       pi is an m-vector, Aj is an m by nj matrix
%       Bj is an mj by nj matrix, bj is an mj vector
% The number of variables nj, equality constraints mj.
% There are three possibilities when solving the subproblem j in generating column:
%   (i): If there are some subproblem j is unbounded, the revised simplex method generate an extreme ray yj at termination,
%        which satisfy yj >= 0, Bj * yj = 0 and (cj - pi * Aj) * yj < 0. In this case, the column and cost to be added to
%        the extremal problem have the form (Aj * yj; 0) and cj * yj, in which 0 is an n-vector
%  (ii): If all subproblems are bounded, then for each j = 1, ... ,n let xj_bar be an optimal extreme point, let xk_bar be
%        such that \delta = (ck - pi * Ak) * xk_bar - pi_bar_k = Min_j {(cj - pi * Aj) * xj_bar - pi_bar_j), ( Dantzig's pivot rule )
%        if \delta < 0, form the new column and its associated cost for the extremal problem as (Ak * xk_bar; ek) and ck * xk_bar,
%        in which ek is a kth unit vector
% (iii): If all subproblems are bounded and if \delta >= 0 ( all reduced costs of nonbasic variables are nonnegative ),
%        terminate the algorithm, sjk solves the extremal problem, xj = sum_k{xjk * sjk} solves the decomposed problem.
%
% The auxiliary problem for the extremal program
%    minimize     0 * s11 + ... +   0  * s1K1 +  0  * s21 + ... +   0  * s2K2 + ...... +  0  * sn1 + ... +   0  * snKn + e * y
%    subject to P11 * s11 + ... + P1K1 * s1K1 + P21 * s21 + ... + P2K2 * s2K2 + ...... + Pn1 * sn1 + ... + PnKn * snKn + I * y = b  (pi)
%                 1 * s11 + ... +   1  * s1K1                                                                                  = 1  (pi_bar_1)
%                                                 1 * s21 + ... +   1  * s2K2                                                  = 1  (pi_bar_2)
%                                                                               ......
%                                                                                          1 * sn1 + ... +   1  * snKn         = 1  (pi_bar_n)
%               s11, ... ,s1K1, s21, ... ,s2K2, ...... ,sn1, ... ,snKn >= 0
%               y >= 0
% where sjk is an Kj-vector for each j = 1, ... ,n or sjk is a scalar for each j = 1, ... ,n, k = 1, ... ,Kj
%       y is an artificial m-vector
%       Wj = {xj1, ... ,xjKj} is the set of all extreme points of the convex polyhedron Sj = { xj | Bj * xj = bj, xj >= 0 }
%       Pjk = Aj * xjk is an m-vector for each j = 1, ... ,n, k = 1, ... ,Kj
%       b is an m-vector
%       e is an m-vector
%       I is an m by m identity matrix
%       pi is an m-vector simplex multipliers or prices
%       pi_bar is an n-vector simplex multipliers or prices
% The number of variables sum_j{Kj} + m, equality constraints m + n.
%
% The modified costs subproblem j of the auxiliary problem for the extremal program
%    minimize   (0 - pi * Aj) * xj
%    subject to Bj * xj = bj
%               xj >= 0
% where xj is a variable nj-vector
%       pi is an m-vector, Aj is an m by nj matrix
%       Bj is an mj by nj matrix, bj is an mj vector
% The number of variables nj, equality constraints mj.
% There are three possibilities when solving the subproblem j in generating column:
%   (i): If there are some subproblem j is unbounded, the revised simplex method generate an extreme ray yj at termination,
%        which satisfy yj >= 0, Bj * yj = 0 and (0 - pi * Aj) * yj < 0. In this case, the column and cost to be added to
%        the extremal problem have the form (Aj * yj; 0) and 0, in which the first 0 is an n-vector, the second 0 is a scalar
%  (ii): If all subproblems are bounded, then for each j = 1, ... ,n let xj_bar be an optimal extreme point, let xk_bar be
%        such that \delta = (0 - pi * Ak) * xk_bar - pi_bar_k = Min_j {(0 - pi * Aj) * xj_bar - pi_bar_j), ( Dantzig's pivot rule )
%        if \delta < 0, form the new column and its associated cost for the extremal problem as (Ak * xk_bar; ek) and 0,
%        in which ek is a kth unit vector
% (iii): If all subproblems are bounded and if \delta >= 0 && e - pi >= 0 ( all reduced costs of nonbasic variables are nonnegative ),
%        then the current basis is optimal, there are two possibilities:
%        (a): If e * y > 0, the linking constraint is violated, the decmposed program is infeasible, terminate the algorithm
%        (b): If e * y = 0, found a basic feasible solution for the extremal program, continue to Phase II
%
% The auxiliary problem for the subproblem j
%    minimize    0 * xj + e * yj
%    subject to Bj * xj + I * yj = bj
%               xj >= 0
%               yj >= 0
% where xj is an nj-vector
%       yj is an artificial mj-vector
%       e is an mj-vector
%       Bj is an mj by nj matrix
%       I is an mj by mj matix
%       bj is an mj-vector
% The number of variables nj + mj, equality constraints mj.
% There are two possibilities when solving the auxliary problem j:
%   (i): If e * yj > 0, the convex polyhedron Sj = { xj | Bj * xj = b, xj >= 0 } is empty, the decomposed program is infeasible,
%        terminate the algorithm
%  (ii): If e * yj = 0, found a basic feasible solution for the Sj, continue to Phase I

A1 = extr.A1; A2 = extr.A2; A3 = extr.A3;
c1 = extr.c1; c2 = extr.c2; c3 = extr.c3;
b = extr.b;
B1 = sub.B1; B2 = sub.B2; B3 = sub.B3;
b1 = sub.b1; b2 = sub.b2; b3 = sub.b3;
[m1, n1] = size(B1); [m2, n2] = size(B2); [m3, n3] = size(B3); % size of subproblems matrix Bj
m = size(b, 1); % number of the coupling constraints
tol = 1e-8; % tolerance for comparing with zero

% Initialize the output variables
xsol = [];
fval = [];
existflag = 0;
iters = 0;

% Iterate in Phase I of the subproblems
% ( find an extreme point x_{j1} for each of the subproblem's feasible region S_{j} (j = 1 ,..., K) )
disp('-------- PHASE I of the subproblems --------');
[x11, sub1] = sub_extr_point(B1, b1);
[x21, sub2] = sub_extr_point(B2, b2);
[x31, sub3] = sub_extr_point(B3, b3);

if ~isempty(find([sub1 sub2 sub3]))
    existflag = 1;
    disp('The decomposed program is infeasible');
    return;
end

% Iterate in Phase I of the extremal program
% ( find an initialbasic feasible solution for the extremal program )
disp('-------- PHASE I of the extremal program --------');
P11 = A1 * x11; P21 = A2 * x21; P31 = A3 * x31;

% Mutiplying both side of the coupling constraints by -1
% to make P11 + ... + Pn1 <= b
if ~isempty(find(b - P11 - P21 - P31 < 0))
    ind = find(b - P11 - P21 - P31 < 0);
    P11(ind) = -P11(ind); P21(ind) = -P21(ind); P31(ind) = -P31(ind);
    b(ind) = -b(ind);
end

cb = [zeros(n, 1); ones(m, 1)];         % coefficients of the objective function for the basic variables
bb = [b; ones(n, 1)];                   % vector of the right-hand side of the constraints
xb = [ones(n, 1); b - P11 - P21 - P31]; % basic feasible solution
Basis = [[P11 P21 P31] eye(m, m);       % matrix of the basic variables
    eye(n, n) zeros(n, m)];
BasisInv = inv(Basis);                  % basis inverse
w = cb' * BasisInv;                     % simplex multiplier
p = w(1 : m); p_bar = w(1 + m : end);

Basic_sub1 = "x11";
Basic_sub2 = "x21";
Basic_sub3 = "x31";
Basic = [Basic_sub1 Basic_sub2 Basic_sub3 "y1" "y2"];

while true
    % Solve all subproblems using prices from the extremal program
    [xsol1, ray1, fval1, sub1] = sub_modi_cost(B1, 0 - p * A1, b1);
    [xsol2, ray2, fval2, sub2] = sub_modi_cost(B2, 0 - p * A2, b2);
    [xsol3, ray3, fval3, sub3] = sub_modi_cost(B3, 0 - p * A3, b3);
    
    col_ray1 = []; col_ray2 = []; col_ray3 = [];
    col_ext1 = []; col_ext2 = []; col_ext3 = [];
    col_arv1 = []; col_arv2 = [];
    % If there are some subproblem unbounded
    if ~isempty(find([sub1 sub2 sub3]))
        if sub1 == 1
            col_ray1 = [A1 * ray1; zeros(n, 1)];
            cost_ray1 = 0;
            
            % Perform the minimum ratio test to find the leaving variable
            h_l = BasisInv * col_ray1; % pivot column
            mrt = find(h_l > 0);
            xb_div_hl = xb(mrt) ./ h_l(mrt);
            [~, rr] = min(xb_div_hl);
            r = mrt(rr);
            
            % Update Basis and prices
            Basic(r) = "y1";
            cb(r) = cost_ray1;
            Basis(:, r) = col_ray1;
            BasisInv = inv(Basis);
            xb = BasisInv * bb;
            w = cb' * BasisInv;
            p = w(1 : m); p_bar = w(1 + m : end);
            
        elseif sub2 == 1
            col_ray2 = [A2 * ray2; zeros(n, 1)];
            cost_ray2 = 0;
            
        elseif sub3 == 1
            col_ray3 = [A3 * ray3; zeros(n, 1)];
            cost_ray3 = 0;
            
            % perform the minimum ratio test to find the leaving variable
            h_l = BasisInv * col_ray3; % pivot column
            mrt = find(h_l > 0);
            xb_div_hl = xb(mrt) ./ h_l(mrt);
            [~, rr] = min(xb_div_hl);
            r = mrt(rr);
            
            % Update Basis and prices
            Basic(r) = "y1";
            cb(r) = cost_ray3;
            Basis(:, r) = col_ray3;
            BasisInv = inv(Basis);
            xb = BasisInv * bb;
            w = cb' * BasisInv;
            p = w(1 : m); p_bar = w(1 + m : end);
        end
    % If all subproblems bounded
    else
        [delta, t] = min([fval1 fval2 fval3] - p_bar);
        
        % If all reduced cost are nonnegative
        if delta >= -tol && all(1 - p >= -tol)
            % If the optimal cost of the auxiliary problem of the extremal program is positive
            if cb' * xb > 0
                existflag = 1;
                disp('The decomposed problem is infeasible');
                return;
            % else, found a basic feasible solution of the extremal program, turn to Phase II
            else
                break;
            end
        % else, pick a basic column to enteer the basis
        elseif delta < tol
            if t == 1
                col_ext1 = [A1 * xsol1; zeros(n, 1)];
                col_ext2(m + t) = 1;
                cost_ext1 = 0;
                
            elseif t == 2
                col_ext2 = [A2 * xsol2; zeros(n, 1)];
                col_ext2(m + t) = 1;
                cost_ext2 = 0;
                
                % perform the minimum ratio test to find the leaving variable
                h_l = BasisInv * col_ext2; % pivot column
                mrt = find(h_l > 0);
                xb_div_hl = xb(mrt) ./ h_l(mrt);
                [~, rr] = min(xb_div_hl);
                r = mrt(rr);
            
                % Update Basis and prices
                Basic(r) = "y1";
                cb(r) = cost_ext2;
                Basis(:, r) = col_ext2;
                BasisInv = inv(Basis);
                xb = BasisInv * bb;
                w = cb' * BasisInv;
                p = w(1 : m); p_bar = w(1 + m : end);
            
            elseif t == 3
                col_ext3 = [A3 * xsol3; zeros(n, 1)];
                col_ext3(m + t) = 1;
                cost_ext3 = 0;
            end
            
        elseif 1 - p < tol
            if 1 - p(1) < tol
                col_arv1 = zeros(m + n, 1);
                col_arv1(1) = 1;
                cost_arv1 = 1;
            elseif 1 - p(2) < tol
                col_arv2 = zeros(m + n, 1);
                col_arv2(2) = 1;
                cost_arv2 = 1;
            end
        end
    end
    
    % Column generation
    
end

% Iterate in Phase II of the extremal program
% ( solve the extremal program by reivsed simplex method )
disp('-------- PHASE II of the extremal program --------');
% generate simplex multipliers and solve the subproblems, iterates and so on


end

function [xsol, existflag] = sub_extr_point(A, b)
% Description: this function is an implementation of the Phase I method for the convex polyhedron
%              P = { x | A * x = b, x >= 0} to find a basic feasible solution
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- b: vector of right-hand side of the constraints (size m x 1)
%
% Output:
% -- xsol: basic feasible solution (size n x 1)
% -- existflag: the reason that the algorithm terminate (0: a basic feasible solution found,
%               1: the convex polyhedron is empty)

[m, n] = size(A); % size of the matrix A
tol = 1e-8;       % tolerance for comparing with zero

% Initialize the output variable
xsol = zeros(n, 1);
existflag = 0;

% Change the right-hand side vector b as nonnegative
if ~isempty(find(b < 0))
    ind = find(b < 0);
    b(ind) = -b(ind);
    A(ind, :) = -A(ind, :);
end

% Introduce artificial variables y1, ... ,ym
Basic = [n + (1 : m)];
NonBasic = setdiff([1 : (n + m)], Basic);
cb = ones(m, 1);      % coefficients of the objective function for the basic variables
xb = b;               % basic feasible solution

% Initialize the simplex tableau
tab = zeros(1 + m, 1 + m + n);
tab(1, 1) = -cb' * xb;
tab(1 + (1 : m), 1) = xb;
tab(1, 1 + (1 : n)) = -cb' * A;
tab(1 + (1 : m), 1 + (1 : n)) = A;
tab(1 + (1 : m), 1 + (n + 1 : n + m)) = eye(m);
cn = tab(1, 1 + (1 : n + m)); % reduced costs

while true
    % Optimality check
    if all(cn >= 0)
        % If the optimal cost associated with the auxiliary problem is positive
        if -tab(1, 1) >= tol
            xsol = [];
            existflag = 1;
            return;
        % else, the auxiliary problem has a zero optimal cost
        elseif abs(tab(1, 1)) <= tol
            [Lia, Locb] = ismember([n + (1 : m)], Basic);
            
            % If the artificial variable leaves from the final optimal basis, return a basic feasible solution
            if all(Lia == 0)
                xsol(Basic) = xb;
                existflag = 0;
                return;
            % else, thr artificial variable is in the final optimal basis and remain at zero level at optimum
            % then the artificial variable leaves the basis and pick a elegible variable as an entering one
            else
                % If there are nb artificial variables in the final optimal basis, pivot nb times
                curarv = 1;
                for i = 1:nnz(Locb)
                    rr = find(Locb, curarv);
                    r = rr(curarv);
                    
                    % If all the entries of the columns BasisInv * A{j}, j = 1, ... ,n is zero, the rth row is redundant
                    if all(tab(1 + r, 1 + (1 : n)) == 0)
                        tab(1 + r, :) = [];
                        Basic(r) = [];
                        xb = tab(2 : end, 1);
                    % else, find a nonzero element of BasisInv * A{j}, j = 1, ... ,n
                    else
                        [~, l] = find(tab(1 + r, 1 + (1 : n)), 1); % x{l} is a nobasic variable of
                        t = find(NonBasic == l);                   % the original problem that enters the basis
                        l = NonBasic(t);                           % index of the entering variable
                        k = Basic(r);                              % index of the leaving variable
                        
                        % Update basic and nonbasic list
                        Basic(r)    = l;
                        NonBasic(t) = k;
                        
                        % Update the basis inverse
                        pvt = tab(1 + r, 1 + l); % pivot element
                        nc = size(tab, 1);       % new size of the number of constraints
                        nonpvtrow = setdiff([1 : nc], 1 + r);
                        tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
                            repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n + m) .* repmat(tab(1 + r, :)', 1, nc - 1)';
                        tab(1 + r, :) = tab(1 + r, :) ./ pvt;
                    end
                    curarv = curarv + 1;
                end
                xsol(Basic) = xb;
                existflag = 0;
                return;
            end
        end
    % If the current basis is not optimal, choose a nonbasic variable to enter the basis   
    else
        [~, l] = min(cn);            % Dantzig's pivot rule
        t = find(NonBasic == l);
        l = NonBasic(t);             % index of the entering variable
        h_l = tab(2 : m + 1, 1 + l); % pivot column BasisInv * A{l}
        
        % Perform minimum ratio test to find the leaving variable
        mrt = find(h_l > 0);
        xb_div_hl = xb(mrt) ./ h_l(mrt);
        p = min(xb_div_hl);
        rr = find(xb_div_hl == p);
        
        % Break the ties in order to aviod stalling
        if length(rr) > 1
            r = mrt(rr);
            k = Basic(r); % index of the leaving variable
            k = min(k);
            r = find(Basic == k);
        else
            r = mrt(rr);
            k = Basic(r); % index of the leaving variable
        end
        
        % Update basic and nonbasic list
        Basic(r)    = l;
        NonBasic(t) = k;
        
        % Update the pivot column BasisInv * A{l}
        pvt = tab(1 + r, 1 + l); % pivot element
        nonpvtrow = setdiff([1 : 1 + m], 1 + r);
        tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
            repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n + m) .* repmat(tab(1 + r, :)', 1, m)';
        tab(1 + r, :) = tab(1 + r, :) ./ pvt;
        
        xb = tab(2 : 1 + m, 1);       % basic feasible solution
        cn = tab(1, 1 + (1 : n + m)); % reduced costs
    end
end
end

function [xsol, ray, fval, existflag] = sub_modi_cost(A, c, b)
% Description: this function is an implementation of the big-M simplex method to solve the modified costs subproblem
%
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the objective function (size 1 x n)
% -- b: vector of right-hand side of the constraints (size m x 1)
%
% Output:
% -- xsol: optimal basic feasible solution (size n x 1)
% -- ray: extreme ray when the problem is unbounded (size n x 1)
% -- fval: optimal objective value
% -- existflag: the reason that the algorithm terminate (0: optimal basic feasible solution found,
%               1: the problem is unbounded)

[m, n] = size(A); % size of the matrix A
M = 100;          % a sufficient large positive constant
tol = 1e-8;       % tolerance for comparing with zero

% Initialize the output variables
xsol = zeros(m + n, 1);
ray = zeros(m + n, 1);
fval = [];
existflag = 0;

% Change the right-hand side vectors to be nonnegative
if ~isempty(find(b < 0))
    ind = find(b < 0);
    b(ind) = -b(ind);
    A(ind, :) = -A(ind, :);
end

% Introduce the artificial variables y1, ... ,ym
Basic = [n + (1 : m)];                    % list of basic variables
NonBasic = setdiff([1 : (n + m)], Basic); % list of nonbasic variables
cb = M * ones(m, 1); % coefficients of the objective function for the basic variables associated with the big-M problem
xb = b;              % basic feasible solution associated with the big-M problem

% Initialize the simplex tableau
tab = zeros(1 + m, 1 + n + m);
tab(1, 1) = -cb' * xb;
tab(1 + (1 : m), 1) = xb;
tab(1, 1 + (1 : n)) = c - cb' * A;
tab(1 + (1 : m), 1 + (1 : n)) = A;
tab(1 + (1 : m), 1 + (n + 1 : n + m)) = eye(m);
cn = tab(1, 1 + (1 : n + m)); % reduced costs

% Iterate in big-M
while true
    % Optimality check
    if all(cn >= -tol)
        xsol(Basic) = xb;
        xarv = xsol([n + (1 : m)]);
        % If the simplex method terminate with a solution with all the artificial variables equal to zero
        if all(abs(xarv) <= tol)
            xsol = xsol([1 : n]);
            ray = [];
            fval = -tab(1, 1);
            existflag = 0;
            return;
        end
    % If the current basis is not optimal, choose another variable to enter the basis
    else
        [~, l] = min(cn); % Dantzig's pivot rule
        t = find(NonBasic == l);
        l = NonBasic(t); % index of the entering variable
        h_l = tab(2 : m + 1, 1 + l); % pivot column BasisInv * A_{l}
        
        mrt = find(h_l > 0);
        % Perform the minimum ratio test to find the leaving variable
        if ~isempty(mrt)
            xb_div_hl = xb(mrt) ./ h_l(mrt);
            p = min(xb_div_hl);
            rr = find(xb_div_hl == p);
            
            % Break the tie in order to avoid stalling
            if length(rr) > 1
                r = mrt(rr);
                k = Basic(r); % index of the leaving variable
                k = min(k);
                r = find(Basic == k);
            else
                r = mrt(rr);
                k = Basic(r); % index of the leaving variable
            end
            
            % Update basic and nonbasic list
            Basic(r) = l;
            NonBasic(t) = k;
            
            % Update the pivot column BasisInv * A_{l}
            pvt = tab(1 + r, 1 + l); % pivot element
            nonpvtrow = setdiff([1 : 1 + m], 1 + r);
            tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
                repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n + m) .* repmat(tab(1 + r, :)', 1, m)';
            tab(1 + r, :) = tab(1 + r, :) ./ pvt;
            
            xb = tab(2 : 1 + m, 1);       % basic feasible solution
            cn = tab(1, 1 + (1 : n + m)); % reduced costs
        % No component of pivot column is positive and all artificial variables are at zero level
        % (since all subproblems feasible region is nonempty, it could only happen that the modified cost subproblem is unbounded)
        else
            xsol = [];
            ray(Basic) = -h_l;
            ray(l) = 1;
            ray = ray([1 : n]);
            fval = -Inf;
            existflag = 1;
            return;
        end
    end
end
end