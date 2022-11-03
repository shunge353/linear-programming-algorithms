function [xsol, fval, existflag, iters, ray] = ftdsa(A, c, b)
% Description: this function is an implementation of a two-phase full tableau dual simplex method
%              reference: page 159 of D.Bersekas and J.N.Tsitsiklis ¡¶Introduction to Linear Optimization¡·
%              url: http://www.athenasc.com/linoptbook.html
%              reference: page 72 of SC Fang and S Puthenpura ¡¶Linear Optimization and Extensions¡·
%              url: https://dl.acm.org/doi/abs/10.5555/153367
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the objective function (size n x 1)
% -- b: vector of right-hand side of the constraints (size m x 1)
%
% Output:
% -- xsol: optimal basic feasible solution (size n x 1)
% -- fval: optimal objective value
% -- existflag: the reason that the algorithm terminate (0: optimal primal basic feasible solution found,
%               1: the primal problem is infeasible, 2: the primal problem is unbounded)
% -- iters: the number of iterations
% -- ray: extreme ray when primal problem unbounded

[m, n] = size(A); % size of the matrix A

% Initialize the output variables
xsol = zeros(n, 1);
fval = [];
existflag = 0;
iters = 0;
ray = zeros(n, 1);

% Find an initial basis
disp('-------- INITIAL BASIS --------');
[~, p] = rref([A b]);
Basic = sort(p);                  % list of basic variables
NonBasic = setdiff([1:n], Basic); % list of nobasic variables
Basis = A(:, Basic);              % matrix of the basic variables
cb = c(Basic);                    % coefficients of the basic variables
cn = c' - cb' * (Basis \ A);      % reduced costs

% Check the feasibility of basic solution to the dual problem
if all(cn >= 0)
    flag = 1;
else
    flag = 0;
end

% Iterate in Phase I
if flag == 0
    disp('-------- PHASE I --------');
    % Initialize the simplex tableau
    bM = Basis * ones(m, 1);
    xb = Basis \ bM;
    tab = zeros(1 + m, 1 + n);
    tab(1, 1) = -cb' * xb;
    tab(1 + (1 : m), 1) = xb;
    tab(1, 1 + (1 : n)) = cn;
    tab(1 + (1 : m), 1 + (1 : n)) = Basis \ A;
    
    % Iterate in full tableau primal simplex method
    while true
        % If the current basis is optimal
        % it is associated with a basic feasible solution to the dual
        if all(cn >= 0)
            flag = 1;
            break;
        % If the current basis is not optimal
        % choose a nonbasic variable to enter the basis
        else
            [~, l] = min(cn);              % Dantzig's pivot rule
            t = find(NonBasic == l);
            l = NonBasic(t);               % index of the entering variable
            h_l = tab(1 + (1 : m), 1 + l); % pivot column BasisInv * A_{l}
            
            mrt = find(h_l > 0);
            if ~isempty(mrt)
                % Perform the minimum ratio test to find the leaving variable
                xb_div_hl = xb(mrt) ./ h_l(mrt);
                p = min(xb_div_hl);
                rr = find(xb_div_hl == p);
                
                % Break the tie in order to aviod stalling
                if length(rr) > 1
                    r = mrt(rr);
                    k = Basic(r); % index of the exiting variable
                    k = min(k);
                    r = find(Basic == k);
                else
                    r = mrt(rr);
                    k = Basic(r); % index of the exiting variable
                end
                
                % Update the basic and nonbasic list
                Basic(r) = l;
                NonBasic(t) = k;
                
                % Update the basic column BasisInv * A_{l}
                pvt = tab(1 + r, 1 + l); % pivot element
                nonpvtrow = setdiff([1 : 1 + m], 1 + r);
                tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
                    repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n) .* repmat(tab(1 + r, :)', 1, m)';
                tab(1 + r, :) = tab(1 + r, :) ./ pvt;
                
                xb = tab(1 + (1 : m), 1); % basic feasible solution
                cn = tab(1, 1 + (1 : n)); % reduced costs
                iters = iters + 1;
            else
                % No component of pivot column is positive
                disp('The LP problem is unbounded');
                existflag = 2;
                iters = iters + 1;
                ray(Basic) = -h_l;
                ray(l) = 1;
                return;
            end
        end
    end
end

% Iterate in Phase II
if flag == 1
    disp('-------- PHASE II --------');
    % Recompute the simplex tableau
    Basis = A(:, Basic);         % matrix of the basic variables
    cb = c(Basic);               % coefficient of the basic variables
    cn = c' - cb' * (Basis \ A); % reduced costs
    xb = Basis \ b;
    tab(1, 1) = -cb' * xb;
    tab(1 + (1 : m), 1) = xb;
    tab(1, 1 + (1 : n)) = cn;
    tab(1 + (1 : m), 1 + (1 : n)) = Basis \ A;
    
    % Iterate in full tablear dual simplex method
    while true
        % If the current basis is optimal
        if all(xb >= 0)
            disp('The LP problem is optimal');
            xsol(Basic) = xb;
            fval = -tab(1, 1);
            existflag = 0;
            return;
        % If the current basis is not optimal
        % choose a nonbasic variable to enter the basis
        else
            [~, r] = min(xb);
            k = Basic(r);                  % index of the exiting variable
            v_r = tab(1 + r, 1 + (1 : n)); % pivot row rth element of BasisInv * A_{i}, i = 1, ... ,n
            
            mrt = find(v_r < 0);
            if ~isempty(mrt)
                % Perform minimum ratio test to find the entering variable
                cn_div_vr = cn(mrt) ./ abs(v_r(mrt));
                p = min(cn_div_vr);
                ll = find(cn_div_vr == p);
                
                % Break the tie in order to aviod stalling
                if length(ll) > 1
                    l = mrt(ll);
                    t = find(NonBasic == l); % index of the entering variable
                    t = min(t);
                    l = NonBasic(t);
                else
                    l = mrt(ll);
                    t = find(NonBasic == l); % index of the entering variable
                    l = NonBasic(t);
                end
                
                % Update the basic and nonbasic list
                Basic(r) = l;
                NonBasic(t) = k;
                
                % Update the basic column BasicInv * A_{l}
                pvt = tab(1 + r, 1 + l); % pivot element
                nonpvtrow = setdiff([1 : 1 + m], 1 + r);
                tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
                    repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n) .* repmat(tab(1 + r, :)', 1, m)';
                tab(1 + r, :) = tab(1 + r, :) ./ pvt;
                
                xb = tab(1 + (1 : m), 1); % basic feasible solution
                cn = tab(1, 1 + (1 : n)); % reduced costs
                iters = iters + 1;
            else
                % No component of the pivot row is negative
                disp('The LP problem is infeasible');
                existflag = 1;
                iters = iters + 1;
                return;
            end
        end
    end
end
end