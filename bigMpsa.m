function [xsol, fval, existflag, iters] = bigMpsa(A, c, b)
% Description: the function is an implementation of big-M full tableau simplex method
%              reference: page 117 and page 135 of D.Bersekas and J.N.Tsitsiklis ¡¶Introduction to Linear Optimization¡·
%              url: http://www.athenasc.com/linoptbook.html
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the objective function (size n x 1)
% -- b: vector of right-hand side of the constraints (size m x 1)
%
% Output:
% -- xsol: optimal basic feasible solution (size n x 1)
% -- fval: optimal objective value
% -- existflag: the reason that the algorithm terminate (0: optimal basic feasible solution found,
%               1: the original problem is infeasible, 2: the original problem is unbounded)
% -- iters: the number of iterations

[m, n] = size(A); % size of the matrix A
M = 100;          % a sufficient large positive constant
tol = 1e-8;       % tolerance for comparing with zero

% Initialize the output variables
xsol = zeros(n + m, 1);
fval = [];
existflag = 0;
iters = 0;

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
tab(1, 1 + (1 : n)) = c' - cb' * A;
tab(1 + (1 : m), 1 + (1 : n)) = A;
tab(1 + (1 : m), 1 + (n + 1 : n + m)) = eye(m);
cn = tab(1, 1 + (1 : n + m)); % reduced costs

% Iterate in big-M
disp('-------- big-M --------');
while true
    % Optimality check
    if all(cn >= 0)
        xsol(Basic) = xb;
        xarv = xsol([n + (1 : m)]);
        % If the simplex method terminate with a solution with all the artificial variables equal to zero
        if all(abs(xarv) <= tol)
            xsol = xsol([1 : n]);
            fval = -tab(1, 1);
            existflag = 0;
            disp('The LP problem is optimal');
            return;
        % else, some of the artificial variables are nonzero at optimal basic feasible solution
        else
            xsol = [];
            fval = [];
            existflag = 1;
            disp('The LP problem is infeasible');
            return;
        end
    % If the current basis is not optimal
    % choose another variable to enter the basis
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
            iters = iters + 1;
        % No component of pivot column is positive
        else
            xsol(Basic) = xb;
            xarv = xsol([n + (1 : m)]);
            % If all the artificial variables are equal to zero
            if all(abs(xarv) <= tol)
                xsol = [];
                fval = -Inf;
                existflag = 2;
                iters = iters + 1;
                disp('The LP problem is unbounded');
                return;
            % else, some of the artificial variables are nonzero
            else
                xsol = [];
                fval = [];
                existflag = 1;
                iters = iters + 1;
                disp('The LP problem is infeasible');
                return;
            end
        end
    end
end
end