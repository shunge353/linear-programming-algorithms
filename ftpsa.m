function [xsol, fval, existflag, iters] = ftpsa(A, c, b)
% Description: the function is an implementation of two-phase full tableau simplex method
%              reference: page 116 of D.Bersekas and J.N.Tsitsiklis ¡¶Introduction to Linear Optimization¡·
%              url: http://www.athenasc.com/linoptbook.html
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the obejctive function (size n x 1)
% -- b: vector of right-hand side of the constraints (size m x 1)
%
% Output:
% -- xsol: optimal basic feasible solution (size n x 1)
% -- fval: optimal objective value
% -- existflag: the reason that the algorithm terminate (0: optimal basic feasible solution found,
%               1: the original problem is infeasible, 2: the original problem is unbounded)
% -- iters: the number of iterations

[m, n] = size(A); % size of the matrix A

% Initialize the output variables
xsol = zeros(n, 1);
fval = 0;
existflag = 0;
iters = 0;

% Change the right-hand side vectors b as nonnegative
if ~isempty(find(b < 0))
    ind = find(b < 0);
    b(ind) = -b(ind);
    A(ind, :) = -A(ind, :);
end

% Introduce artificial variables y1, ... ,ym
Basic = [n + (1 : m)];
NonBasic = setdiff([1 : (n + m)], Basic);
fb = ones(m, 1); % coefficients of the objective function for the basic variables associated with the auxiliary problem
xb = b;          % basic feasible solution associated with the auxiliary problem

% Initialize simplex tableau
tab = zeros(1 + m, 1 + n + m);
tab(1, 1) = -fb' * xb;
tab(1 + (1 : m), 1) = xb;
tab(1, 1 + (1 : n)) = -fb' * A;
tab(1 + (1 : m), 1 + (1 : n)) = A;
tab(1 + (1 : m), 1 + (n + 1 : n+m)) = eye(m);
cn = tab(1, 1 + (1 : n + m)); % reduced costs

% Iterate in Phase I
disp('-------- PHASE I --------');
while true
    % Optimality check
    if all(cn >= 0)
        % If the optimal cost associated with the auxiliary problem is positive
        if -tab(1,1) >= eps
            disp('The LP problem is infeasible');
            existflag = 1;
            return;
        % the auxiliary problem has a zero optimal cost
        elseif abs(tab(1,1)) <= eps
            [Lia, Locb] = ismember([n + (1 : m)], Basic);
            
            % If the artificial variable is leaving from the final optimal basis, Phase I terminated and Phase II start
            if all(Lia == 0)
                tab(:, 1 + (n + 1 : n + m)) = [];   % the columns asociated with artificial variables are eliminated
                NonBasic = setdiff([1 : n], Basic); % non basic variable associate with the original problem
                iters = iters + 1;
                break;
            % the artificial variable is in the final optimal basis and remain at zero level at optimum
            % then the artificial variable is leaving the basis and pick a elegible original variable as an entering one
            else
                % If there are nb artificial variables in the final optimal basis, pivot nb times
                curarv = 1; % the current artificial variable in the final optimal basis
                for i = 1:nnz(Locb)
                    rr = find(Locb, curarv);
                    r = rr(curarv);
                    
                    % If all the entries of the columns BasisInv * A_{j}, j = 1, ... ,n is zero , the rth row is redundant
                    if all(tab(1 + r, 1 + (1 : n)) == 0)
                        tab(1 + r, :) = [];
                        Basic(r) = [];
                        iters = iters + 1;
                    % find a nonzero element of A_{j}, j = 1, ... ,n,
                    else
                        [~, l] = find(tab(1 + r, 1 + (1 : n)), 1); % x_{l} is a nonbasic variable of
                        t = find(NonBasic == l);                   % original problem which enters the basis
                        l = NonBasicList(t);                       % index of the entering variable
                        k = Basic(r);                              % index of the leaving variable
                        
                        % Update the basic and nonbasic list
                        Basic(r)    = l;
                        NonBasic(t) = k;
                        
                        % Update the pivot column BasisInv * A_{l}
                        pvt  = tab(1 + r, 1 + l); % pivot element
                        nc = size(tab, 1);        % new size of the number of constraints
                        % tab(i, :) = tab(i, :) - (tab(i, 1 + l) / pvt) * tab(1 + r, :); % i ~= 1 + r, to be vectorized later ...
                        nonpvtrow = setdiff([1 : nc], 1 + r);
                        tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
                            repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n + m) .* repmat(tab(1 + r, :)', 1, nc - 1)';
                        tab(1 + r, :) = tab(1 + r, :) ./ pvt;
                        iters = iters + 1;
                    end
                    curarv = curarv + 1;
                end
                tab(:, 1 + (n + 1 : n + m)) = [];   % the columns asociated with artificial variables are eliminated
                NonBasic = setdiff([1 : n], Basic); % non basic variable associate with the original problem
                break;
            end
        end
    % If the current basis is not optimal
    % choose a nonbasic variable to enter the basis
    else
        [~, l] = min(cn);            % Dantzig's pivot rule
        t = find(NonBasic == l);
        l = NonBasic(t);             % index of the entering variable
        h_l = tab(2 : m + 1, 1 + l); % pivot column BasisInv * A_{l}
        
        % Perform the minimum ratio test to find the leaving variable
        mrt = find(h_l > 0);
        xb_div_hl = xb(mrt) ./ h_l(mrt);
        p = min(xb_div_hl);
        rr = find(xb_div_hl == p);
        
        % Break the ties in order to avoid stalling
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
        
        % Update the pivot column BasisInv * A_{l}
        pvt = tab(1 + r, 1 + l); % pivot element
        nonpvtrow = setdiff([1 : 1 + m], 1 + r);
        tab(nonpvtrow, :) = tab(nonpvtrow, :) - ...
            repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n + m) .* repmat(tab(1 + r, :)', 1, m)';
        tab(1 + r, :) = tab(1 + r, :) ./ pvt;
        
        xb = tab(2 : 1 + m, 1);       % basic feasible solution
        cn = tab(1, 1 + (1 : n + m)); % reduced costs
        iters = iters + 1;
    end
end

% Recompute the zeroth row of tableau
m = size(tab, 1) - 1; % the number of constraints of the new tableau
xb = tab(1 + (1 : m), 1);
cb = c(Basic); % coefficients of the objective function for the basic variables
cn = c' - cb' * tab(1 + (1 : m), 1 + (1 : n)); % reduced costs
tab(1, 1) = -cb' * xb;
tab(1, 1 + (1 : n)) = cn;

% Iterate in Phase II
disp('-------- PHASE II --------');
while true
    % Optimality check
    if all(cn >= 0)
        % The original LP problem is optimal
        disp('The LP problem is optimal');
        existflag = 0;
        xsol(Basic) = xb;
        fval = -tab(1,1);
        return;
    % If the current basi is not optimal
    % choose another variable to the the basis
    else
        [~, l] = min(cn); % Dantzig's pivot rule
        t = find(NonBasic == l);
        l = NonBasic(t);             % index of the entering variable
        h_l = tab(2 : m + 1, 1 + l); % pivot column BasisInv * A_{l}
        
        mrt = find(h_l > 0);
        if ~isempty(mrt)
            % Perform the minimum ratio test to find the leaving variable
            xb_div_hl = xb(mrt) ./ h_l(mrt);
            p = min(xb_div_hl);
            rr = find(xb_div_hl == p);
            
            % Break the ties in order to avoid stalling
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
                repmat(tab(nonpvtrow, 1 + l) ./ pvt, 1, 1 + n) .* repmat(tab(1 + r, :)', 1, m)';
            tab(1 + r, :) = tab(1 + r, :) ./ pvt;
            
            xb = tab(2 : 1 + m, 1);   % basic feasible solution
            cn = tab(1, 1 + (1 : n)); % reduced costs
            iters = iters + 1;
        else
            % No component of pivot column is positive
            disp('The LP problem is unbounded');
            existflag = 2;
            iters = iters + 1;
            return;
        end
    end
end
end