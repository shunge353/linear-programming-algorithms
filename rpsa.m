function [xsol, fval, existflag, iters] = rpsa(A, c, b, Eqin, MinMaxLP, c0, reinv)
% Description: this function is an implementation of the revised primal simplex algorithm
%
% Input:
% -- A: matrix of coefficients of the constraints (size m x n)
% -- c: vector of coefficients of the objective function (size n x 1)
% -- b: vector of the right-hand side of the constraints (size m x 1)
% -- Eqin:vectro of the type of the constraints (size m x 1)
% -- MinMaxLP: the type of optimization (optional: default value -1 - minimization)
% -- c0: constant term of the objective function (optional: default value 0)
% -- reinv: every reinv number of iterations, the basis inverse is re-computed from scratch
%    (optional: default value 80)
%
% Output:
% -- xsol: the solution found (size m x 1)
% -- fval: the value of the objective function at the solution x
% -- existflag: the reason that the algorithm terminated (0: optimal solution found,
%    1: the LP problem is infeasible, 2: the LP problem is unbounded)
% -- iters: the number of iterations

% initialize output variables
xsol = [];
fval = [];
existflag = 0;
iters = 0;

[m, n] = size(A); % size of the matrix A
counter = 1; % number of re-computation of basis inverse

% set default value to missing inouts
if ~exist('MinMaxLP')
    MinMaxLP = -1;
end
if ~exist('c0')
    c0 = 0;
end
if ~exist('reinv')
    reinv = 80;
end

% if the type of optimization is maximization, then multiply vector c and constant c0 by -1
if MinMaxLP == 1
    c = -c;
    c0 = -c0;
end

% transform a general LP problem to its standard form
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

% check the consistency of linear systems Ax == b
if rank(A) ~= rank([A b])
    existflag = 1;
    return;
end

disp('---- I N I T I A L    B A S I S ----');

% find an invertible basis using rref
[R, p] = rref([A b]);
BasicList = sort(p);
NonBasicList = setdiff([1:n], BasicList);

% delete redundant constraints if exists
idx = all(R, 2);
A(idx, :) = 0;
b(idx) = 0;

[m, n] = size(A); % new size of the matrix A

Basis = A(:, BasicList);       % matrix of the basic variable
BasisInv = inv(Basis);         % basis inverse
xb = BasisInv * b;             % basic solution

% check if the Phase I is needed
if all(xb >= 0)
    flag = 1;
else
    flag = 0;
end

% Phase I
if flag == 0
    disp('---- P H A S E    I ----');
    % find the index of the minimum element of vector xb
    [~, minidx] = min(xb);
    % create a copy of vector c for Phase I
    c2 = zeros(n, 1);
    c2(n + 1) = 1;
    % create a vector d
    d = -Basis * ones(m, 1);
    % re-define matrix A since the artificial variable needs to be added for Phase I
    A(:, n + 1) = d;
    % compute the new basic and nonbasic variables
    NonBasicList(1, length(NonBasicList) + 1) = BasicList(1, minidx);
    BasicList(1, minidx) = n + 1;
    
    [m, n] = size(A); % new size of matrix A
    
    Basis = A(:, BasicList);       % matrix of the new basic variable
    fb = c2(BasicList);            % coefficients of the objective function for the new basic variables
    NonBasis = A(:, NonBasicList); % matrix of the new nonbasic variables
    fn = c2(NonBasicList);         % coefficients of the objective function for the new nonbasic variables
    BasisInv = inv(Basis);         % basis inverse
    xb = BasisInv * b;             % basic solution
    w = fb' * BasisInv;            % simplex multiplier
    sn = fn' - w * NonBasis;       % reduced costs
    
    % iterate in Phase I
    while true
        % optimality test for Phase I
        if all(sn >= 0)
            % if the artificial variable is leaving from the basic list, Phase I terminated and Phase II start
            [~, col] = find(NonBasicList == n);
            if ~isempty(col)
                NonBasicList(col) = [];
                A(:, n) = [];
                flag = 1;
                iters = iters + 1;
                break;
            else
            % if the artificial variable is in the basic list but reduced costs sn >= 0, and if it is zero at optimum,
            % then the artificial variable is leaving basic list and pick a random but acceptable variable as an entering one
            % The basis produced is feasible for Phase II
                [~, col] = find(BasicList == n);
                if ~isempty(col)
                    if xb(col) == 0
                        A(:, n) = [];
                        for i = 1:length(NonBasicList)
                            l = NonBasicList(i);
                            h_l = BasisInv * A(:, l); % pivot column
                            if h_l(col) ~= 0
                                % update the basic and nonbasic lists
                                BasicList(col) = l;
                                NonBasicList(i) = [];
                                flag = 1;
                                iters = iters + 1;
                                break;
                            end
                        end
                        break;
                    else
                    % the LP problem is infeasible
                        disp('The LP problem is infeasible');
                        existflag = 1;
                        iters = iters + 1;
                        return;
                    end
                end
            end
        else
        % the LP problem of Phase I is not optimal
        % find the entering variable using the selected pivoting rule
            [~, t] = min(sn);         % Dantzig
            l = NonBasicList(t);      % index of the entering variable
            h_l = BasisInv * A(:, l); % pivot column
            % perform the minimum ratio test to find the leaving variable
            mrt = find(h_l > 0);
            xb_div_hl = xb(mrt) ./ h_l(mrt);
            [p, r] = min(xb_div_hl);
            rr = find(xb_div_hl == p);
            % break the ties in order to avoid stalling
            if length(rr) > 1
                r = mrt(rr);
                k = BasicList(r); % index of the leaving variable
                k = max(k);
                r = find(BasicList == k);
            else
                r = mrt(r);
                k = BasicList(r); % index of the leaving variable
            end
            % update basic and nonbasic list
            BasicList(r) = l;
            NonBasicList(t) = k;
            
            Basis = A(:, BasicList);       % matrix of the new basic variable
            fb = c2(BasicList);            % coefficients of the objective function for the new basic variables
            NonBasis = A(:, NonBasicList); % matrix of the new nonbasic variables
            fn = c2(NonBasicList);         % coefficients of the objective function for the new nonbasic variables
            
            % basis inverse
            if iters == counter * reinv
                % recompute the inverse of the basis from scratch every reinv iterations
                BasisInv = inv(Basis);
                counter = counter + 1;
            else
                % Product Form of the Inverse basis update method
                eta = -h_l / h_l(r);
                eta(r) = 1 / h_l(r);
                EInv = speye(length(BasisInv));
                EInv(:, r) = eta;
                BasisInv = EInv * BasisInv;
            end
            
            xb = BasisInv * b;       % basic solution
            w = fb' * BasisInv;      % simplex multiplier
            sn = fn' - w * NonBasis; % reduced costs
            iters = iters + 1;
        end
    end
end

% Phase II
if flag == 1
    disp('---- P H A S E    II ----');
    % matrix A does not contain the artificial variable
    % Here we solve the original LP problem
    [m, n] = size(A); % new size of matrix A
    
    Basis = A(:, BasicList);       % matrix of the basic variable
    cb = c(BasicList);             % coefficients of the objective function for the basic variables
    NonBasis = A(:, NonBasicList); % matrix of the nonbasic variables
    cn = c(NonBasicList);          % coefficients of the objective function for the nonbasic variables
    BasisInv = inv(Basis);         % basis inverse
    xb = BasisInv * b;             % basic solution
    w = cb' * BasisInv;            % simplex multiplier
    sn = cn' - w * NonBasis;       % reduced costs
    
    % iterate in Phase II
    while true
        % optimality test for Phase II
        if all(sn >= 0)
            % the problem is optimal
            % calculate the value of the objective function
            if MinMaxLP == 1
                fval = -(c(BasicList)' * xb + c0);
            else
                fval = c(BasicList)' * xb + c0;
            end
            existflag = 0;
            xsol = xb;
            iters = iters + 1;
            disp('The LP problem is optimal');
            return;
        else
            % the LP problem of Phase II is not optimal
            % find the entering variable using Dantzig's rule
            [~, t] = min(sn);
            l = NonBasicList(t);      % index of the entering variable
            h_l = BasisInv * A(:, l); % pivot column
            
            mrt = find(h_l > 0);
            if ~isempty(mrt)
                % perform the minimum ratio test to find the leaving variable
                xb_div_hl = xb(mrt) ./ h_l(mrt);
                [p, r] = min(xb_div_hl);
                rr = find(xb_div_hl == p);
                % break the ties in order to avoid stalling
                if length(rr) > 1
                    r = mrt(rr);
                    k = BasicList(r); % index of the leaving variable
                    k = max(k);
                    r = find(BasicList == k);
                else
                    r = mrt(r);
                    k = BasicList(r); % index of the leaving variable
                end
                % update basic and nonbasic list
                BasicList(r) = l;
                NonBasicList(t) = k;
                
                Basis = A(:, BasicList);       % matrix of the new basic variable
                cb = c(BasicList);             % coefficients of the objective function for the new basic variables
                NonBasis = A(:, NonBasicList); % matrix of the new nonbasic variables
                cn = c(NonBasicList);          % coefficients of the objective function for the new nonbasic variables
                
                % basis inverse
                if iters == counter * reinv
                    % recompute the inverse of the basis from scratch every reinv iterations
                    BasisInv = inv(Basis);
                    counter = counter + 1;
                else
                    % Product Form of the Inverse basis update method
                    eta = -h_l / h_l(r);
                    eta(r) = 1 / h_l(r);
                    EInv = speye(length(BasisInv));
                    EInv(:, r) = eta;
                    BasisInv = EInv * BasisInv;
                end
                
                xb = BasisInv * b;       % basic solution
                w = cb' * BasisInv;      % simplex multiplier
                sn = cn' - w * NonBasis; % reduced costs
                iters = iters + 1;
            else
                % the LP problem is unbounded
                existflag = 2;
                iters = iters + 1;
                disp('The LP problem is unbounded');
                return;
            end
        end
    end
end
end