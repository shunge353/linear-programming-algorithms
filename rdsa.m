function [xsol, fval, existflag, iters] = rdsa(A, c, b, Eqin, MinMaxLP, c0, reinv)
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

% set default value to missing inputs
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

disp('---- I N I T I A L    B A S I S ----');
% find an invertiable basis using rref
[R, p] = rref([A b]);
BasicList = sort(p);
NonBasicList = setdiff([1:n], BasicList);

% delete redundant constraints if exists
ind = all(R, 2);
A(ind, :) = [];
b(ind) = [];

[m, n] = size(A); % new size of the matrix A

Basis = A(:, BasicList); % matrix of the basic variable
cb = c(BasicList); % coefficients of the objective function for the basic variable
NonBasis = A(:, NonBasicList); % matrix of the nonbasic variable
cn = c(NonBasicList); % coefficients of the objective function for the nonbasic variable
BasisInv = inv(Basis); % basis inverse
xb = BasisInv * b; % basic solution
w = cb' * BasisInv; % simplex multiplier
sn = cn' - w * NonBasis; % reduced costs

% check if the current basis is dual feasible
if all(sn >= 0)
    flag = 1;
else
    flag = 0;
end

% modified big-M method
if flag == 0
    disp('---- D U A L    W I T H    B I G - M    M E T H O D ----');
    % find the entering variable
    [p, t] = min(sn);
    rr = find(sn == p);
    % break the ties in order to aviod stalling
    if length(rr) > 1
        l = NonBasicList(rr);
        l = max(l);
    else
        l = NonBasicList(t);
    end
    % add a constraint and an artificial variable
    A(m + 1, :) = zeros(1, n);
    A(m + 1, NonBasicList) = 1;
    A = [A [zeros(m, 1); 1]];
    c = [c; 0];
    b = [b; 0];
    bM = [zeros(m, 1); 1]; % right-hand side of the big-M method
    % compute the new basic and nonbasic variables
    BasicList(m + 1) = l;
    NonBasicList(NonBasicList == l) = n + 1;
    
    artificialVariableInN = 1; % artificial variable is in the nonbaic list
    
    Basis = A(:, BasicList); % matrix of the new basic variable
    cb = c(BasicList); % coefficients of the objective function for the new basic variable
    NonBasis = A(:, NonBasicList); % matrix of the new nonbasic variable
    cn = c(NonBasicList); % coefficients of the objective function for the new nonbasic variable
    BasisInv = inv(Basis); % new basis inverse
    xb = BasisInv * b; % basic solution for the original problem
    xbM = BasisInv * bM; % basic solution for the big-M problem
    w = cb' * BasisInv; % new simplex multiplier
    sn = cn' - w * NonBasis; % new reduced costs
    
    % iterate in big-M method
    while true
        % optimality test for big_M method
        if all(xbM > 0)
            col = find(BasicList == (n + 1));
            if ~isempty(col)
                % if the artificial variable is in the basic list, the problem is optimal
                if MinMaxLP == 1
                    fval = -(cb' * xb + c0);
                else
                    fval = cb' * xb + c0;
                end
                existflag = 0;
                xsol = xb;
                iters = iters + 1;
                disp('The LP problem is optimal');
                return;
            else
                % if the artificial variable is not in the basic list, check if the reduced cost of the artificial variable
                % is equal to zero
                w = cb' * BasisInv; % simplex multiplier
                s = c' - w * A; % reduced costs
                if s(n + 1) == 0
                    % if the reduced cost of the artificial variable is equal to zero, the problem is optimal
                    if MinMaxLP == 1
                        fval = -(cb' * xb + c0);
                    else
                        fval = cb' * xb + x0;
                    end
                    existflag = 0;
                    xsol = xb;
                    iters = iters + 1;
                    disp('The LP problem is optimal');
                    return;
                else
                    % the problem is unbounded
                    existflag = 2;
                    iters = iters + 1;
                    disp('The LP rpoblem is unbouned');
                    return;
                end
            end
        elseif all(xbM >= 0)
            row = find(xbM == 0);
            if ~isempty(row)
                if all(xb(row) >= 0)
                    col = find(BasicList == (n + 1));
                    if ( ~isempty(col) && xbM(col) == 0 && xb(col) == 0 ) || isempty(col)
                        % if the artificial variable is in the basic list, check if the reduced cost of the artificial
                        % variable is equal to zero
                        w = cb' * BasisInv; % simplex multiplier
                        s = c' - w * A; % reduced costs
                        if s(n + 1) == 0
                            % if the reduced costs of the artificial variable is equal to zero, the problem is optimal
                            if MinMaxLP == 1
                                fval = -(cb' * xb + c0);
                            else
                                fval = cb' * xb + c0;
                            end
                            existflag = 0;
                            xsol = xb;
                            iters = iters + 1;
                            disp('The LP problem is optimal');
                            return;
                        else
                            % the problem is unbounded
                            existflag = 2;
                            iters = iters + 1;
                            diap('The LP problem is unbounded');
                            return;
                        end
                    elseif ( ~isempty(col) && xbM(col) > 0 ) || xb(col) > 0
                        % the problem is optimal
                        if MinMaxLP == 1
                            fval = -(cb' * xb + c0);
                        else
                            fval = cb' * xb + c0;
                        end
                        existflag = 0;
                        xsol = xb;
                        iters = iters + 1;
                        disp('The LP rpoblem is optimal');
                        return;
                    end
                else
                    % the artificial variable left the nonbasic list
                    artificialVariableInN = 0;
                end
            end
        end
        
        % if the aritificial variable is in the nonbasic list, use xbM to select the leaving variable
        if artificialVariableInN == 1
            % find the leaving variable
            mrt = find(xbM < 0);
            [a, r] = min(xbM(mrt));
            rr = find(xbM(mrt) == a);
            % break the ties in order to avoid stalling
            if length(rr) > 1
                r = max(rr);
                k = BasicList(r);
                k = min(k);
                r = find(BasicList == k);
            else
                r = mrt(rr);
                k = BasicList(r);
            end
        else
            % if the artificial variable is not in the nonbasic list, use xb to select the leaving variable
            mrt = xb(row) < 0;
            mrt = row(mrt);
            [a, r] = min(xb(mrt));
            rr = find(xb(mrt) == a);
            % break the ties in order to stalling
            if length(rr) > 1
                r = mrt(rr);
                k = BasicList(r);
                k = min(k);
                r = find(BasicList == k);
            else
                r = mrt(rr);
                k = BasicList(r);
            end
        end
        
        HrN = BasisInv(r, :) * NonBasis;
        mrt = find(HrN < 0);
        % if there is not any candidate to enter the basic list, the problem is infeasible
        if isempty(mrt)
            existflag = 1;
            iters = iters + 1;
            disp('The LP problem is infeasible');
            return;
        end
        
        % perform the minimum ratio test to select the entering variable
        [a, t] = min(-sn(mrt) ./ HrN(mrt));
        rr = find(-sn(mrt) ./ HrN(mrt) == a);
        % break the ties in order to avoid stalling
        if length(rr) > 1
            l = NonBasicList(mrt(rr));
            l = max(l);
        else
            l = NonBasicList(mrt(t));
        end
        
        h_l = BasisInv * A(:, l); % pivot column
        % check if the problem is unbounded
        if all(h_l <= 0)
            existflag = 2;
            iters = iters + 1;
            disp('The LP problem is unbounded');
            return;
        end
        
        % update basic and nonbasic list
        BasicList(r) = l;
        NonBasicList(t) = k;
        
        Basis = A(:, BasicList); % matrix of the basic variable
        cb = c(BasicList); % coefficients of the objective function for the basic variable
        NonBasis = A(:, NonBasicList); % matrix of the nonbasic variable
        cn = c(NonBasicList); % coefficients of the objective function for the nonbasic variable
        
        % basis inverse
        if iters == reinv * counter
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
        
        xb = BasisInv * b; % basic solution for the original problem
        xbM = BasisInv * bM; % basic solution for the big-M problem
        w = cb' * BasisInv; % simplex multiplier
        sn = cn' - w * NonBasis; % reduced costs
        iters = iters + 1;
    end
end

% dual simplex method
if flag == 1
    disp('---- D U A L    S I M P L E X    M E T H O D----');
    while true
    % optimality test
    if all(xb >= 0)
        % the problem is optimal
        if MinMaxLP == 1
            fval = -(cb' * xb + c0);
        else
            fval = cb' * xb + c0;
        end
        existflag = 0;
        xsol = xb;
        iters = iters + 1;
        disp('The LP problem is optimal');
        return;
    end
    
    % find the leaving variable
    mrt = find(xb < 0);
    [a, r] = min(xb(mrt));
    rr = find(xb(mrt) == a);
    % break the ties in order to avoid stalling
    if length(rr) > 1
        r = mrt(rr);
        k = BasicList(r);
        k = min(k);
        r = find(BasicList == k);
    else
        r = mrt(rr);
        k = BasicList(r);
    end
    
    HrN = BasisInv(r, :) * NonBasis;
    mrt = find(HrN < 0);
    % if there is not any candidate to enter the basic list, then the problems is infeasible
    if isempty(mrt)
        existflag = 1;
        iters = iters + 1;
        disp('The LP problem is infeasible');
        return;
    end
    
    % perform the minimum ratio test to select the entering variable
    [a, t] = min(-sn(mrt) ./ HrN(mrt));
    rr = find(-sn(mrt) ./ HrN(mrt) == a);
    % break the ties in order to avoid stalling
    if length(rr) > 1
        l = NonBasicList(mrt(rr));
        l = max(l);
    else
        l = NonBasicList(mrt(t));
    end
    
    h_l = BasisInv * A(:, l); % pivot column
    % check if the problem is unbounded
    if all(h_l <= 0)
        existflag = 2;
        iters = iters + 1;
        disp('The LP problem is unbounded');
        return;
    end
    
    % update the basic and nonbasic list
    BasicList(r) = l;
    NonBasicList(t) = k;
    
    Basis = A(:, BasicList); % matrix of the basic variable
    cb = c(BasicList); % coefficients of the objective function for the basic variable
    NonBasis = A(:, NonBasicList); % matrix of the nonbasic variable
    cn = c(NonBasicList); % coefficients of the objective function for the nonbasic variable
    
    % basis inverse
    if iters == reinv * counter
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
      
    xb = BasisInv * b; % basic solution for the original problem
    w = cb' * BasisInv; % simplex multiplier
    sn = cn' - w * NonBasis; % reduced costs
    iters = iters + 1;
    end
end
end
