%% Test for the function ftdsa
% optimal
clc; clear; close all;

A = [2 -1 -2 -3;
     0 -2 -2 -6];
c = [1; 0; 1; 3];
b = [-4; -6];

[xsol, fval, existflag, iters, ray] = ftdsa(A, c, b);

%% Test for the function ftdsa
% Exercise 3.17 (optimal)
clc; clear; close all;

A = [1  3 0  4 1;
     1  2 0 -3 1;
    -1 -4 3  0 0];
c = [2; 3; 3; 1; -2];
b = [2; 2; 1];

[xsol, fval, existflag, iters, ray] = ftdsa(A, c, b);

%% Test for ftpsa
% Infeasible
clc; clear; close all;

A = [1  1  1;
     0 -2 -3];
c = [0; -1; -2];
b = [5; 3];

[xsol, fval, existflag, iters, ray] = ftdsa(A, c, b);

%% Test for ftpsa
% unbounded
clc; clear; close all;

A = [2  4 -1  0 0;
     1  1  0 -1 0;
    -2 -7  0  0 1];
c = [3; -2; 0; 0; 0];
b = [8; 3; -15];

[xsol, fval, existflag, iters, ray] = ftdsa(A, c, b);