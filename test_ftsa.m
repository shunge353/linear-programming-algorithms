%% Test for the function ftsa
% Example 3.8
clc; clear; close all;

A = [ 1 2 3 0;
     -1 2 6 0;
      0 4 9 0;
      0 0 3 1];
c = [1; 1; 1; 0];
b = [3; 2; 5; 1];

[xsol, fval, existflag, iters] = ftsa(A, c, b);

%% Test for the function ftsa
% Example 3.5
clc; clear; close all;

A = [1 2 2 1 0 0;
     2 1 2 0 1 0;
     2 2 1 0 0 1];
c = [-10; -12; -12; 0; 0; 0];
b = [20; 20; 20];

[xsol, fval, existflag, iters] = ftsa(A, c, b);

%% Test for ftsa
% Infeasible
clc; clear; close all;

A = [1  1  1;
     0 -2 -3];
c = [0; -1; -2];
b = [5; 3];

[xsol, fval, existflag, iters] = ftsa(A, c, b);

%% Test for ftsa
% unbounded
clc; clear; close all;

A = [2  4 -1  0 0;
     1  1  0 -1 0;
    -2 -7  0  0 1];
c = [3; -2; 0; 0; 0];
b = [8; 3; -15];

[xsol, fval, existflag, iters] = ftsa(A, c, b);