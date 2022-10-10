%% Test for MPC
clc; clear; close all;

A = [2  0  2  3;
     0 -2 -2 -6];
c = [1; 0; -1; -3];
b = [10; -6];
Eqin = [0; 0];

[xsol, fval, exisflag, iters] = ipdipm(A, c, b, Eqin);

%% Test for MPC
clc; clear; close all;

A = [6 1 -2;
     1 1  1;
     6 4 -2];
c = [5; 2; -4];
b = [5; 4; 10];
Eqin = [1; -1; 1];

[xsol, fval, exisflag, iters] = ipdipm(A, c, b, Eqin);