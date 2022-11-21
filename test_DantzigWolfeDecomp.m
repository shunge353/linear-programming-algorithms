%% Test for the function DantzigWolfeDecomp
% data from page 292 Example 10.5 & 10.6 of George B.Dantzig and Mukund N.Thapa
% ¡¶Linear Programming2 Theory and Extensions¡·.
% there are three subproblems, sub1 and sub2 have bounded feasible region S1 and S2,
% whereas sub3 has unbouned feasible region S3.
clc; clear; close all;

A1 = [3 2 1 6 5 4;
    1 8 3 7 1 4];
A2 = [8 5 7 3 4 1;
    5 2 5 3 2 6];
A3 = [1 2;
    3 4];
c1 = [1 2 3 4 5 6];
c2 = [1 2 3 4 5 6];
c3 = [7 -10];
B1 = [1 1 1 0 0 0;
    0 0 0 1 1 1;
    1 0 0 1 0 0;
    0 1 0 0 1 0;
    0 0 1 0 0 1];
B2 = [1 1 1 0 0 0;
    0 0 0 1 1 1;
    1 0 0 1 0 0;
    0 1 0 0 1 0;
    0 0 1 0 0 1];
B3 = [1 -1];
b = [64; 63];
b1 = [3; 4; 2; 1; 4];
b2 = [4; 5; 3; 3; 3];
b3 = [1];

extr.A1 = A1; extr.A2 = A2; extr.A3 = A3;
extr.c1 = c1; extr.c2 = c2; extr.c3 = c3;
sub.B1 = B1; sub.B2 = B2; sub.B3 = B3;
extr.b = b; sub.b1 = b1; sub.b2 = b2; sub.b3 = b3;
n = 3;

% Solve the decomposed problem by Dantzig-Wolfe decomposition
[xsol, fval, existflag, iters] = DantzigWolfeDecomp(extr, sub, n);

% Solve the decomposed problem by cplex
n1 = size(c1, 1);
n2 = size(c2, 1);
n3 = size(c3, 1);
x1 = sdpvar(n1, 1);
x2 = sdpvar(n2, 1);
x3 = sdpvar(n3, 1);
obj = c1' * x1 + c2' * x2 + c3' * x3;
const = [];
const = [const, A1 * x1 + A2 * x2 + A3 * x3 == b];
const = [const, B1 * x1 == b1];
const = [const, B2 * x2 == b2];
const = [const, B3 * x3 == b3];
const = [const, x1 >= 0];
const = [const, x2 >= 0];
const = [const, x3 >= 0];
opts = sdpsettings('solver', 'cplex', 'verbose', 2);
diag = optimize(const, obj, opts);