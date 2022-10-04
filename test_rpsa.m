%% Test for rpsa
clc; clear; close all;

A = [ 6  1 -2;
      1  1  1;
      6  4 -2];
c = [5; 2; -4];
b = [5; 4; 10];
Eqin = [1; -1; 1];

[xsol, fval, existflag, iters] = rpsa(A, c, b, Eqin);

%% Test for rpsa
clc; clear; close all;

A = [2  0  2  3;
     0 -2 -2 -6];
c = [1; 0; -1; -3];
b = [10; -6];
Eqin = [0; 0];

[xsol, fval, existflag, iters] = rpsa(A, c, b, Eqin);

%% Test for rpsa
clc; clear; close all;

A = [1  1  1;
     0 -2 -3];
c = [0; -1; -2];
b = [5; 3];
Eqin = [0; 0];

[xsol, fval, existflag, iters] = rpsa(A, c, b, Eqin);

%% Test for rpsa
clc; clear; close all;

A = [2  4;
     1  1;
    -2 -7];
c = [3; -2];
b = [8; 3; -15];
Eqin = [1; 1; -1];

[xsol, fval, existflag, iters] = rpsa(A, c, b, Eqin); 