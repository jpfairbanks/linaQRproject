% Domenic Carr & James Fairbanks
% MATH6643 Project

% Testing Script for QR by Givens

display('new')

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = GivensQR(A)

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = HouseholderQR(A)

A = [2 1 5 4; -1 3 1 2; 2 0 1 1; -1 0 0 -3];
display(A);
ind = ones(1,1);

[Qt,R] = SparseGivens(A,ind);
Qt'
R
A = [2 1 5 4; -1 3 1 2; 2 0 1 1; -1 0 0 -3];
[Q,R] = qr(A)


B = [2 1 5 4 -2; 1 -3 1 2 3; 0 2 0 1 -4; 4 -1 0 0 -3; 3 -2 0 0 1]

[Qt,R] = SparseGivens(B,[1,2]);
Qt'
R

[Q,R] = qr(B)
