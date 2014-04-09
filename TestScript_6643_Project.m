% Domenic Carr & James Fairbanks
% MATH6643 Project

% Testing Script for QR by Givens

display('new')

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = GivensQR(A)

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = HouseholderQR(A)

A = [2 1 5; -1 3 0; 2 0 1];
ind = ones(1,1);

[Qt,R] = SparseGivens(A,ind);
Qt'
R
A = [2 1 5; -1 3 0; 2 0 1];
[Q,R] = qr(A)

