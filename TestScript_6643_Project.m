% Domenic Carr & James Fairbanks
% MATH6643 Project

% Testing Script for QR by Givens

display('new')

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = GivensQR(A);
assert(qrerror(A,Q,R) < 1e-10, 'GivensQR failed');

A = [2 1 5; -1 3 -1; 2 1 1];
[Q,R] = HouseholderQR(A);
assert(qrerror(A,Q,R) < 1e-10, 'HouseholderQR failed');

A = [2 1 5 4; -1 3 1 2; 2 0 1 1; -1 0 0 -3];
display(A);
ind = ones(1,1);

[Qt,R] = SparseGivens(A,ind);
Qt'*R
A = [2 1 5 4; -1 3 1 2; 2 0 1 1; -1 0 0 -3]
assert(qrerror(A,Qt',R) < 1e-10, 'SparseGivens failed');
[Q,R] = qr(A);

B = [2 1 5 4 -2; 1 -3 1 2 3; 0 2 0 1 -4; 4 -1 0 0 -3; 3 -2 0 0 1]

[Qt,R] = SparseGivens(B,[1,2]);
display('Qttranspose');
Qt'
R

[Q,R] = qr(B);

% SparseQR fails for nonsquare A
m = 8;
n = 4;
A = rand(m,n)
[Q, R] = qr(A)
deltaA = zeros(m,n);
deltaA(m,n/2) = 1
B = R + Q'*deltaA
[Qt, R] = SparseGivens(B, [n/2])
display('R is not upper triangular')
QtR = Qt*R
assert(norm(triu(R)-R,'fro') < 1e-10, 'R is not upper triangular')

% cannot recover B from Qt*R
m = 8;
n = 8;
A = rand(m,n)
[Q, R] = qr(A)
deltaA = zeros(m,n);
deltaA(m,n/2) = 1
[Q1,R1] = qr(A+deltaA)
B = R + Q'*deltaA
[Qt, R] = SparseGivens(B, [n/2])
%test that Qt and R satisfy orthogonality and upper triangularity
assert(norm(triu(R)-R,'fro') < 1e-10, 'R is not upper triangular')
assert(norm(Qt'*Qt - eye(m,m)) < 1e-10, 'sparseQR is not orthogonal')
B
% none of these equal B -- now they do
QR  = Qt'*R
% none of these equal R -- now they do
QtB = Qt*B
assert(qrerror(B, Qt', R) < 1e-10, 'sparse update failed')