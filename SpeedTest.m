% Domenic Carr & James Fairbanks
% MATH6643 Project

% Top-level Script for Testing Speed of Updates
clear all

%%%% Parameters to Change
m = 100;
n = 40;
k = 1;

%%%% Start Work
v = RandVec(n,k);
A = rand(m,n);
Acopy =  A;
Acopy2 = A;
X = zeros(m,n);
X(:,v) = rand(m,k);
[Q,R] = HouseholderQR(A);
Qcopy = Q;
Rcopy = R;
tic;
R = R + Q'*X;
time1 = toc();
tic
[Qt,R] = SparseGivens(R,v');
time2 = toc();
tic
Q = Q*Qt';
time3 = toc();
timeSparseLazy = time1+time2+time3;
tic
[Q_full, R_full] = HouseholderQR(Acopy+X);
timeFull = toc();
% tic
% [Q_full2, R_full2] = qr(Acopy2+X);
% timeMatlab = toc();
% test = triu(R)-R;
% norm(test,'fro')
assert(norm(triu(R)-R,'fro') < 1e-10, 'R is not upper triangular')
assert(norm(Q*Q' - eye(m,m)) < 1e-10, 'sparseQR is not orthogonal')
assert(qrerror(A+X, Q, R) < 1e-10, 'sparse update failed')

timeSparseEager = 0;
for i=1:1:k
    tic
    Rcopy(:,v(i)) = Rcopy(:,v(i)) + Qcopy'*X(:,v(i));
    [Qt,Rcopy] = SparseGivensSingle(Rcopy,v(i));
    Qcopy = Qcopy*Qt';
    time = toc();
    timeSparseEager = timeSparseEager + time;
end
assert(norm(triu(Rcopy)-Rcopy,'fro') < 1e-10, 'R is not upper triangular')
assert(norm(Qt'*Qt - eye(m,m)) < 1e-10, 'sparseQR is not orthogonal')
assert(qrerror(Acopy+X, Qcopy, Rcopy) < 1e-10, 'sparse update failed')    
