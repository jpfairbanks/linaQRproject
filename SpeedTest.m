% Domenic Carr & James Fairbanks
% MATH6643 Project

% Top-level Script for Testing Speed of Updates
clear all

%%%% Parameters to Change
m = 110;
n = 24;
k = 4;

%%%% Start Work
v = RandVec(n,k);
A = rand(m,n);
Acopy =  A;
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
HouseholderQR(Acopy+X)
timeFull = toc();
% test = triu(R)-R;
% norm(test,'fro')
assert(norm(triu(R)-R,'fro') < 1e-10, 'R is not upper triangular')
assert(norm(Q*Q' - eye(m,m)) < 1e-10, 'sparseQR is not orthogonal')
assert(qrerror(A+X, Q, R) < 1e-10, 'sparse update failed')

timeIter = 0;
for i=1:1:k
    X_new = zeros(m,n);
    X_new(:,v(i)) = X(:,v(i));
    tic
    Rcopy = Rcopy + Qcopy'*X_new;
    [Qt,Rcopy] = SparseGivensSingle(Rcopy,v(i));
    Qcopy = Qcopy*Qt';
    time = toc();
    timeIter = timeIter + time;
end
timeSparseEager = timeIter;
assert(norm(triu(Rcopy)-Rcopy,'fro') < 1e-10, 'R is not upper triangular')
assert(norm(Qt'*Qt - eye(m,m)) < 1e-10, 'sparseQR is not orthogonal')
assert(qrerror(Acopy+X, Qcopy, Rcopy) < 1e-10, 'sparse update failed')    
