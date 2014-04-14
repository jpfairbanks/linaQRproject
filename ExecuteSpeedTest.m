function [ tFull,tEager,tLazy ] = ExecuteSpeedTest( m,n,k )
%ExecuteSpeedTest computes the times to perform QR factorizations
%       using full, sparse-eager, and sparse-lazy methods
%
%
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
    tLazy = time1+time2+time3;

    tic
    [Q_full, R_full] = HouseholderQR(Acopy+X);
    tFull = toc();

    tEager = 0;
    for i=1:1:k
        tic
        Rcopy(:,v(i)) = Rcopy(:,v(i)) + Qcopy'*X(:,v(i));
        [Qt,Rcopy] = SparseGivensSingle(Rcopy,v(i));
        Qcopy = Qcopy*Qt';
        time = toc();
        tEager = tEager + time;
    end

end