function [ Q,R ] = HouseholderQR( A )
%HouseholderQR computes the QR factorization of A
%
%       Follows Algorithm 5.2.1 from Golub & Van Loan
%
%       A is an mxn matrix with m > n
%
%       Q is an mxm matrix
%       R is an mxn matrix
%

% be sure both inputs are okay
try
    [m,n] = size(A);
    if (m < n)
        error('dimensions of A are not correct! (require: m >= n)');
    end

    % passed all necessary conditions... now start computations
    
    % compute QR factorization and overwrite into A
    for j = 1:n
        [v,beta] = house(A(j:m,j));
        A(j:m,j:n) = (eye(m-j+1)-beta*(v*v'))*A(j:m,j:n);
        if j < m
            A(j+1:m,j) = v(2:m-j+1);
        end
    end    
    R = triu(A);    
    Q = eye(m);
    for j = n:-1:1
        v = vertcat(1,A(j+1:m,j));
        beta = 2/(1+norm(A(j+1:m,j),2)^2);
        Q(j:m,j:m) = Q(j:m,j:m)-beta*v*(v'*Q(j:m,j:m));
    end
    
catch err    
    throw(err);    
end

end