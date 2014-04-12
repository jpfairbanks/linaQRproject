function [ Q, R ] = GivensQR( A )
%QRbyGivens computes the QR factorization of A
%
%       Uses ideas from Algorithm 5.2.4 from Golub & Van Loan
%
%       Q is an mxm matrix
%       R is an mxn matrix 
%
%       A is an mxn matrix with m >= n
%

try
    [m,n] = size(A);
    if (m < n)
        error('dimensions of A are not correct! (require: m >= n)');
    end
    
    Q = eye(m);
    R = A;
    
    for j = 1:n
        for i = m:-1:j+1
            [c,s] = Givens(R(i-1,j),R(i,j));
            R(i-1:i,j:n) = [c -s; s c]*R(i-1:i,j:n);
            Q(:,[i-1,i]) = Q(:,[i-1,i])*[c s; -s c];            
        end
    end

catch err
    throw(err);
end
    
end