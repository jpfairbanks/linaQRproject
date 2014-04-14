function [ Qt, R ] = SparseGivensSingle( A,ind )
%SparseGivensSingle returns Q' and R for a given matrix A (nearly upper
%       triangular) such that Q'R = R
%
%       Uses ideas from Algorithm 5.2.4 from Golub & Van Loan
%
%       Q' is an mxm matrix
%       R is an mxn matrix 
%
%       A is an mxn matrix with m >= n (should be nearly triangular)
%       ind is a non-decreasing vector of column indices that need zeroing
%
%       This function differs from SparseGivens in that it only allows for
%       processing one column and one fill-in diagonal at a time
%

try
    [m,n] = size(A);
    if (m < n)
        error('dimensions of A are not correct! (require: m >= n)');
    end
    
    Qt = eye(m);
    R = A;   
    
    if (length(ind) > 1)
        error('too many index values(require: length(ind) == 1)');
    end
    
    if (ind(1) == n)
        if  m == n
            return;
        end
    end
    
    for j = 1:1:length(ind)
        ndx = ind(j);
        % this loop removes non-zero column entries
        for i = m:-1:ndx+1
            [c,s] = Givens(R(i-1,ndx),R(i,ndx));
            R(i-1:i,ndx:n) = [c -s; s c]*R(i-1:i,ndx:n);
            Qt([i-1,i],:) = [c -s; s c]*Qt([i-1,i],:);            
        end

        % this takes care of fill-ins
        for i = ndx+1:1:n
            if i+1 > m
                break;
            else
                [c,s] = Givens(R(i,i),R(i+1,i));
                R(i:i+1,i:n) = [c -s; s c]*R(i:i+1,i:n);
                Qt(i:i+1,:) = [c -s; s c]*Qt(i:i+1,:);
            end
        end
    end

catch err
    throw(err);
end
    
end