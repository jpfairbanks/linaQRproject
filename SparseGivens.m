function [ Qt, R ] = SparseGivens( A,ind )
%SparseGivens returns Q' and R for a given matrix A (nearly upper
%       triangular) such that Q'R = R
%
%       Uses ideas from Algorithm 5.2.4 from Golub & Van Loan
%
%       Q' is an mxm matrix
%       R is an mxn matrix 
%
%       A is an mxn matrix with m >= n (should be nearly triangular)
%       ind is a vector of column indices that need zeroing
%

try
    [m,n] = size(A);
    if (m < n)
        error('dimensions of A are not correct! (require: m >= n)');
    end
    
    Qt = eye(m);
    R = A;
    
    for j = 1:size(ind)
        ndx = ind(j);
        % this loop removes non-zero column entries
        for i = m:-1:ndx+1
            G = eye(m);
            [c,s] = Givens(R(i-1,ndx),R(i,ndx));
            R(i-1:i,ndx:n) = [c -s; s c]*R(i-1:i,ndx:n);
            G([i-1,i],[i-1,i]) = [c s; -s c];
            Qt = G*Qt;            
        end
        
        % this takes care of fill-ins
        if j ~= size(ind)
            % this loop removes fill-in between non-zero columns
            for i = ndx+1:1:ind(j+1)-1
                G = eye(m);
                [c,s] = Givens(R(i,i),R(i+1,i));
                R(i:i+1,i:n) = [c -s; s c]*R(i:i+1,i:n);
                G([i,i+1],[i,i+1]) = [c s; -s c];
                Qt = G*Qt;            
            end
        else
            % this loop removes all remaining fill-ins
            for i = n-size(ind):1:n-1
                for k = i+1:1:n
                    G = eye(m);
                    [c,s] = Givens(R(i,i),R(k,i));
                    R([i,k],i:n) = [c -s; s c]*R([i,k],i:n);
                    G([i,k],[i,k]) = [c s; -s c];
                    Qt = G*Qt;
                end
            end
        end
    end

catch err
    throw(err);
end
    
end