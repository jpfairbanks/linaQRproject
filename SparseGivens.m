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
%       ind is a non-decreasing vector of column indices that need zeroing
%

try
    [m,n] = size(A);
    if (m < n)
        error('dimensions of A are not correct! (require: m >= n)');
    end
    
    Qt = eye(m);
    R = A;
    
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
%         display('oh no we created fill in');
%         display(R);
        % this takes care of fill-ins
        if j ~= length(ind)
            % this loop removes fill-in between non-zero columns
            for i = ndx+1:1:ind(j+1)-1
                if m > n
                    val = min(j,m-i);
                else
                    val = min(j,n-i);
                end
                for k = i+1:1:i+val
                    [c,s] = Givens(R(i,i),R(k,i));
                    R([i,k],i:n) = [c -s; s c]*R([i,k],i:n);
                    Qt([i,k],:) = [c -s; s c]*Qt([i,k],:);
                end            
            end
%             display('did we catch all of the fill in')
%             display(R);
        else
%             display('fillin else case')
%             display(R);
            % this loop removes all remaining fill-ins
            if m == n
                for i = ndx+1:1:n-1
                    val = min(j,n-i);
                    for k = i+1:1:i+val
                        [c,s] = Givens(R(i,i),R(k,i));
                        R([i,k],i:n) = [c -s; s c]*R([i,k],i:n);
                        Qt([i,k],:) = [c -s; s c]*Qt([i,k],:);
                    end
                end
            else
                for i = ndx+1:1:n
                    val = min(j,m-i);
                    for k = i+1:1:i+val
                        [c,s] = Givens(R(i,i),R(k,i));
                        R([i,k],i:n) = [c -s; s c]*R([i,k],i:n);
                        Qt([i,k],:) = [c -s; s c]*Qt([i,k],:);
                    end
                end
            end
        end
    end

catch err
    throw(err);
end
    
end