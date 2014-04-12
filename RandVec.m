function [ x ] = RandVec( n,k )
%RandVec returns a random vector with k non-repeating
%       random values chosen in the span 1->n       
%
%       n is the upper bound
%       k is the number of elements to be chosen between 1 and n
%
%       x is a k-vector with random non-repeating values between a and n
%
%
    
    % choose k non-repeating values in the span 1 -> n
    x = zeros(1,k);
    p = zeros(1,n);
    p(1:n) = 1:n;
    for i = n:-1:n-k+1
        q = ceil(i*rand);
        x(1,n-i+1) = p(q);
        p(q:i-1) = p(q+1:i);
    end
    
    x = sort(x);

end