function [ c,s ] = Givens( a,b )
%Givens returns cos and sin values for Givens rotation matrix
%       based on a and b
%
%       Follows Algorithm 5.1.3 from Golub & Van Loan
%
%       a is a scalar
%       b is a scalar
%
%       c is a scalar (cosine value of rotation)
%       s is a scalar (sine value of rotation)
%

    if b == 0
        c = 1;
        s = 0;
    else
        if abs(b) > abs(a)
            t = -a/b;
            s = 1/sqrt(1+t*t);
            c = s*t;
        else
            t = -b/a;
            c = 1/sqrt(1+t*t);
            s = c*t;
        end
    end

end