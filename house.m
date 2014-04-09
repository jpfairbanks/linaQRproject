function [ v,beta ] = house( x )
%house returns Householder vector v and beta such that
%       P = I - beta*v'*v
%
%       Follows Algorithm 5.1.1 from Golub & Van Loan
%
%       x is an mx1 vector
%
%       v is an mx1 vector
%       beta is a real number
%

    m = length(x);
    if m > 1
        sigma = x(2:m)'*x(2:m);
        v = vertcat(1,x(2:m));
    
        if ((sigma == 0) && (x(1) >= 0))
            beta = 0;
        elseif ((sigma == 0) && (x(1) < 0))
            beta = -2;
        else
            mu = sqrt(x(1)*x(1)+sigma);        
            if (x(1) <= 0)
                v(1) = x(1) - mu;
            else
                v(1) = -sigma/(x(1)+mu);
            end

            beta = 2*v(1)*v(1)/(sigma+v(1)*v(1));
            v = v/v(1);
        end
    else
        v = x;
        beta = 0;
    end

end