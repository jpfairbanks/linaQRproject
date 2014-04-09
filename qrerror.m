function [r] = qrerror(A,Q,R)
  r = norm(A-(Q*R), 'fro');
end