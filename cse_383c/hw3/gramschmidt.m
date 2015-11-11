function [Q, R] = gramschmidt(A, flag)

if (nargin == 1) flag = false; end;
    
[m, n] = size(A);
R = zeros(m, n);
Q = zeros(m,n);

if (flag)

    for j=1:n
      q = A(:,j);
      for i=1:j-1
          R(i,j) = Q(:,i)'*A(:,j);
          q = q - R(i,j)*Q(:,i);
          
      end
      R(j, j) = norm(q);
      Q(:,j) = q/R(j,j);
      
    end
    
    
    
else
    
   for j=1:n
      q = A(:,j);
      for i=1:j-1
          R(i,j) = Q(:,i)'*q;
          q = q - R(i,j)*Q(:,i);
          
      end
      R(j, j) = norm(q);
      Q(:,j) = q/R(j,j);
      
    end
    
end