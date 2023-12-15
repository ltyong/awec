function [Z] = update_Z(C,Z,alpha,gamma,F)

[~,n] = size(C);
I = eye(n);
max_iter = 10;
%zr = 10e-11;

L_c = diag(sum(C))-C;
A = L_c+alpha*I;
diff = gamma*L2_distance_1(F',F');

for i = 1:n
%index = find(C(:,i)>0);
% Ii = I(index,i);
% si = S(index,i);
b = 2*alpha*I(:,i)- diff(:,i);
rho = 1.5;
yita = 10;
q = ones(n,1);
t = 0;
Z(:,i) = ones(n,1)/n;
%obj = Z(:,i)'*A*Z(:,i)-Z(:,i)'*b;
while t <max_iter
t = t + 1;
p = Z(:,i)-1/yita*(A'*Z(:,i)+q);
temp = p -1/yita*q-(A*p-b)/yita; 
Z(:,i) = EProjSimplex_new(temp);
yita = yita*rho;
q = q+yita*(Z(:,i)-p);
%obj = [obj;Z(:,i)'*A*Z(:,i)-Z(:,i)'*b];

end
end
%Z = Z+I;
Z=(Z+Z')/2;

end
