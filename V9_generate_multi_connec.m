function [KHL] = V9_generate_multi_connec(KHN,n)

KHL{1} = KHN;
KHL{1} = normalized(KHL{1});
for i = 2:n
KHL{i} = KHL{i-1}'* KHL{1};

KHL_P = KHL{i}(KHL{i} > 0);
mean_value = mean(KHL_P);
std_value = std(KHL_P);
KHL{i}(KHL{i} < mean_value - std_value/2) = 0;
KHL{i} = normalized(KHL{i});
end

end

function Y = normalized(K)
eye_matrix = 1 - eye(size(K));
K = K .* eye_matrix;
c_diag = sum(K, 2);

c_diag(c_diag < 10^(-10)) = 10^(-10);
D = diag(sqrt(c_diag.^(-1)));
Y = D *K *D;
end