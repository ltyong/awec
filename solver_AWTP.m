function [C,Z] = solver_AWTP(CA,A,lambda,gamma,K)
n = size(CA,1);
t = 0;
e = 1e-2;
max_iter = 100;
I = eye(n);

E = zeros(n,n);
C = I;
F = C;
mu = 0.1;
Y1 = ones(n);
Y2 = ones(n);
rho = 1.2;
max_mu = 1e8;
[~,order] = size(A);

alpha = ones(1,order)./order;

Z = C;
L = diag(sum(Z)) - Z;
[H0,~,~] = eig1(L,K,0);
H = H0(:,1:K);
gam = 1.1;
D = diag(sum(C));
pre_com = inv(2 * order*I + (mu + mu) * I);
while t < max_iter
    t = t + 1;
    
    % update C

    P1 = CA - E + Y1 / mu;
    P2 = F - Y2 / mu;
    alphaA = zeros(n,n);
    for i = 1:order
    alphaA = alpha(i)*A{i}+alphaA;
    end
    G = diag(1./sqrt(diag(D)))*Z;
    C = pre_com * (mu * P1 + mu * P2 + G*G' + 2*alphaA);
    D = diag(sum(C));

    % update E
    E = mu * (CA - C) + Y1;

    E = E /(lambda + mu);
    
    % update F
    F = C + Y2 / mu;
    F = min(max((F + F') / 2,0),1);
    
    % update Z
    Z = update_Z(C,Z,gamma,gam,H);

    sumTrZHZZ = 0;
    suma = 0;
    for i = 1:order
     
        TrAtA(i) = trace(A{i}'*A{i});
        TrCtA(i) = trace(C'*A{i});
        TrCAAA= TrCtA(i)/TrAtA(i);
        sumTrZHZZ = sumTrZHZZ +TrCAAA;
        
        suma= suma+1/TrAtA(i);
    end

    beta = (sumTrZHZZ-1)*2/suma;
    
    for i =1:order
        alpha(i) = (2*TrCtA(i)-beta)/(2*TrAtA(i));
    end
    
     L_z = diag(sum(Z))-Z;
     %phi =L_z;
     [H, ~, ev]=eig1(L_z, K, 0);
    
     % update Y
     
    Y1t = Y1;
    residual1 = CA - C - E;
    
    Y2t = Y2;
    residual2 = C - F;

     stopC = max(max(abs(residual1),abs(residual2)));
    if stopC < e
        break;
      else
        Y1 = Y1 + mu*residual1;
        Y2 = Y2 + mu*residual2;
        mu = min(max_mu,mu*rho);
    end
 
end

end