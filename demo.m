
clear;
addpath('./ClusteringMeasure');
n = 20;                   
data = 'appendicitis_base_clustering.mat';                   
load(data);       
gt = y;
K = length(unique(gt));
lambda_all = [0.01,0.02,0.04,0.08,0.1,0.2];
gamma_all = [0.1,0.5,1,5,10,50];
for lambda_ind = 1:length(lambda_all)
    lambda = lambda_all(lambda_ind);
for gamma_ind =1:length(gamma_all)
    gamma = gamma_all(gamma_ind);

parfor i = 1:10

[a,b] = size(E);
zz = RandStream('mt19937ar','Seed',i);
RandStream.setGlobalStream(zz);
indx = randperm(b);
EC_end = E(:,indx(1:n));
M = Gbe(EC_end); 
CA = M*M'./n;
NNRate = 0.5;
[KHN] = V9_LocalKernelCalculation(CA , NNRate, K);    
A = V9_generate_multi_connec(KHN,2);

[C,Z] = solver_AWTP(CA,A,lambda,gamma,K);
 [U,S,V] = svd(Z,'econ');
S = diag(S);
r = sum(S>1e-4*S(1));
U = U(:,1:r);S = S(1:r);
U = U*diag(sqrt(S));
U = normr(U);
L = (U*U').^4;
results_S  = spectral_clustering(L,K);


s = squareform(Z - diag(diag(Z)),'tovector');
d = 1 - s;
results_H = cluster(linkage(d,'average'),'maxclust',K);
if min(results_H) == 0
    results_H = results_H + 1;
end

NMI_S(i)= compute_nmi(results_S,gt);
ARI_S(i)= RandIndex(results_S,gt);
 Fscore_S(i)= compute_f(results_S,gt);
 
NMI_H(i)= compute_nmi(results_H,gt);
ARI_H(i)= RandIndex(results_H,gt);
Fscore_H(i)= compute_f(results_H,gt);

end

nmi_S_final(lambda_ind,gamma_ind) = mean(NMI_S);
ari_S_final(lambda_ind,gamma_ind) = mean(ARI_S);
F_S_final(lambda_ind,gamma_ind) = mean(Fscore_S);

nmi_z_final(lambda_ind,gamma_ind) = mean(NMI_H);
ari_z_final(lambda_ind,gamma_ind) = mean(ARI_H);
F_z_final(lambda_ind,gamma_ind) = mean(Fscore_H);
end
end

