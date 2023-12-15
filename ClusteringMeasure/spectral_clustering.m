function [ids] = spectral_clustering(W, k)

D = diag(1./sqrt(sum(W, 2)));
W = D * W * D;
[U, s, V] = svd(W);
V = U(:, 1 : k);
V = normr(V);
ids = kmeans(V, k);
%ids = kmeans(V, k, 'emptyaction', 'singleton', 'replicates', 100, 'display', 'off');

end
