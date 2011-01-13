function X = scheinerman_tucker_latent_features(A,d,epsilon)
n = length(A);
D = zeros(n);
Dold = ones(n);
X = ones(d,n);
Xold = zeros(d,n);

if nargin < 3
    epsilon = 10^-8;
end
k=1;
while norm(D-Dold,'fro')>epsilon
    Xold = X;
    Dold = D;
    [U,Lambda] = eigs(A+D,d,'la');
    Lambda = max(zeros(d),Lambda);
    X =  Lambda.^(.5)*U';
    D = diag(diag(X'*X));
    k=k+1;
end
k
end
