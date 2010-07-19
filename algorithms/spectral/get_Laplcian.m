function L = get_Laplcian(W)
% get normalized graph laplacian of adjecency matrix W
% L = I - D^(-1)*W

D=sum(W)+sum(W,2)';
D=diag(1./D);
L = eye(length(D)) - D*W;