function L = get_laplacian(W)
% get normalized graph laplacian of adjecency matrix W
% L = I - D^(-1)*W

% D=sum(W)+sum(W,2)';
% D=diag(1./D);
% L = eye(length(D)) - D*W;

D=diag(sum(W)+sum(W,2)');
L=D-W;