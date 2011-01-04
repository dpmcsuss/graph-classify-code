function L = get_laplacian(A)
% get normalized graph laplacian of adjecency matrix A
% L = I - D^(-1)*A

% D=sum(A)+sum(A,2)';
% D=diag(1./D);
% L = eye(length(D)) - D*A;

D=diag(sum(A)+sum(A,2)');
L=D-A;