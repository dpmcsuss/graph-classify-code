function [ U, S, V ] = graph_svd( L )
%graph_svd Computes the SVD of the Laplacian of an adjacancy matrix or
%matrices
%   Returns the standard U,S,V matrix or matrices for each Laplacian matrix
%   in L. Ie L can be either a 2d square matrix or a 3d array of square
%   matrices.
%   Q: Do we need to keep U or V or should we just return a vector of
%   singular values?
d =  size(L);
if numel(d) == 2
    [U,S,V] = svd(L);
    return;
end
U = zeros(d);
S = U;
V = U;
for k=1:d(3)
    [U(:,:,k),S(:,:,k),V(:,:,k)] = svd(L(:,:,k)); 
end

end

