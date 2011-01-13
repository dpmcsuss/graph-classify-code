function [ latIn, latOut ] = estimate_latent_features_SVD( A, dIn, dOut )
%estimate_latent_features_SVD Returns estimated in and out latent features
%   A is an (or many) nxn (or nxnxk) adjacency matrix). dIn and dOut
%   optionally indicate the dimensions of the in and out feature vectors.
%   The default is n for each. The feature vectors are estimated using the
%   SVD method. latIn is dIn x n x k and latOut is dOut x n x k.

sz =  size(A);

if nargin == 1
    dIn = sz(1);
    dOut = sz(1);
end

if nargin == 2
    dOut = dIn;
end

[U,S,V] = graph_svd(graph_laplacian(A));

subspaceIn = 1:dIn;
subspaceOut = 1:dOut;

if length(sz) == 2
    latIn = S.^.5 * V(:,subspaceIn);
    latOut = S.^.5 * U(:,subspaceOut);
    
    return;
end

latIn = zeros(sz(1),dIn,sz(3));
latOut = zeros(sz(1),dOut,sz(3));
for k=1:sz(3)
    latIn(:,:,k) = S(:,:,k).^.5 * V(:,subspaceIn,k);
    latOut(:,:,k) = S(:,:,k).^.5 * U(:,subspaceOut,k);
end
    

end

