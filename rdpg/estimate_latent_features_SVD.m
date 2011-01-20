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

% subspaceIn = 1:dIn;
% subspaceOut = 1:dOut;

if length(sz) == 2
    latIn =  V(:,1:dIn)*(S(1:dIn,1:dIn).^.5);
    latOut = U(:,1:dOut)*(S(1:dOut,1:dOut).^.5);
    
    return;
end

latIn = zeros(sz(1),dIn,sz(3));
latOut = zeros(sz(1),dOut,sz(3));
for k=1:sz(3)
     latIn(:,:,k) =  V(:,1:dIn,k)*(S(1:dIn,1:dIn,k).^.5);
     latOut(:,:,k) = U(:,1:dOut,k)*(S(1:dOut,1:dOut,k).^.5);
end
    

end

