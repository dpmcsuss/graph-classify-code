function constants = get_constants(adjacency_matrices,ys,directed,loopy)
% get constants for doing graph classification

constants.ys  = ys;             % class labels for each subject
constants.y0  = find(ys==0);    % list of class 0 subjects
constants.y1  = find(ys==1);    % list of class 1 subjects
siz   = size(adjacency_matrices);       % output is n x n x s
n   = siz(1);         % # vertices
constants.s0  = length(constants.y0);   % # of class 0 examples
constants.s1  = length(constants.y1);   % # of class 1 examples
constants.s   = siz(3);         % # training samples

dummy=ones(n);
if nargin==4,
    if directed && loopy
        d=n^2;
        constants.inds=1:n^2;
    elseif ~directed && loopy
        d=(n^2+n)/2;
        constants.inds=find(triu(dummy));
    elseif ~directed && ~loopy
        d=n*(n-1)/2;
        constants.inds=find(triu(dummy,+1));        
    elseif directed && ~loopy
        d=n*(n-1);
        constants.inds=1:n^2;
        constants.inds(1:n+1:n^2)=[];
    end
else % default is that graphs are directed and loopy
    d=n^2;
end

constants.d=d;
constants.p=2^d;
constants.n=n;