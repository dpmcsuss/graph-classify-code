function G = get_constants(As,ys)
% get constants for doing graph classification

G.ys  = ys;             % class labels for each subject
G.y0  = find(ys==0);    % list of class 0 subjects
G.y1  = find(ys==1);    % list of class 1 subjects
siz   = size(As);       % output is n x n x s
G.n   = siz(1);         % # vertices
G.s0  = length(G.y0);   % # of class 0 examples
G.s1  = length(G.y1);   % # of class 1 examples
G.s   = siz(3);         % # training samples