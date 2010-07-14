function constants = get_constants(adjacency_matrices,ys)
% get constants for doing graph classification

constants.ys  = ys;             % class labels for each subject
constants.y0  = find(ys==0);    % list of class 0 subjects
constants.y1  = find(ys==1);    % list of class 1 subjects
siz   = size(adjacency_matrices);       % output is n x n x s
constants.n   = siz(1);         % # vertices
constants.s0  = length(constants.y0);   % # of class 0 examples
constants.s1  = length(constants.y1);   % # of class 1 examples
constants.s   = siz(3);         % # training samples