function d_pval = run_get_fisher_pvals(As,constants)

d_pval = get_Fisher_pvals(...
    As(:,:,constants.y0),...    % class 0 graphs
    As(:,:,constants.y1),...    % class 1 graphs
    constants.d,...             % # of edges
    constants.s0,...            % # of class 0 samples
    constants.s1);              % # of class 1 samples
