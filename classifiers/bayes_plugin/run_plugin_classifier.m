function [Lstats params subspace] = run_plugin_classifier(Atrn,Gtrn,Atst,Gtst,subspace)
% runs plugin classifiers using the various prespecified subspaces

params = get_ind_edge_params(Atrn,Gtrn);

subspace = get_subspace_indices(params,subspace,Atrn,Gtrn);

Lstats = plugin_classify(Atst,Gtst.ys,params,subspace);

end