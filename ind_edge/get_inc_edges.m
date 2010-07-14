function ind = get_inc_edges(delhat,num_edges)
% finds the top num_edges edges  
% if num_edges isn't specified, it is set to # of vertices

if nargin==1, num_edges=sqrt(numel(delhat)); end 
[delhatsort delind]=sort(delhat(:),'descend');
ind  = delind(1:num_edges);