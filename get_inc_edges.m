function ind = get_inc_edges(delhat,Ninc)
% finds the top G.Ninc edges.  
% if G.Ninc isn't specified, it is set to # of vertices

if nargin==1, Ninc=numel(delhat)/10; end 
[delhatsort delind]=sort(delhat(:),'descend');
ind  = delind(1:Ninc);