function ind = get_max_edges(delhat,num_vertices)
% finds the num_vertices vertices with max degree
% if num_vertices is unspecified, n/10 is used

deg = sum(delhat,1) + sum(delhat,2)';
[degsort IX] = sort(deg,'descend');
n=length(deg);
if nargin==1, num_vertices=round(n/10); end
ind_max_mat=zeros(n);
ind_max_mat(IX(1:num_vertices),IX(1:num_vertices))=1;
ind=find(ind_max_mat);