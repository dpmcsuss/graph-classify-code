function ind = get_max_edges(delhat,Nmax)
% finds the constants.Nmax vertices with max degree
% if constants.Nmax is unspecified, n/10 is used

deg = sum(delhat,1) + sum(delhat,2)';
[degsort IX] = sort(deg,'descend');
n=length(deg);
if nargin==1, Nmax=round(n/10); end
ind_max_mat=zeros(n);
ind_max_mat(IX(1:Nmax),IX(1:Nmax))=1;
ind=find(ind_max_mat);