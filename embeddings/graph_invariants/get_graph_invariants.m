function x = get_graph_invariants(As)
% this function computes graph invariants for all adjacency matrices in As
% As is an array of s adjacency matrices, each n x n
% x is a p x s matrix,  where p is the number of graph invariants computed
% per graph

siz=size(As);       
x=zeros(5,siz(3));          % pre-allocate memory

for i=1:siz(3);
    A=As(:,:,i);            % for each adjacency matrix
    x(1,i)=sum(A(:));       % size
    x(2,i)=max(sum(A));     % max-degree
    x(3,i)=maxeval(A);      % max eigenvalue (upper bound on max average degree)
    x(4,i)=scan1(A,siz(2)); % scan statistic 1
    x(5,i)=numtri(A);       % number of triangles
    x(6,i)=clustcoef(A,siz(2)); % clustering coefficient
    x(7,i)=apl(A,siz(2));   % average path length
%     x(8,i)=scan2(A);        % scan statistic 2
%     x(9,i)=scan3(A);        % scan statistic 3
end

x(5,:)=x(5,:)/siz(1);

vars=var(x,[],2);
rm_var=find(vars<1e-4);
x(rm_var,:)=[];

if ~isempty(rm_var), 
    disp(['no variance in invars ' num2str(rm_var') ', so i removed them']); 
end

end

%% get number of trianlges
function tri=numtri(A)

S=sparse(A);
T=trace(S^3)/6;
tri=T(1,1);

end


%% returns the maximum number of edges in a induced neighborhood (k=1) subgraph
function T=scan1(G,n)

T=0;
for i=1:n
    ind=find(G(:,i));
    Gsub=G([ind;i],[ind;i]);
    tem=sum(sum(Gsub))/2;
    if tem>T
        T=tem;
    end
end

end

%% get max eval
function T=maxeval(A)

OPTS.disp = 0;
T=eigs(double(A),1,'lm',OPTS);

end

%% get clustering coefficient
function ci = clustcoef(A,n)

ci = zeros(1,n);
for k = 1:n
    neighbours = find(A(:,k))';
    neighbours(neighbours==k) = []; % self link deleted
    a = length(neighbours);
    if a == 0; ci(k) = 0; continue; end
    if a < 2, continue; end
    E = 0;
    for ks = 1:a
        k1 = neighbours(ks);
        for k2 = neighbours(ks+1:a)
            if A(k1,k2) || A(k2,k1)
                E = E + 1;
            end
        end
    end
    ci(k) = 2 * E / (a * (a-1));
end
ci = mean(ci);

end        


%% get average path length
function x = apl(B,num_vertices)

if nargin==1, num_vertices=numel(B); end % assume directed graph (doesn't matter much though probably
B(B==0)=Inf;
C=ones(size(B));
while any(C(:))
    C=B;
    B=min(B,squeeze(min(repmat(B,[1 1 length(B)])+...
        repmat(permute(B,[1 3 2]),[1 length(B) 1]),[],1)));
    C=B-C;
end
B(logical(eye(length(B))))=0;
B=B(:);
B(isinf(B))=[];
x=sum(B)/num_vertices;

end







        
