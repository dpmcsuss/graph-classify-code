function [ L ] = graph_laplacian( A )
%graph_laplacian Comput scaled sign positive graph laplacian
%   compute the scaled sign positive graph laplacian of each graph (call it
%   L_i).  one obtains this by filling in the diagonal of each adjacency
%   matrix by the degree of that vertex, divided by n-1 (this is the 
%   expected weight of the self loop).
d = size(A);
L=A;
if numel(d)==2
for k=1:d(1)
    L(k,k) = sum(A(k,:))/(d(1)-1);
end
end
if numel(d) == 3
   for g=1:d(3)
       for k=1:d(1)
           L(k,k,g) = sum(A(k,:,g))/(d(1)-1);
       end
   end   
end
end

