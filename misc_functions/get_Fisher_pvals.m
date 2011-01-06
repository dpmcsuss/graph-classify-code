function pval = get_Fisher_pvals(A0,A1,d,s0,s1)
% compute Fisher's exact test value for the two classes
% INPUT:
%   A0: array of s0 adjacency matrices, each n x n
%   A1: array of s1 adjacency matrices, each n x n
%   d:  n^2

siz=size(A0);
if length(siz)==3
    A0=reshape(A0,[d s0]);
    A1=reshape(A1,[d s1]);
else
    d=siz(1);
end

for i=1:d
    ss0=sum(A0(i,:));
    x(1,1)=ss0;
    x(1,2)=s0-ss0;
    
    ss1=sum(A1(i,:));
    x(2,1)=ss1;
    x(2,2)=s1-ss1;
    
    pvalue=FisherExactTest22(x);
    pval(i)=pvalue(3);
end

n=sqrt(d);
if length(siz)==3
    pval=reshape(pval,n,n);
end