function [ Pvalue ConTable ]=FisherExactTest22( x,y )
% Fisher's Exact Probability Test for contigency table 2*2.
% Explanation: http://mathworld.wolfram.com/FishersExactTest.html
% Also refer to: http://www.mathworks.com/matlabcentral/fileexchange/15434
% Inputs:
%   Choice (1) one variable  x - 2x2 data matrix 
%   Choice (2) two variable  x, y for binary vector with the same length.(No missing value)
% Outputs:
%   (1) four p-values organized as: 
%   Pvalue = [ left tail , Right tail, 2-tails, Mid-p correction for 2-tails ] 
% Created by Lowell Guangdi 2010/01/28

if nargin > 1
   ConTable = zeros(2);
   for p = 1:length(x)
       if x( p ) == 0 && y( p ) == 0, ConTable(1,1) = ConTable(1,1) + 1;
       elseif x( p ) == 0 && y( p ) == 1 ,ConTable(1,2) = ConTable(1,2) + 1;
       elseif x( p ) == 1 && y( p ) == 0 ,ConTable(2,1) = ConTable(2,1) + 1;
       else ConTable(2,2) = ConTable(2,2) + 1;
       end
   end
   x = ConTable;
else
   ConTable = x; 
end
Rs=sum(x,2); Cs=sum(x); N=sum(Rs); 
%Rearrange the matrix if necessary.
if ~issorted(Rs),x=flipud(x); Rs=sort(Rs); end
if ~issorted(Cs),x=fliplr(x); Cs=sort(Cs); end
A=0:1:min(Rs(1),Cs(1)); 
z=[A;Rs(2)-Cs(1)+A;Rs(1)-A;Cs(1)-A;];
np=zeros(1,length(A)); lz=log(z);
np(1)=sum(gammaln([Rs(2)+1 Cs(2)+1])-gammaln([N+1 z(2)+1]));
f=sum(lz(3:4,1:end-1))-sum(lz(1:2,2:end));
np(2:end)=np(1)+cumsum(f);
np=exp(np); W=x(1)+1;
%now compute the 1-tailed p-values and the 2-tailed p-value
if x(1)<=round(A(end)/2) %choose direction
    P=[sum(np(1:W)) sum(np(W:end)) sum(np(np<=np(W)))];
else
    P=[sum(np(W:end)) sum(np(1:W)) sum(np(np<=np(W)))];
end
Pvalue = [ P 0.5*np(W)+sum(np(np<np(W))) ];
end