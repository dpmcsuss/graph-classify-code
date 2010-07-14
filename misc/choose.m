function nck=choose(n,k)
%assumes k~n

nck=1;
if k==n
    return;
end
i=0;
while i<n-k
    nck=nck*(n-i);
    i=i+1;
end
nck=nck/factorial(round(i));