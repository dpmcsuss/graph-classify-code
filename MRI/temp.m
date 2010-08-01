a=poissrnd(1e6,1);
for t=1:1000
    any(a>1);
    all(a<=1);
end