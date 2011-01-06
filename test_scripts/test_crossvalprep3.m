clear, clc

n = 10;
m = 3;
s = 1000;

p = 0.5;
q = 0.35;
q0 = p-q;
q1 = p+q;

E0 = p.*ones(n); E0(1:m,1:m) = q0;
E1 = p.*ones(n); E1(1:m,1:m) = q1;

As = nan(n,n,s);
As(:,:,1:s/2)=repmat(E0,[1 1 s/2]) > rand(n,n,s/2);
As(:,:,s/2+1:s)=repmat(E1,[1 1 s/2]) > rand(n,n,s/2);

targs=[zeros(1,s/2) ones(1,s/2)];

Lstar = 1-binocdf(floor(m^2/2),m^2,p-q)

constants = get_constants(As,targs);

%%

xval.s0_trn=round(constants.s0*0.9);
xval.s1_trn=round(constants.s1*0.9);
xval.num_iters=100;

subspace(1).name='manual';
subspace(1).indices=find(E0-E1);
subspace(1).size=length(subspace(1).indices);
Lhats = classify_par(As,constants,xval,subspace);
mean(Lhats)

% xval.num_repeats=100;
% subspace.maual=find(E0-E1);
% Lhats = classify_holdout(sim{i}.As,constants,xval,subspace);
% mean(Lhats.maual)
