function sim = gen_ind_edge_data(M)
% takes as input some parameters, and sample data and computes Lstar

n = M.n;
s = M.s;
E0 = M.E0;
E1 = M.E1;

As = nan(n,n,s);
As(:,:,1:s/2)=repmat(M.E0,[1 1 s/2]) > rand(n,n,s/2);
As(:,:,s/2+1:s)=repmat(M.E1,[1 1 s/2]) > rand(n,n,s/2);

sim.As=As;
sim.targs=[zeros(1,s/2) ones(1,s/2)];

sim.params.E0=E0;
sim.params.E1=E1;

sim.params.lnE0  = log(E0);
sim.params.ln1E0 = log(1-E0);
sim.params.lnE1  = log(E1);
sim.params.ln1E1 = log(1-E1);

sim.params.lnprior0 = log(1/2);
sim.params.lnprior1 = log(1/2);

sim.params.d_pos = abs(E0-E1);           % position difference
sim.params.d_opt = abs(E0./sqrt(E0.*(1-E0)) - E1./sqrt(E1.*(1-E1))); % optimal difference