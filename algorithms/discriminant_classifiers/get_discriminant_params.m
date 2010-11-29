function params = get_QDA_params(x,inds)



mu0=mean(x(:,inds.y0trn),2);
mu1=mean(x(:,inds.y1trn),2);

Sig0=cov(x(:,inds.y0trn)');
Sig1=cov(x(:,inds.y1trn)');

ytrn=[inds.y0trn inds.y1trn];
Sig=cov(x(:,ytrn)');

params.mu0=mu0;
params.mu1=mu1;
params.InvSig=inv(Sig);
params.InvSig0=inv(Sig0);
params.InvSig1=inv(Sig1);

s = length(ytrn);
params.lnprior0 = log(length(inds.y0trn)/s);
params.lnprior1 = log(length(inds.y1trn)/s);
