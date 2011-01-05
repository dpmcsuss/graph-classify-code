clear, clc,

s=100;
s0 = s/2;
s1 = s-s0;

k=5;

mu0=zeros(1,k);
mu1=0.1*ones(1,k); %[1 1 4 0 0];
Sig=eye(k);
c=0.5;
Sig(2:k+1:k^2)=c*ones(4,1);
Sig(k+1:k+1:k^2)=c*ones(4,1);

d=0.25;
Sig(3:k+1:4*k)=d*ones(3,1);
Sig(2*k+1:k+1:5*k)=d*ones(3,1);

Sig0=Sig;

imagesc(Sig0)
L=chol(Sig0);

Sig1=eye(k);

inds.y0trn  = 1:s0/2;
inds.y1trn  = s0+1:s-s1/2;
inds.y0tst  = s0/2+1:s0;
inds.y1tst  = s0+s1/2+1:s;
inds.ytst   = [inds.y0tst inds.y1tst];
inds.s_tst  = length(inds.ytst);


x=randn(k,s);
for i=1:s
    if any(i==[inds.y0trn inds.y0tst])
        x(:,i)=x(:,i)'*L;
    else
        x(:,i)=x(:,i)'*Sig1;
    end
end

x(:,s0+1:s)=x(:,s0+1:s)+repmat(mu1',1,s1);
y=[zeros(s0,1); ones(s1,1)];



% % use correct params
params.mu0=mu0';
params.mu1=mu1';
params.InvSig=inv(Sig);
params.InvSig0=inv(L);
params.InvSig1=eye(k);



params.lnprior0 = log(s0/s);
params.lnprior1 = log(s1/s);

%%
discrim.LDA=[];
discrim.QDA=[];


Lhat_true = discriminant_classifiers(x,y',params,discrim)


% estimate params
est_params = get_QDA_params(x,inds);
Lhat_est = discriminant_classifiers(x,y',est_params,discrim)
