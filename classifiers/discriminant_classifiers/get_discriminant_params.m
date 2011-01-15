function params = get_discriminant_params(x,inds,discrim)



mu0=mean(x(:,inds.y0trn),2);
mu1=mean(x(:,inds.y1trn),2);
siz=size(x);
diag_ind=1:siz(1)+1:(siz(1))^2;


if isfield(discrim,'LDA') || isfield(discrim,'dLDA')
    Sig=cov(x(:,inds.ytrn)');
end

if isfield(discrim,'LDA')
    params.InvSig=inv(Sig);
end

if isfield(discrim,'dLDA')
    InvdSig=zeros(siz(1));
    InvdSig(diag_ind)=diag(Sig).^-1;
    params.InvdSig=InvdSig;
    if any(isinf(params.InvdSig))
        error('something is constant that should not be');
    end
end


% if isfield(discrim,'dLDA');
%     Sig=cov(x(:,ytrn)');
%     params.InvSig=inv(Sig);
% end



if isfield(discrim,'QDA')  || isfield(discrim,'dQDA')
    Sig0=cov(x(:,inds.y0trn)');
    Sig1=cov(x(:,inds.y1trn)');
    if isfield(discrim,'QDA') 
        params.InvSig0=inv(Sig0);
        params.InvSig1=inv(Sig1);
    end
    if isfield( discrim,'dQDA' )
        params.InvdSig0 = diag(diag(Sig0).^-1);
        params.InvdSig1 = diag(diag(Sig1).^-1);
    end
end


params.mu0=mu0;
params.mu1=mu1;

s = length(inds.ytrn);
params.lnprior0 = log(length(inds.y0trn)/s);
params.lnprior1 = log(length(inds.y1trn)/s);
