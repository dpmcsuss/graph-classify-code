function P = Bin_get_params(As,G,i)

if nargin==3
    y0 = G.y0; y0(y0==i)=[];
    y1 = G.y1; y1(y1==i)=[];
else
    y0=G.y0;
    y1=G.y1;
end

% deltahat
P.E0    = mean(As(:,:,y0),3);
P.E0(P.E0==0) = 1/G.s;
P.E0(P.E0==1) = 1-1/G.s;

P.E1    = mean(As(:,:,y1),3);
P.E1(P.E1==0) = 1/G.s;
P.E1(P.E1==1) = 1-1/G.s;

P.del   = abs(P.E0-P.E1);

P.lnE0  = log(P.E0);
P.ln1E0 = log(1-P.E0);
P.lnE1  = log(P.E1);
P.ln1E1 = log(1-P.E1);

% 2x2 block model
if isfield(G,'b')
    P.B0 = get_blocks(P.E0,G);
    P.B1 = get_blocks(P.E1,G);
end