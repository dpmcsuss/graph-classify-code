function [Y,W,A,B,R,P]=mixsepcatrain(...
	X,latd,M,...
	fn_g,fn_dg,fn_d2g,...
	InitR,...
	varargin)
%-------------------------------------------------------
% [Y,W,A,R,P]=mixsepcatrain(X,d,M,fn_g,fn_dg,fn_d2g,varargin)
%
% Train mixture of SePCA model
% 
% -- input --
% X             [D x N] data matrix
% M             #.SePCA models in the mixture
% fn_g          func [g] in the model
% fn_dg         1st gradient of [g]
% fn_d2g        2nd gradient of [g]
% -- output --
% Y             [M]-cell: [q x N] low dimensional representation of data
% W             [M]-cell: [D x q] PC axes
% A             [M]-cell: [q] array, the alpha's
% B             [M], beta of each loacl SePCA (dummy)
% R             [M x N], membership
% P             [M], prior of each local SePCA
%-------------------------------------------------------

% Created Jun Li, Oct, 2009
% Modified Jun Li, 8 Aug, 2010
% - revise comments

% todo: use new parameter parser and checker

% house keeping
i = 1;
while(i<length(varargin))
	switch upper(varargin{i})
	case 'REPORT'
		op_Report = varargin{i+1};
  case 'INITA'
    op_InitA = varargin{i+1};
  case 'INITB'
    op_InitB = varargin{i+1};
	case 'MAXITERS'
		op_nMaxIters = varargin{i+1};
	case 'TOL'
		op_Tol = varargin{i+1};
	case 'CUT'
		op_Cut = varargin{i+1};
	case 'VIRFACT'
		op_VirFact = varargin{i+1};
	case 'MAXCUTSEQ'
		op_nMaxCutSeq = varargin{i+1};
	end
	i=i+2;
end
if ~exist('op_InitB','var')
	op_InitB = 1;
end

% consts and init
optmopt_fastW = optimset('Display', 'notify', 'GradObj', 'on', 'Hessian','on', 'TolFun', 0.05, 'TolX', op_Tol.W);
optmopt_fastY = optimset('Display', 'notify', 'GradObj', 'on', 'Hessian','on', 'TolFun', 0.1, 'TolX', op_Tol.Y);
[D,N] = size(X);
Y = cell(1,M);
W = cell(1,M);
A = cell(1,M);
for m=1:M
	Y{m} = zeros(latd,N);
	W{m} = rand(D,latd); 
	A{m} = ones(latd,1)*op_InitA;
	Y{m}(1,:) = 1;
	A{m}(1) = 0;
end
R = InitR;
d = ones(1,M)*latd;          % local dim.
B = ones(1,M)*op_InitB;      % local beta.
P = sum(R,2)./sum(R(:));     % mix.prop.
% records
ObjBase = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,op_VirFact);
MemberShipRecord=cell(1,op_nMaxIters.R+1);
MemberShipRecord{1} = R;

% EM-Loop for R and others
cIter_R  =0;
bCnvg_R  = false;
ObjVal_R = NaN;
while  cIter_R<op_nMaxIters.R && ~bCnvg_R
	ObjVal_OldR = ObjVal_R;
	
	% M-step: Pi's
	P = sum(R,2)./sum(R(:));

	% M-step: A, W_MP, Y_MP
	for m=1:M % for each local model
		r_m = R(m,:)./sum(R(m,:)); % n-vec, norm. wei. for each point, *diff* from R_m
		Y{m}(2:end,:)=0;  % each time of membership realloc. redo the mixtures

		% Optimise a local model
		cIter_LocMix = 0; 
		bCnvg_LocMix = false;
		ObjVal_LocMix = NaN;
		while cIter_LocMix<op_nMaxIters.LocMix && ~bCnvg_LocMix
			nCutInARow = 0;
			ObjVal_OldLocMix = ObjVal_LocMix;
			j=1;
			while j<=d(m) % for each latent dim
				W{m}(:,j) = rand(D,1);

				cIter_A  = 0;
				bCnvg_A  = false;
				ObjVal_A = NaN;
				while cIter_A<op_nMaxIters.A && ~bCnvg_A
					bCut=false;
					ObjVal_OldA = ObjVal_A;

					cIter_YW  = 0;
					bCnvg_YW  = false;
					ObjVal_YW = NaN;
					while cIter_YW<op_nMaxIters.YW && ~bCnvg_YW
						% save the old values
						y_j = Y{m}(j,:);
						w_j = W{m}(:,j);
						ObjVal_OldYW = ObjVal_YW;
						% optimise y_j, it takes the same procedure as in a local model
						if j>1 % add a mean
							D1 = W{m}'*X;
							Y{m}(j,:) = fminunc(@fobjY_bat,Y{m}(j,:),optmopt_fastY);
						else
							Y{m}(j,:) = 1;
						end
						% optimise w_j
						if j>1
							W{m}(:,j) = fminunc(@fobjW_bat,W{m}(:,j),optmopt_fastW);
						else
							W{m}(:,j) = fminunc(@fobjMu_bat,W{m}(:,j),optmopt_fastW);
						end
						% eval improv. and judge cnvg.
						ObjVal_YW = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,op_VirFact);
						bCnvg_YW = all(abs(y_j-Y{m}(j,:))<op_Tol.Y);
						bCnvg_YW = bCnvg_YW || (all(abs(w_j-W{m}(:,j))<op_Tol.W));
						bCnvg_YW = bCnvg_YW || (ObjVal_OldYW - ObjVal_YW < op_Tol.Obj);
						cIter_YW = cIter_YW+1;
					end % YW-Loop
					
					if j>1
						if W{m}(:,j)'*W{m}(:,j) < D/op_Cut  % cutoff? <=> a_j>xxx, but num. stab.
							bCut = true;
							bCnvg_A = true;
						else
							old_A = A{m}(j);
							A{m}(j)  = D/(W{m}(:,j)'*W{m}(:,j));
							ObjVal_A = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,op_VirFact);
							bCnvg_A = (ObjVal_OldA - ObjVal_A < op_Tol.Obj) || (abs(old_A-A{m}(j))<op_Tol.A && cIter_A>0);
						end % IF update A or cut-off
					else
						bCnvg_A = true;
						ObjVal_A = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,op_VirFact);
					end
					
					cIter_A = cIter_A+1;
					
					fprintf(1,'[%8.2f]mx %d(%d):dim %2d:A(%2d) %9.6f:E %12.2f\n', ...
						toc,m,cIter_LocMix+1, j, cIter_A, ...
						max(A{m}(j), inf*bCut), ...
						ObjVal_A-ObjBase);
				end % WLOOP: A 

				if bCut
					Old_A = A;
					Old_W = W;
					Old_Y = Y;
					Old_d = d;
					
					nCutInARow = nCutInARow+1;
					if nCutInARow < op_nMaxCutSeq
						vid=[1:j-1, j+1:d(m)];
						A{m}=A{m}(vid);
						W{m}=W{m}(:,vid); % discard cut-off axes
						Y{m}=Y{m}(vid,:); % discard cut-off coeff
						j=j-1;
						d(m)=d(m)-1;
					else
						vid=[1:j-1];						
						A{m}=A{m}(vid);
						W{m}=W{m}(:,vid); % discard cut-off axes
						Y{m}=Y{m}(vid,:); % discard cut-off coeff
						d(m)=j-1;
						j=j-1;
					end
				else
					nCutInARow = 0;
				end % IF: cut off happened
				
				% .. some report ..
				fprintf('[I %d:M %d]this dim %d, mix dim %d\n', cIter_R,m,j,d(m));
				% .. repo.done ..
				j=j+1;
			end % FLOOP: j-dim 

			ObjVal_LocMix = ObjVal_A; % the lastest objval doen't change after A-Loop
			bCnvg_LocMix  = (ObjVal_OldLocMix-ObjVal_LocMix < op_Tol.Obj);
			cIter_LocMix  = cIter_LocMix + 1;
		end % WLOOP: m-mixture EM

	end % FLOOP: each local mixture

	% E-step: R_{n,m}'s
	ROld = R;
	LocLogProb = zeros(M,N);
	for m=1:M
		LocLogProb(m,:) = LocMixLkhd(X,Y{m},W{m},fn_g);
	end
	for m=1:M
		% 1/(1+exp(L_1-L_i)+...+exp(L_N-L_i))
		R(m,:) = 1./(sum(exp(LocLogProb-repmat(LocLogProb(m,:),M,1)),1));
	end
	
	% .. record ..
	MemberShipRecord{cIter_R+1} = R;
	if exist('op_Report', 'var')
		tmpfname = sprintf('%s_I%02d',op_Report.Perfix,cIter_R+1);
		save(tmpfname,'Y','W','A','R','P','d');
	end
	% .. reco.done ..
	
	ObjVal_R = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,op_VirFact);
	bCnvg_R = (ObjVal_OldR - ObjVal_R < op_Tol.Obj);
	bCnvg_R = bCnvg_R || nnz(abs(ROld-R)>=op_Tol.R(1))<op_Tol.R(2);
	cIter_R = cIter_R+1;
end % WLOOP: R EM 

function [E Gr Gr2] = fobjY_bat(t)
	% obj val
	Y2 = Y{m}; Y2(j,:)=t;
	T  = W{m}*Y2;
	L = sum(T.*X+fn_g(T),1) - .5*t.*t*B(m);
	E = -sum(L);
	% 1st deri
	d1 = D1(j,:);
	d2 = W{m}(:,j)'*fn_dg(T);
	dY = d1 + d2 - t*B(m);
	Gr = -dY;
	% 2nd deri
	d2 = (W{m}(:,j).*W{m}(:,j))'*fn_d2g(T);
	dY2 = d2-B(m);
	Gr2 = spdiags(-dY2', 0, N,N);
end % FUNC. batch optimise Y obj

function [E, Gr, Gr2]=fobjW_bat(t)
W2 = W{m}; W2(:,j) = t;
vf = op_VirFact; % virtually data points (for enough supp.)
% obj val
T = W2*Y{m};
L = sum(r_m .* sum( T.*X+fn_g(T), 1) )*vf -.5*A{m}(j)*t'*t;
% 1st deri
w_yj = Y{m}(j,:).*r_m;
dW = (X+fn_dg(T))*w_yj'*vf - A{m}(j)*t;
% 2nd deri
dg_T2 = fn_d2g(T) * (Y{m}(j,:) .* w_yj)';
dW2 = dg_T2*vf-A{m}(j);

E  = -L;
Gr = -dW;
Gr2 = spdiags(-dW2, 0, D,D);
end

function [E, Gr, Gr2]=fobjMu_bat(t)
T = repmat(t,1,N);
vf = op_VirFact; % virtually data points (for enough supp.)
% obj val
L = sum(r_m .* sum( T.*X+fn_g(T), 1) )*vf;
% 1st deri
dT = (X+fn_dg(T))*r_m'*vf;
% 2nd deri
dT2 = fn_d2g(T) * r_m';

E  = -L;
Gr = -dT;
Gr2 = spdiags(-dT2, 0, D,D);
end

end % main function

function L = NegMargLkhd(X,Y,W,A,P,R,d,M,fn_g,vf)
D = size(X,1);
l2pi = log(2*pi);
L_mx = zeros(1,M);
for m=1:M % for each mixture, compute the log-likehood for each data point
	if ~all(A{m}(2:end)>0)
		L_mx(m)=NaN;
		continue;
	end
	T = W{m}*Y{m};
	L_x_Wy = sum(T.*X +fn_g(T),1);
	L_y    = -.5*sum(Y{m}.*Y{m},1) - d(m)/2*l2pi;
	L_W = -.5*sum(sum((W{m}.*W{m})*diag(A{m}),1)) ...
		- .5*D*d(m)*l2pi ...
		+ .5*D*sum(log(A{m}(2:end)));
	L_P = log(P(m));
	L_mx(m) = sum(R(m,:).*(L_x_Wy + L_y))*vf + sum(R(m,:))*(L_W+L_P);
end % for each mixture
L = - sum(L_mx);
end % func MargLkhd 

function L = LocMixLkhd(X,Y,W,fn_g) % Use log-likelihood instead
T = W*Y;
L_x_Wy = sum(T.*X +fn_g(T),1);
L = L_x_Wy; 
end % FUNC LocMixLkhd
