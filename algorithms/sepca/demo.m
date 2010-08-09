% data
tic;
load FlipRate010_Set01
N=size(X,2);

M = 1;  
d = 15;
maxiter.R = 1;
maxiter.LocMix = M;
maxiter.A = 99;
maxiter.YW = 5;
tol.Y = 0.05;
tol.W = 0.05;
tol.A = 0.03;
tol.Obj = 1;
tol.R = [.05 1]; 
THMAX_A = 100;
MAX_CUT_SEQ=3; 
fn_g   = @g_BernExp;
fn_dg  = @g_deriv_BernExp;
fn_d2g = @g_deriv2_BernExp;
fn_h   = @(x)(zeros(size(x))); % not used for Bern.distrib.
InitR = ones(1,N);

rand('twister', 1);
[Y W A B R P] = mixsepcatrain(...
	X,d,M, ...               % model para.
	fn_g, fn_dg, fn_d2g, ... % exp.fam.
	InitR, ...               % init.guess
	'MaxIters', maxiter, ... % nui.para.
	'Tol', tol, ...
	'Cut', THMAX_A, ...
	'InitA', 1, ...
	'InitB', 1, ...
	'VirFact', 20, ...
	'MaxCutSeq',MAX_CUT_SEQ);

Xr = W{1}*Y{1};
subplot(1,2,1); imagesc(X);  title('Data')
subplot(1,2,2); imagesc(Xr); title('Reconstruction')
colormap gray;
