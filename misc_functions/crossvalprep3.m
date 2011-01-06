function [Atrn Gtrn Atst Gtst inds] = crossvalprep3(As,constants,xval,i)
% this function prepares data for cross-validation using either:
% (1) hold-out data, or
% (2) using all permutations of leave-two-out (one from each class)
% depending on whether num_iters is sufficiently small to do all possible
% permuations
%
% INPUT:
%   As:         adjacency matrices
%   constants:  number of samples, etc.
%   xval:       structure containing xval parameters including
%       s0_trn: # of training samples for class 0
%       s1_trn: # of training samples for class 1
%       num_iters: (optional) max # of iters to perform
%       s0_tst: (optional) # of testing samples for class 0
%       s1_tst: (optional) # of testing samples for class 1
%
% OUTPUT:
%   Atrn:       adjacency matrices for training
%   Gtrn:       constatnts for training data
%   Atst:       adj. mat.'s for testing
%   Gtst:       constants for testing data
%   inds:       collection of indices (for graph invariant approach)


% max_iters=choose(constants.s0,xval.s0_trn)*choose(constants.s1,xval.s1_trn);
if ~isfield(xval,'num_iters'),
    xval.num_iters=constants.s; %max_iters;
end

if xval.num_iters==constants.s % leave-one-out
    y0trn=constants.y0;
    y1trn=constants.y1;
    if i<=constants.s0
        y0trn(i)=[];
        y0tst=i;
        y1tst=[];
    else
        y1trn(i-constants.s0)=[];
        y1tst=i-constants.s0;
        y0tst=[];
    end
elseif xval.num_iters==choose(constants.s0,xval.s0_trn)*choose(constants.s1,xval.s1_trn) && nargin==4 % leave-two-out
    
    smax=max(constants.s0,constants.s1);
    
    ind0tst=ceil(i/smax);
    ind1tst=mod(i,smax); if ind1tst==0, ind1tst=smax; end
    
    y0tst=constants.y0(ind0tst);
    y1tst=constants.y1(ind1tst);
    
    y0trn=constants.y0; y0trn(ind0tst)=[];
    y1trn=constants.y1; y1trn(ind1tst)=[];
else
    if ~isfield(xval,'s0_tst'), xval.s0_tst=constants.s0-xval.s0_trn; end
    if ~isfield(xval,'s1_tst'), xval.s1_tst=constants.s1-xval.s1_trn; end
    
    ind0  = randperm(constants.s0);
    ind1  = randperm(constants.s1);
    
    y0trn = constants.y0(ind0(1:xval.s0_trn));
    y0tst = constants.y0(ind0(xval.s0_trn+1:xval.s0_trn+xval.s0_tst));
    
    y1trn = constants.y1(ind1(1:xval.s1_trn));
    y1tst = constants.y1(ind1(xval.s1_trn+1:xval.s1_trn+xval.s1_tst));
    
end

try inds.ytrn=vertcat(y0trn,y1trn);
catch exception, inds.ytrn=horzcat(y0trn,y1trn); end
try inds.ytst=vertcat(y0tst,y1tst);
catch exception, inds.ytst=horzcat(y0tst,y1tst); end


Atrn = As(:,:,inds.ytrn);
Atst = As(:,:,inds.ytst);

ytrn = [zeros(1,length(y0trn)) ones(1,length(y1trn))];
ytst = [zeros(1,length(y0tst)) ones(1,length(y1tst))];

Gtrn = get_constants(Atrn,ytrn);
Gtst = get_constants(Atst,ytst);

inds.y0trn = y0trn;
inds.y1trn = y1trn;
inds.y0tst = y0tst;
inds.y1tst = y1tst;
inds.s_tst = length(inds.ytst);