function [Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,s0_trn,s1_trn,s0_tst,s1_tst)
% this function prepares data for cross-validation using hold-out data
% 
% INPUT:
%   As:         adjacency matrices
%   constants:  number of samples, etc.
%   s0_trn:     # of training samples for class 0
%   s1_trn:     # of training samples for class 1
%   s0_tst:     (optional) # of testing samples for class 0
%   s1_tst:     (optional) # of testing samples for class 1
% 
% OUTPUT:
%   Atrn:       adjacency matrices for training
%   Gtrn:       constatnts for training data
%   Atst:       adj. mat.'s for testing
%   Gtst:       constants for testing data
%   inds:       collection of indices (for graph invariant approach)

if nargin<5, 
    s0_tst=constants.s0-s0_trn;
    s1_tst=constants.s1-s1_trn;
end

ind0  = randperm(constants.s0);
ind1  = randperm(constants.s1);

y0trn = constants.y0(ind0(1:s0_trn));
y0tst = constants.y0(ind0(s0_trn+1:s0_trn+s0_tst));

y1trn = constants.y1(ind1(1:s1_trn));
y1tst = constants.y1(ind1(s1_trn+1:s1_trn+s1_tst));

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