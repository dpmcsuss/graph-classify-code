function [Atrn Gtrn Atst Gtst inds] = crossvalprep(As,constants,trn_s0,trn_s1)
% this function prepares data to be run by various classification
% algorithms
% INPUT:
%   As:         adjacency matrices
%   constants:  number of samples, etc.
%   trn_s0:     how many samples to use for class 0
%   trn_s1:     how many samples to use for class 1
% OUTPUT:
%   Atrn:       adjacency matrices for training
%   Gtrn:       constatnts for training data
%   Atst:       adj. mat.'s for testing
%   Gtst:       constants for testing data
%   inds:       collection of indices (for graph invariant approach)

tst_s0=constants.s0-trn_s0;
tst_s1=constants.s1-trn_s1;

ind0  = randperm(constants.s0);
ind1  = randperm(constants.s1);

y0trn = constants.y0(ind0(1:trn_s0));
y0tst = constants.y0(ind0(trn_s0+1:trn_s0+tst_s0));

y1trn = constants.y1(ind1(1:trn_s1));
y1tst = constants.y1(ind1(trn_s1+1:trn_s1+tst_s1));

Atrn = As(:,:,[y0trn y1trn]);
Atst = As(:,:,[y0tst y1tst]);

ytrn = [zeros(1,length(y0trn)) ones(1,length(y1trn))];
ytst = [zeros(1,length(y0tst)) ones(1,length(y1tst))];

Gtrn = get_constants(Atrn,ytrn);
Gtst = get_constants(Atst,ytst);

inds.y0trn = y0trn;
inds.y1trn = y1trn;
inds.y0tst = y0tst;
inds.y1tst = y1tst;
inds.ytst  = [inds.y0tst inds.y1tst];
inds.s_tst = length(inds.ytst);