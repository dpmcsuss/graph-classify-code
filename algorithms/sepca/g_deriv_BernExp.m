function dg = g_deriv_BernExp(theta)
%
% g_deriv_BernExp.m
% The derivative of the g-functional to the input, where the g-functional is
% the one in exp. family form for Bernoulli distribution (see also g_BernExp.m) 

dg = - 1./(1+exp(-theta));
