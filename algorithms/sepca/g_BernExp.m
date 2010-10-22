function g = g_BernExp(theta)
%
% g_BernExp.m
% The g-functional in exp. family form for Bernoulli distribution.

g = -log(1+exp(theta));
