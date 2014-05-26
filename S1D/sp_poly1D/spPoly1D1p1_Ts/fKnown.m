% fKnown.m
% Sarah Geneser
% 02-09-04
%
% forcing function for the example problem

function [f] = fKnown(x)
f = -2*pi*(exp(x).*sin(2*pi*x) + 2*pi*exp(x).*cos(2*pi*x));
