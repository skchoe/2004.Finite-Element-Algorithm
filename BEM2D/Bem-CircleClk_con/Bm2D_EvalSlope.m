%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spEvalAvSlope.m 
% Computation of average slope by given curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = Bm2D_EvalSlope(domvec, imgvec)

  Nd = size(domvec, 1)
  Ni = size(imgvec, 1)

  N = min(Nd, Ni);
    if N <=0 N=1; end
  db = conv2log10(domvec(N, 1));
  da = conv2log10(domvec(1, 1));

  ib = conv2log10(imgvec(N, 1));
  ia = conv2log10(imgvec(1, 1));

  diffdom = db - da;
  diffimg = ib - ia;
  

  s = diffimg/diffdom;

return

function lg = conv2log10(val)

  if val < 1e-15
    lg = -log10(abs(val));
  else
    lg = log10(abs(val));
  end
return



% Restriction: N1 > N2 > N3 - relative slope
% function [s1, s2, s3] = spEvalSlope(imgvec1, imgvec2, imgvec3)
% 
% imgvec1
% imgvec2
% imgvec3
%   N1 = size(imgvec1, 1)
%   N2 = size(imgvec2, 1)
%   N3 = size(imgvec3, 1)
%   
%   s1 = 1;
%   
%   sum1 = 0;
%   sum2 = 0;
%   for i = 1:N2-1
%     d1 = imgvec1(i+1, 1) - imgvec1(i, 1);
%     sum1 = sum1 + d1;
%     d2 = imgvec2(i+1, 1) - imgvec2(i, 1);
%     sum2 = sum2 + d2;
%     
%   end
%   
%   s2 = sum2 / sum1;
% 
%   sum1 = 0;
%   sum3 = 0;
%   for i = 1:N2-1
%     d1 = imgvec1(i+1, 1) - imgvec1(i, 1)
%     sum1 = sum1 + d1;
%     d3 = imgvec3(i+1, 1) - imgvec3(i, 1)
%     sum3 = sum3 + d3;
%   end
%   
%   s3 = sum3 / sum1;
% 
% return

% This module is to compute slope of imgvec w.r.t domvec
% function slope = spEvalSlope(domvec, imgvec)
% 
% length = size(domvec,1);
% lengthtmp = size(imgvec, 1);
% if length ~= lengthtmp
%   exit;
% else
%     
%   slopeset = size(length-1, 1);
%   slopesum = 0;
%   counter = length - 1;
%   
%   for i = 1:length-1
%     nom = domvec(i+1, 1) - domvec(i, 1)
%     denom = imgvec(i+1, 1) - imgvec(i, 1)
%     if nom > 1e-15
%       slopesum = slopesum + denom / nom
% %      slopeset(i, 1) = slopesum;
%     else 
%       counter = counter - 1;
%     end
%   end
%     
%   slope = slopesum/counter
% end
% return
