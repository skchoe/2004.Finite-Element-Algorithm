% [dgq, dgrq, dglq] = QuadratureDegree(udegree)

function [dgq, dgrq, dglq] = QuadratureDegree(udegree)

fltQ = (udegree+1)/2;

intQ = round(fltQ);
if intQ < fltQ
  dgq = intQ + 1;
else
  dgq = intQ;
end

if dgq < (udegree+2)/2
  dgrq = dgq + 1;
else
  dgrq = dgq;
end

if dgrq < (udegree+3)/2
  dglq = dgrq + 1;
else
  dglq = dgrq;
end

return