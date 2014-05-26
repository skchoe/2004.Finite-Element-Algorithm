% function [SolInside, geomInside, SolOuterBD, geomOuterBD, SolInnerBD, ...
%	  geomInnerBD] = tmpDivUnknowns(NsolIn, NsolBd, BdyDP, v, geom)
% this function separates the inner boundary, outer boundary and 
% inside points into different vectors 
% it also divides up the geometric location vector into the above
% mentioned categories


function [SolInside, geomInside, SolOuterBD, geomOuterBD, SolInnerBD, ...
	  geomInnerBD] = tmpDivUnknowns(Ninside, Nouter, boundaryConds, ...
					solution, geom) 

SolInside = sparse(Ninside, 1);
SolOuterBD = sparse(Nouter, 1);
SolInnerBD = boundaryConds;

%take care of the inside points
for i = 1:Ninside
  SolInside(i, 1) = solution(i, 1);
  geomInside(i,1:2) = geom(i,2:3);
end

%take care of the outer boundary points
for i = Ninside+1:Ninside+Nouter
  SolOuterBD(i-Ninside, 1) = solution(i, 1);
  geomOuterBD(i-Ninside, 1:2) = geom(i,2:3);
end

%take care of the inner boundary points
for i = Ninside+Nouter+1:length(geom)
  geomInnerBD(i-Ninside-Nouter,1:2) = geom(i,2:3);
end
