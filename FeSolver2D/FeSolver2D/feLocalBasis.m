% p : Linear/Quadratic/Cubic

function coefmat = feLocalBasis(bsOdr)


  NumLx = 0;

  switch bsOdr

    case 'L'
        NumLx = 3;
        
        Onode = [-1,-1;...
                  1,-1;...
                 -1, 1];
             
        Lb = zeros(NumLx, NumLx);
        
        for s = 1:NumLx
            sx = Onode(s, 1);
            sy = Onode(s, 2);

            evalvec = [1 sx sy];
            Lb(s, :) = evalvec;
        end
%             Lb = [1, -1, -1;...
%                   1, 1, -1;...
%                   1, -1, 1];

        rhs = eye(NumLx);
    case 'Q'
        
        NumLx = 6;
        
        Onode = [-1,-1;...
                  1,-1;...
                 -1, 1;...
                  0,-1;...
                  0, 0;...
                 -1, 0];

        Lb = zeros(NumLx, NumLx);
        
        for s = 1:NumLx
            sx = Onode(s, 1);
            sy = Onode(s, 2);

            evalvec = [1 sx sy sx^2 sx*sy sy^2];
            Lb(s, :) = evalvec;
        end
        
%         Lb = [1, -1, -1, 1, 1, 1;...
%               1, 1, -1, 1, -1, 1;...
%               1, -1, 1, 1, -1, 1;...
%               1, 0, -1, 0, 0, 1;...
%               1, 0, 0, 0, 0, 0;...
%               1, -1, 0, 1, 0, 0];
        rhs = eye(NumLx);

    case 'C'
        NumLx = 10;
        stdpos = 1/3;
        Onode = [-1, -1;...
                  1, -1;...
                 -1, 1;...
                 -stdpos, -1;...
                 stdpos, -1;...
                 stdpos, -stdpos;...
                 -stdpos, stdpos;...
                 -1, stdpos;...
                 -1, -stdpos;...
                 -stdpos, -stdpos];

        Lb = zeros(NumLx, NumLx);
        
        for s = 1:NumLx
            sx = Onode(s, 1);
            sy = Onode(s, 2);

            evalvec = [1 sx sy sx^2 sx*sy sy^2 sx^3 sx^2*sy sx*sy^2 sy^3];
            Lb(s, :) = evalvec;
        end

        rhs = eye(NumLx);
  end

%   EIG = eig(Lb);
  iLb = inv(Lb);
  coefmat = iLb * rhs;
    
% TEST% %   sx = -1;
% % % % %   sy = -1;
% % % % %   textvec = [1 sx sy sx^2 sx*sy sy^2 sx^3 sx^2*sy sx*sy^2 sy^3]
% % % % %   aswer = textvec*coe
return              