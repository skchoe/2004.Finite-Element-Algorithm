


% evaluate solution function, 1st, 2nd, kth derivative defined at deriv_r.
function ueval = ANNEvalUserSolution(left, right, nth, deriv_l, deriv_r, ueval)
    
    len = right - left;
    if len < 1e-15, exit; end

%     %---------------------------------------------------------------
%     % Example 1: u(r,theta) = - cos(theta)
%     for i = 1 : nth
%         theta = (i-1) * 2*pi / nth;
%         sth = - cos(theta);
%         if deriv_l == 0   %   Evaluation of original function
%             ueval(i, 1) = sth;
%         elseif deriv_l == 1   %   Evaluation of derivative
%             ueval(i, 1) = 0;
%         elseif deriv_l == 2   %   Evaluation of 2nd derivative
%             ueval(i, 1) = 0;
%         else
%         end
%     end
% 
%     for i = 1 : nth
%         theta = (i-1) * 2*pi / nth;
%         sth = - cos(theta);
%         if deriv_r == 0   %   Evaluation of original function
%             ueval(i, 2) = sth;
%         elseif deriv_r == 1   %   Evaluation of derivative
%             ueval(i, 2) = 0;
%         elseif deriv_r == 2   %   Evaluation of 2nd derivative
%             ueval(i, 2) = 0;
%         else
%         end
%     end
    
%     %---------------------------------------------------------------
%     % Example 2: u(r,theta) = - r cos(\theta)
%     for i = 1 : nth
%         theta = (i-1) * 2*pi / nth;
%         sth = cos(theta);
%         if deriv_l == 0   %   Evaluation of original function
%             ueval(i, 1) = - left * sth;
%         elseif deriv_l == 1   %   Evaluation of derivative
%             ueval(i, 1) = - sth;
%         elseif deriv_l == 2   %   Evaluation of 2nd derivative
%             ueval(i, 1) = 0;
%         else
%         end
%     end
% 
%     for i = 1 : nth
%         theta = (i-1) * 2*pi / nth;
%         sth = cos(theta);
%         if deriv_r == 0   %   Evaluation of original function
%             ueval(i, 2) = - right * sth;
%         elseif deriv_r == 1   %   Evaluation of derivative
%             ueval(i, 2) = - sth;
%         elseif deriv_r == 2   %   Evaluation of 2nd derivative
%             ueval(i, 2) = 0;
%         else
%         end
%     end


    %---------------------------------------------------------------
    % Example 3: 
    %     S = 2 * pi * (xvec(i, 1) - left) / len - pi;
    %     zMat_sol(:, i) = sin(S) * cos(yvec);
    for i = 1 : nth
        theta = (i-1) * 2*pi / nth;
        sth = cos(theta);
        if deriv_l == 0   %   Evaluation of original function
            ueval(i, 1) = 0;
        elseif deriv_l == 1   %   Evaluation of derivative
            ueval(i, 1) = - 2*pi*sth/len;
        elseif deriv_l == 2   %   Evaluation of 2nd derivative
            ueval(i, 1) = 0;
        else
        end
    end

    for i = 1 : nth
        theta = (i-1) * 2*pi / nth;
        sth = cos(theta);
        if deriv_r == 0   %   Evaluation of original function
            ueval(i, 2) = 0;
        elseif deriv_r == 1   %   Evaluation of derivative
            ueval(i, 2) = - 2*pi*sth/len;
        elseif deriv_r == 2   %   Evaluation of 2nd derivative
            ueval(i, 2) = 0;
        else
        end
    end

return