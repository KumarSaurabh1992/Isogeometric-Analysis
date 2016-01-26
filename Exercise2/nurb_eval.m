function S = nurb_eval(nurb, coeffs, dimS, u)

size_u = size(u,2);
% Preallocate S
S = zeros(dimS,size_u);
dim = numel(nurb.order);
if (dim == 1)
    % Perform two bspline_eval calls to determine the numerator and
    % the denominator of Equation (4).
    num = bspline_eval(nurb,coeffs,dimS,u);
    den = bspline_eval(nurb,nurb.coeffs(4,:),1,u);
    for i = 1:size(num,2)
        S(:,i) = num(:,i)./den(i);
    end
    
elseif (dim == 2)
    % Almost the same code as in the dim == 1 case.
     num = bspline_eval(nurb,coeffs,dimS,u);
    den = bspline_eval(nurb,nurb.coeffs(4,:,:),1,u);
    for i = 1:size(num,2)
        S(:,i) = num(:,i)/den(i);
    end
    
end
end