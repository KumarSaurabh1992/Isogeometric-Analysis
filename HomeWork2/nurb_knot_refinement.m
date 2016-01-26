function new_nurb = nurb_knot_refinement (nurb, level)

% This function inserts in every non-degenerate knot-interval knots. The
% count of inserted knots per interval is given by level.

% The general problem is to search P_i' and w_i' such that
%   sum N_i w_i P_i / sum N_i w_i = sum N_i' w_i' P_i'/ sum N_i' w_i'
% where N_i are the bspline basisfunctions of the original NURB and N_i'
% the bspline basisfunctions of the new one. The problem can easily be
% split into two equation systems:
%   sum N_i' w_i' = sum N_i w_i
%   sum N_i' w_i' P_i' = sum N_i w_i P_i
% which are independent since we store w_i' P_i' premultiplied in the nurb
% structure. Now one could transform this into a linear equation system by
% evaluating the equations at sufficient u's to ensure that the matrix has
% full rank. In contrast to this here we persue a least-square approch,
% where we consider the weak formulation of these equations
%   int N_j'(u) sum N_i'(u) w_i' = int N_j'(u) sum N_i(u) w_i
% which is fed to the isogeometric solver using a dummy nurb as integration
% domain.

dim = numel(nurb.order);
assert(dim == numel(nurb.number));

if (dim == 1)
    xknots = unique(nurb.knots{1});
    newknots = [];
    for i=1:size(xknots,2)-1
        deltaX = xknots(i+1)-xknots(i);
        newknots = [newknots ([1:level].*deltaX ./ (level+1) + xknots(i)) ];
    end
    new_nurb = bspline_knot_refinement(nurb,nurb.coeffs,4,newknots,1);
elseif (dim == 2)
    new_nurb.order = nurb.order;
    
    %remove double entries in the knot arrays
    xknots = unique(nurb.knots{1});
    yknots = unique(nurb.knots{2});
    new_nurb.knots{1} = nurb.knots{1};
    for i=1:size(xknots,2)-1
        deltaX = xknots(i+1)-xknots(i);
        new_nurb.knots{1} = [new_nurb.knots{1} ...
            ([1:level].*deltaX ./ (level+1) + xknots(i)) ];
    end
    new_nurb.knots{2} = nurb.knots{2};
    for i=1:size(yknots,2)-1
        deltaY = yknots(i+1)-yknots(i);
        new_nurb.knots{2} = [new_nurb.knots{2} ...
            ([1:level].*deltaY ./ (level+1) + yknots(i)) ];
    end
    new_nurb.knots{1} = sort(new_nurb.knots{1});
    new_nurb.knots{2} = sort(new_nurb.knots{2});
    
    new_nurb.number(1) = nurb.number(1) + ...
        level*(size(xknots,2)-1);
    new_nurb.number(2) = nurb.number(2) + ...
        level*(size(yknots,2)-1);
    
    % dummy values, might result in an overhead to generate trash numbers,
    % but it is an easy way to realize the knot refinement
    new_nurb.coeffs = zeros(4,new_nurb.number(1),new_nurb.number(2));
    for i=1:new_nurb.number(1)
        for j=1:new_nurb.number(2)
            new_nurb.coeffs(1,i,j) = i-1;
            new_nurb.coeffs(2,i,j) = j-1;
            new_nurb.coeffs(4,i,j) = 1.0;
        end
    end
    
    % now determine the coefficients by using the isogeometric solver to
    % solve a least-square problem that corresponds to the weak formulation
    rhs_w = @(ndof, quad_points, quad_weight, jac_det, F, v, dV) ...
        rhs_bsplinev(nurb,nurb.coeffs(4,:,:), ndof, quad_points, ...
        quad_weight, jac_det, F, v, dV);
    
    [matW, rhsW] = assemble_matrix_2d(true,1,new_nurb,@blk_uv,rhs_w);
    if isempty(ver('Octave'))
        % Use a CG method to solve the system in combination with a incomplete
        % Cholesky preconditioner
        cholW = ichol(matW);
        w = reshape(pcg(matW,rhsW,1e-8,20,cholW,cholW'), 1, ...
            new_nurb.number(1), new_nurb.number(2));
    else
        % Octave does not provide ichol
        w = reshape(matW\rhsW, 1, ...
            new_nurb.number(1), new_nurb.number(2));
    end
    
    % We can assume that the NURBS surface is embedded in R^3 since the
    % the forth entry is reserved for the weights and therefore enough memory
    % is allocated.
    dimS=3;

    rhs_coeff = @(ndof, quad_points, quad_weight, jac_det, F, v, dV) ...
        rhs_bsplinev(nurb,nurb.coeffs, ndof, quad_points, ...
        quad_weight, jac_det, F, v, dV);
    
    [matCoeff, rhsCoeff] = assemble_matrix_2d(true,dimS,new_nurb,@blk_uv, ...
        rhs_coeff);
    
    if isempty(ver('Octave'))
        cholCoeff = ichol(matCoeff);
        new_nurb.coeffs(1:dimS,:,:) = reshape(pcg(matCoeff,rhsCoeff,1e-8,20, ...
            cholCoeff,cholCoeff'), dimS, new_nurb.number(1), new_nurb.number(2));
    else
        new_nurb.coeffs(1:dimS,:,:) = reshape(matCoeff\rhsCoeff, dimS, ...
            new_nurb.number(1), new_nurb.number(2));
    end
    new_nurb.coeffs(4,:,:) = w;
end
draw_nurb_surf(new_nurb,[50 50],@(xi,eta) xi^2);


end
