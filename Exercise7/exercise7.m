function exercise7()
addpath('../Exercise1/');
addpath('../Exercise2/');
addpath('../Exercise3/');

% we want to simulate on 25 elements
refinement_level=5;
ref_nurb = nurb_knot_refinement(generate_testnurb(),refinement_level-1);

% degrees of freedom per shape function

ndof =1;
[mat, rhs] = assemble_matrix_2d(false, ndof, ref_nurb, @blk_dudv, @rhs_testfun);

n = ref_nurb.number(1);
m = ref_nurb.number(2);

all_dofs = [1 : n*m];

% equation numbers that correspond to weighting functions with support on
% the boundary with eta=0
boundary = 1:n:m*(n-1);
% boundary with eta=1
boundary = union(boundary,n:n:n*m);
% boundary with xi=0
boundary = union(boundary,1:n);
% boundary with xi=1
boundary = union(boundary,(m-1)*n+1:1:n*m);

% subtract the boundary set from all degrees
inner_dofs = setdiff(all_dofs,boundary);

% initialize solution vector to zero, including the boundary terms
solution = zeros(n*m,1);

% solve the linear system only for the inner degrees

solution(inner_dofs) = mat(inner_dofs,inner_dofs) \ rhs(inner_dofs);

% reshape the solution into a coefficients array
coeffs = reshape(solution,ndof,n,m);
% remember that the nurb_eval routine needs the coefficients premultiplied
% with the weights
for i=1:ndof
    coeffs(i,:,:) = coeffs(i,:,:).*ref_nurb.coeffs(4,:,:);
end

% show a surface plot of the solution
draw_nurb_surf(ref_nurb, [30 30], ...
    @(xi,eta) nurb_eval(ref_nurb, coeffs, ndof, [xi; eta]));

end

