function heat_eq_example()

addpath('../Exercise1/');
addpath('../Exercise2/');
addpath('../Exercise3/');

% we want to simulate on 25 elements
refinement_level=5;
ref_nurb = nurb_knot_refinement(generate_testnurb(),refinement_level-1);
ndof =1;
%Mass and Stiffness Matrix
delta_t = 1;
lambda = 0.01;

[mat, rhs] = assemble_matrix_2d(false, ndof, ref_nurb, @blk_dudv, @rhs_testfun);
[mat1, rhs] = assemble_matrix_2d(false, ndof, ref_nurb, @blk_uv, @rhs_testfun);
mat = mat*(delta_t*lambda);
mat = mat + mat1;
n = ref_nurb.number(1);
m = ref_nurb.number(2);

all_dofs = [1 : n*m];

% equation numbers that correspond to weighting functions with support on
% the boundary with eta=0
boundary1 = 1:n;
% boundary with eta=1
 boundary2 = (m-1)*n+1:1:n*m;
% boundary = union(boundary,1:n:m*(n-1));
% % boundary with eta=1
% boundary = union(boundary,n:n:n*m);
% inner_dofs = setdiff(all_dofs,boundary);

% initialize solution vector to zero, including the boundary terms
solution = ones(n*m,1)*30;
% solution(1:n) = -10; %eta = 0;
% solution((m-1)*n:1:n*m) = 100; %eta = 1
% % solve the linear system only for the inner degrees
% solution(inner_dofs) = mat(inner_dofs,inner_dofs) \ rhs(inner_dofs);
% 
% % reshape the solution into a coefficients array
% coeffs = reshape(solution,ndof,n,m);
% % remember that the nurb_eval routine needs the coefficients premultiplied
% % with the weights
% solution = ones(n*m,1)*30;

for i = 1:size(boundary1,2)
    mat(boundary1(1,i),:) = 0;
    mat(boundary1(1,i),boundary1(1,i)) = 1;
    mat1(boundary1(1,i),:) = 0;
    mat1(boundary1(1,i),boundary1(1,i)) = 1;
    solution(boundary1(1,i)) = -10;
end
for i = 1:size(boundary2,2)
    mat(boundary2(1,i),:) = 0;
    mat(boundary2(1,i),boundary2(1,i)) = 1;
    mat1(boundary2(1,i),:) = 0;
    mat1(boundary2(1,i),boundary2(1,i)) = 1;
    solution(boundary2(1,i)) = 100;
end


for time = 1:delta_t:15
    
    old_solution = solution;
%     solution = zeros(n*m,1);
%     solution(1:n) = -10; %eta = 0;
%     solution((m-1)*n:1:n*m) = 100; %eta = 1
%     rhs = old_solution - mat*solution;
    
    solution = mat\(mat1*old_solution);
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


end