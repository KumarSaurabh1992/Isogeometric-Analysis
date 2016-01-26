function matloc = blk_uv(ndof, quad_points, quad_weight, jac_det, ...
    F, U, dU, V, dV)

% Did you forget to use 'squeeze' to remove the singleton dimensions of the
% physical derivatives array in assemble_matrix_2d?
assert(ndims(dU)==2)
assert(ndims(dV)==2)

matloc = zeros(ndof, ndof);
nquad = size(quad_weight,2);

for idof=1:ndof
    for j=1:nquad
        matloc(idof,idof) = matloc(idof,idof) + ...
            quad_weight(j) * jac_det(j) * U(j) * V(j);
    end
end

end
