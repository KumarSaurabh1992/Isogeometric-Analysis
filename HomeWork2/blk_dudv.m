function matloc = blk_dudv(ndof, quad_points, quad_weight, jac_det, ...
    F, U, dU, V, dV)

matloc = zeros(ndof, ndof);
nquad = size(quad_weight,2);

for idof=1:ndof
    for j=1:nquad
        matloc(idof,idof) = matloc(idof,idof) + ...
            quad_weight(j) * abs(jac_det(j)) * dot(dU(:,j), dV(:,j));
    end
end

end
