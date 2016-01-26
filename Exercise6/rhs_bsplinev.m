function rhsloc = rhs_bsplinev(nurb, coeffs, ndof, quad_points, quad_weight, ...
    jac_det, F, v, dV)

rhsloc = zeros(ndof, 1);
nquad = size(quad_points,2);

S = bspline_eval(nurb,coeffs,ndof,quad_points);

for idof=1:ndof
    for j=1:nquad
        rhsloc(idof) = rhsloc(idof) + ...
            quad_weight(j) * jac_det(j) * v(j) * S(idof,j);
    end
end

end
