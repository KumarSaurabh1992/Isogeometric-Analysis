function rhsloc = rhs_testfun(ndof, quad_points, quad_weight, ...
    jac_det, F, V, dV)

rhsloc = zeros(ndof, 1);
nquad = size(quad_points,2);

for idof=1:ndof
    for j=1:nquad
        rhsloc(idof) = rhsloc(idof) + ...
            quad_weight(j) * jac_det(j) * V(j) * val_rhs(F(1,j),F(2,j));
    end
end

end

function val = val_rhs(x,y)
  val = ((8 - 9*(sqrt(x^2+y^2)))/(x^2+y^2))*sin(2*atan(y/x));
%val = 0;
end
