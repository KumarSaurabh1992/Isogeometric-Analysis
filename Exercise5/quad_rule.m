function [quad_weight, quad_points] = quad_rule()

% fifth-order accurate Gaussian quadrature (c.f. Hughes page 142)

tmp1 = sqrt(3.0 ./ 5.0);
tmp2 = 5.0 ./ 9.0;
tmp3 = 8.0 ./ 9.0;
quad_weight = [tmp2 tmp3 tmp2];
quad_points = [-tmp1 0.0 tmp1];

end