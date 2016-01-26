function nurb = generate_circle(radius, refine)

nurb.number = [ 9 ];
nurb.order = [ 3 ];
nurb.knots{1} = [ 0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
nurb.coeffs = zeros(4,nurb.number);
s2 = 1.0/sqrt(2);
nurb.coeffs(1,:) = [-1 -s2  0  s2 1  s2  0 -s2 -1] * radius;
nurb.coeffs(2,:) = [ 0 -s2 -1 -s2 0  s2  1  s2  0] * radius;
nurb.coeffs(4,:) = [ 1  s2  1  s2 1  s2  1  s2  1];

if (refine > 0)
  nurb = nurb_knot_refinement(nurb,refine);
end
end
