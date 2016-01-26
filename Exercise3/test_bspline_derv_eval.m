function [] = test_bspline_derv_eval()

% Testing 1D-implementation: partition of unity
bspline1D.coeffs(1,:) = [1.00000   0.707106781186548   0.00000];
bspline1D.number = [ 3 ];
bspline1D.order = [ 3 ];
bspline1D.knots{1} = [0 0 0 1 1 1];
coeffs1D = ones(1,bspline1D.number(1));
u1D = [0.0, 0.25, 0.5, 0.75, 1.0];
dS1D = bspline_derv_eval(bspline1D, coeffs1D, 1, u1D);
assert(all(all(all(dS1D == 0.0))));

% Testing 2D-implementation: partition of unity
bspline2D = generate_testnurb();
coeffs2D (1,:,:) = ones(bspline2D.number(1),bspline2D.number(2));
u2D = [0.0, 0.25, 0.5, 0.75, 1.0; 0.0, 0.25, 0.5, 0.75, 1.0];
dS2D = bspline_derv_eval(bspline2D, coeffs2D, 1, u2D);
assert (all(all(all(dS2D == 0.0))));

end