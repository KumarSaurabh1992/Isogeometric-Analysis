function [] = test_nurb_eval()
addpath('../Exercise1');
% Testing 1D-implementation: partition of unity
nurb1D.coeffs(1,:) = [1.00000   0.707106781186548   0.00000];
nurb1D.coeffs(4,:) = [1.00000   0.707106781186548   1.00000];
nurb1D.number = [ 3 ];
nurb1D.order = [ 3 ];
nurb1D.knots{1} = [0 0 0 1 1 1];
% remember that coefficients are stored premultiplied
coeffs1D = ones(1,nurb1D.number(1)) .* nurb1D.coeffs(4,:);
u1D = [0.0, 0.25, 0.5, 0.75, 1.0];
S1D = nurb_eval(nurb1D, coeffs1D, 1, u1D);
assert(all(all(S1D == 1.0)));

% Testing 2D-implementation: partition of unity
nurb2D = generate_testnurb();
% remember that coefficients are stored premultiplied
coeffs2D = ones(1,nurb2D.number(1),nurb2D.number(2)) .* nurb2D.coeffs(4,:,:);
u2D = [0.0, 0.25, 0.5, 0.75, 1.0; 0.0, 0.25, 0.5, 0.75, 1.0];
S2D = nurb_eval(nurb2D, coeffs2D, 1, u2D);
assert (all(all(S2D == 1.0)));

end