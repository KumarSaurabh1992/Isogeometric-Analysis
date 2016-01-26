function nurb = generate_helix(torusradius, ncurl, refine)

nurb.number = [ 9 ];
nurb.order = [ 3 ];
p = nurb.order - 1;
nurb.knots{1} = [ 0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
nurb.coeffs = zeros(4,nurb.number);
s2 = 1.0/sqrt(2);
nurb.coeffs(1,1:nurb.number) = [-1 -s2 0 s2 1  s2  0 -s2 -1];
nurb.coeffs(3,1:nurb.number) = [ 0  s2 1 s2 0 -s2 -1 -s2  0];
nurb.coeffs(4,1:nurb.number) = [ 1  s2 1 s2 1  s2  1  s2  1];
% nurb.coeffs(4,:) = ones(2*nurb.number,1);
if (refine > 0)
  nurb = nurb_knot_refinement(nurb, refine);
end
nurb.coeffs_new = zeros(4,ncurl*(nurb.number - 1) + 1);
nurb.coeffs_new(1,1:nurb.number) = nurb.coeffs(1,1:nurb.number);
nurb.coeffs_new(3,1:nurb.number) = nurb.coeffs(3,1:nurb.number);
nurb.coeffs_new(4,1:nurb.number) = nurb.coeffs(4,1:nurb.number);

nurb.coeffs = nurb.coeffs_new;
% Remove premultiplied weights from control point coordinates
nurb.coeffs(1,1:nurb.number) = nurb.coeffs(1,1:nurb.number) ./ nurb.coeffs(4,1:nurb.number);
nurb.coeffs(3,1:nurb.number) = nurb.coeffs(3,1:nurb.number) ./ nurb.coeffs(4,1:nurb.number);
delta_Y = 1/(ncurl*(nurb.number - 1));
nurb.coeffs(2,:) = 1:-delta_Y:0;
% TODO: Add code to bend helix into circular shape.
U = nurb.knots{1};
U_tilda = U(p:nurb.number+2);
U_tilda = U_tilda./ncurl;
% Add premultiplied weights to new control point coordinates
x = nurb.coeffs(1,2:nurb.number);
z = nurb.coeffs(3,2:nurb.number);
w = nurb.coeffs(4,2:nurb.number);
x_tilda = x;
z_tilda = z;
w_tilda = w;
U = U_tilda;
for i = 2:ncurl
    x = [x,x_tilda];
    z = [z,z_tilda];
    w = [w,w_tilda];
    U = [U,(U_tilda(3:size(U_tilda,2)) +(i - 1)/ncurl)];
end
nurb.coeffs(1,:) = [nurb.coeffs(1,1),x];
nurb.coeffs(3,:) = [nurb.coeffs(3,1),z];
nurb.coeffs(4,:) = [nurb.coeffs(4,1),w];
nurb.knots{1} = [0,U,1];
r = torusradius + nurb.coeffs(1,:);
theta = 2*pi*nurb.coeffs(2,:);
nurb.coeffs(1,:) = r.*cos(theta);
nurb.coeffs(2,:) = r.*sin(theta);
nurb.coeffs(1,:) = nurb.coeffs(1,:) .* nurb.coeffs(4,:);
nurb.coeffs(2,:) = nurb.coeffs(2,:) .* nurb.coeffs(4,:);
nurb.coeffs(3,:) = nurb.coeffs(3,:) .* nurb.coeffs(4,:);
nurb.number = [(nurb.number - 1)*ncurl + 1];
end

