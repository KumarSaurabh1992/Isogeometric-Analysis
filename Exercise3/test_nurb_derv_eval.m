function test_nurb_derv_eval()
addpath('../Exercise1');
addpath('../Exercise2');
nurb = generate_testnurb();

% 5000 evaluation points should be sufficient
size_u = 5000;

% use unifrnd to generate u
a = 0;
b = 1;
du = 1/(size_u - 1);

% eval the derivatives at u
u(1,:) = a:du:b;
u(2,:) = a:du:b;
dS = nurb_derv_eval(nurb,nurb.coeffs,2,u);
area = 0.0;
% calculate the area according to equation (9)
for i = 1:size_u
    area = area+det(dS(:,:,i));
end

area = area/size_u
% with 5000 points the error is in most cases smaller than 0.02
assert(abs(area-pi*3/4) < 0.02);

end