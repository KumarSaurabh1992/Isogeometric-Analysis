function [] = draw_fig(nurb,points)
deltaX = 1/(points(1)-1);
%deltaY = 1/(points(2)-1);

[X] = meshgrid(0:deltaX:1,1);

u = zeros(1, points(1));
u(1,:) = reshape(X,1,[]);
S = nurb_eval(nurb,nurb.coeffs,3,u);
plot3(S(1,:),S(2,:),S(3,:),'LineWidth',2);
xlabel('x');
ylabel('y');
zlabel('z');
end