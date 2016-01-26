function [ nurb_area ] = integrate_hw01( nurb_circle,nurb_helix )

%remove double entries in the knot arrays
xknots = unique(nurb_circle.knots{1});
yknots = unique(nurb_helix.knots{1});

[quad_weight_1d, quad_points_1d] = quad_rule();
nquad_1d = size(quad_weight_1d, 2);

nquad = nquad_1d * nquad_1d;
quad_weight = zeros(1,nquad);
quad_points = zeros(2,nquad);

nurb_area = 0.0;
% Loop over all elements, which are spaced from knot to knot
for ie = 1:size(xknots,2)-1
    xi1 = xknots(ie);
    del_Xi = xknots(ie+1) - xi1;
    for k = 1:nquad
        l = floor((k-1)/nquad_1d) + 1;
        quad_points(1,k) = 0.5*(del_Xi) *quad_points_1d(l) + 0.5*(xknots(ie+1) + xi1);
    end
    u_circle = quad_points(1,:);
     F_circle = nurb_eval(nurb_circle, nurb_circle.coeffs, 3, u_circle);
        dF_circle = nurb_derv_eval(nurb_circle,nurb_circle.coeffs,...
            3,u_circle);
    for je = 1:size(yknots,2)-1
        % build quad rules
        eta1 = yknots(je);
        del_Eta = yknots(je+1) - eta1;
        area = del_Xi * del_Eta;
        for k = 1:nquad
            l = floor((k-1)/nquad_1d) + 1;
            m = mod(k-1,nquad_1d) + 1;
            quad_weight(1,k) = quad_weight_1d(l) * quad_weight_1d(m) * area/4;
            quad_points(2,k) = 0.5*(del_Eta)*quad_points_1d(m) + ...
                0.5*(yknots(je+1) + eta1);
        end
        
        
        u_helix  = quad_points(2,:);
        
       
        F_helix = nurb_eval(nurb_helix, nurb_helix.coeffs, 3, u_helix);
        dF_helix = nurb_derv_eval(nurb_helix,nurb_helix.coeffs,3,...
            u_helix);
        
        cross_prod = cross(dF_circle,dF_helix);
        difference = F_helix-F_circle;
        val = zeros(1,nquad);
        for i=1:nquad
            val(1,i)=  dot(cross_prod(:,1,i), difference(:,i)) / ((norm(difference(:,i)))^3 );
        end
        nurb_area = nurb_area + sum(quad_weight.*val);
    end
    
end

nurb_area = (nurb_area)/(4*pi);
end

