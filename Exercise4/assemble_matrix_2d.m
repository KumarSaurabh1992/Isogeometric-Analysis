function [mat rhs] = assemble_matrix_2d(use_bspline, ndof, nurb, ...
    element_matrix, rhs_function)
addpath('../Exercise1/');
addpath('../Exercise2/');
addpath('../Exercise3/');
% don't worry about this for the moment
%assert(isa(element_matrix, 'function_handle'));
%assert(isa(use_bspline, 'logical'));

% number of shape functions per element
nshp = nurb.order(1) * nurb.order(2);
alldof = ndof*nurb.number(1)*nurb.number(2);
mat = spalloc(alldof,alldof,alldof*nshp);
rhs = zeros(alldof,1);

%remove double entries in the knot arrays
xknots = unique(nurb.knots{1});
yknots = unique(nurb.knots{2});

[quad_weight_1d, quad_points_1d] = quad_rule();
nquad_1d = size(quad_weight_1d, 2);

nquad = nquad_1d * nquad_1d;
quad_weight = zeros(1,nquad);
quad_points = zeros(2,nquad);

nurb_area = 0.0;

% Loop over all elements, which are spaced from knot to knot
for ie = 1:size(xknots,2)-1
    for je = 1:size(yknots,2)-1
        % build quad rules
        xi1 = xknots(ie);
        eta1 = yknots(je);
        deltaXi = xknots(ie+1) - xi1;
        deltaEta = yknots(je+1) - eta1;
        % area of the element in the reference domain
        area = deltaXi * deltaEta;
        % Compute the weights according to equation (4).
        % It is easier if you use the function 'reshape', which makes
        % it possible to avoid the integer division and the modulo operation.
%         quad_weight = reshape(quad_weight_1d'*quad_weight_1d,nquad,1)*area/4;
        
        
        % Compute the quadrature points according to formula (5).
        % If possible use the functions 'repmat' and 'reshape'.
        k = 1;
        for l = 1:nquad_1d
            for m = 1:nquad_1d
                quad_weight(1,k) = quad_weight_1d(l)*quad_weight_1d(m)*area/4; 
                 quad_points(1,k) = 0.5*(deltaXi)*quad_points_1d(l) + 0.5*(xknots(ie+1) + xi1);
                 quad_points(2,k) = 0.5*(deltaEta)*quad_points_1d(m) + 0.5*(yknots(je+1) + eta1);
                 k = k + 1;
            end
        end
%         quad_points(1,:) = YOUR_CODE;
%         quad_points(2,:) = YOUR_CODE;
        
        % evaluate the geometry mapping and its derivatives
%         coeffs2D = ones(2,nurb.number(1),nurb.number(2)) .* nurb.coeffs(4,:,:);
        u2D =quad_points;
        F = nurb_eval(nurb, nurb.coeffs, 2, u2D);
        dF = nurb_derv_eval(nurb,nurb.coeffs,2,u2D);
        jac_det = zeros(1,nquad);
        for i=1:nquad
            % evaluate the absolute value of the Jacobian's determinant
            jac_det(i) = det(dF(:,:,i));
        end
        
        % Compute the NURBS' area
        nurb_area = nurb_area + sum(quad_weight.*jac_det);
    end
end

% print nurb_area for testing purposes
nurb_area

end