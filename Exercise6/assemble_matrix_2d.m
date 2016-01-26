function [mat rhs] = assemble_matrix_2d(use_bspline, ndof, nurb, ...
    element_matrix, rhs_function)

% preamble and quadrature part of the code
...
    addpath('../Exercise1/');
addpath('../Exercise2/');
addpath('../Exercise3/');
% don't worry about this for the moment
assert(isa(element_matrix, 'function_handle'));
assert(isa(use_bspline, 'logical'));
assert(isa(rhs_function,'function_handle'));
%number of shape functions per element
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
        
        % print nurb_area for testing purposes
        nurb_area;
        
        
        
        % evaluate nurb basis functions and their derivatives
        %         for l = 1:nquad
        coeffs = zeros(nshp, nurb.number(1), nurb.number(2));
        p = nurb.order(1) - 1;
        q = nurb.order(2) - 1;
        i0 = bspline_findspan(nurb.number(1),p,xi1,nurb.knots{1});
        j0 = bspline_findspan(nurb.number(2),q,eta1,nurb.knots{2});
        connectivity = zeros(ndof,nshp);
        %         coeffs((j0-1)*(p+1)+i0,:,:) = 1;
        %         coeffs(:,:,:) = 1;
        %         coeffs(6,1,1) = 1;
        %         coeffs(5,2,1) = 1;
        %         coeffs(4,1,2) = 1;
        %         coeffs(3,2,2) = 1;
        %         coeffs(2,1,3) = 1;
        %         coeffs(1,2,3) = 1;
        l = 1;
        for j = j0-q:j0-q+q
            for i = i0-p:i0-p+p
                coeffs(l,i,j) = 1;
                l = l + 1;
            end
        end
        %         coeffs(1,2,3) = 1;
        %         coeffs(2,1,3) = 1;
        %         coeffs(3,2,2) = 1;
        %         coeffs(4,1,2) = 1;
        %         coeffs(5,2,1) = 1;
        %         coeffs(6,1,1) = 1;
        for k = 1:ndof
            for j = 1:nurb.order(2)
                for i = 1:nurb.order(1)
                    connectivity(k,(j-1)*(p+1)+i) = ...
                        ndof*((j0-q+j-2)*nurb.number(1)+i0-p+i-2)+k;
                    
                end
            end
        end
        
        % assemble coeffs and connectivity
        
        
        %         YOUR_CODE
        %
        if (use_bspline)
            % use B-Splines basis functions
            S = bspline_eval(nurb,coeffs, nshp, u2D);
            % derivatives in the reference domain
            dS_ref = bspline_derv_eval(nurb,coeffs,nshp,u2D);
        else
            % use NURBS basis functions
            
            % remember that nurb_(derv_)eval needs the coefficients
            % premultiplied with the weight
            % Remark: This step is only introduced to allow to test for
            % partition of unity of the NURBS basis functions. From an
            % algebraic point of view it does not matter whether we
            % premultiply the coefficients, it only introduces overhead in
            % the sense that the solution of the linear equation system
            % also needs to be mutiplied again with the weights before
            % plotting.
            for i=1:nshp
                coeffs(i,:,:) = coeffs(i,:,:) .* nurb.coeffs(4,:,:);
            end
            S = nurb_eval(nurb,coeffs, nshp, u2D);
            % derivatives in the reference domain
            dS_ref = nurb_derv_eval(nurb,coeffs,nshp,u2D);
        end
        dS_phys = zeros(nshp,2,nquad);
        
        %         % derivatives in the physical domain
%         for i=1:nquad
%             for j=1:nshp
%                 % compute dS_phys with the help of equation (3)
%                 %                 Use the '/' operator to compute the left inverse of the jacobian.
%                 K = dS_phys(j,:,i);
%                 dS_phys(j,:,i) = K + reshape(((inv(dF(:,:,i)))'*(dS_ref(j,:,i)')),1,2);
%             end
%         end

        for j = 1:nshp
            K = reshape(dS_ref(j,:,:),2,nquad);
            L = zeros(2,nquad);
            for i = 1:nquad
                L(:,i) = (inv(dF(:,:,i)))'*K(:,i);
            end
            dS_phys(j,:,:) = L;
        end
%
        % check partition of unity
       
        testSum = sum(S,1);
        testSumDervRef = sum(dS_ref,1);
        testSumDervPhys = sum(dS_phys,1);
        for i=1:nquad
            assert(abs(testSum(1,i) - 1.0) < 10^-10);
            for j=1:2
                assert(abs(testSumDervRef(1,j,i)) < 10^-10);
                assert(abs(testSumDervPhys(1,j,i)) < 10^-10);
            end
        end
%         dS_phys(:,:,:) = 0;
        
        ...
            
    for ishp = 1:nshp
        for jshp = 1:nshp
            matloc = element_matrix(ndof, quad_points, quad_weight,...
                jac_det,F, S(ishp,:), reshape(dS_phys(ishp,:,:),2,nquad), S(jshp,:), reshape(dS_phys(jshp,:,:),2,nquad));
             assert(isdiag(matloc) == 1);
            for l = 1:ndof
                mat(connectivity(l,ishp),connectivity(l,jshp)) = ...
                    mat(connectivity(l,ishp),connectivity(l,jshp)) + matloc(l,l);
            end
        end
        rhsloc = rhs_function(ndof, quad_points, quad_weight, jac_det, F,S(ishp,:), reshape(dS_phys(ishp,:,:),2,nquad));
        for l = 1:ndof
            rhs(connectivity(l,ishp)) = rhs(connectivity(l,ishp)) + rhsloc(l);
        end
    end
    
    
    end
end
end

