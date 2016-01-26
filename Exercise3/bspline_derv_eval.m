function dS = bspline_derv_eval(bspline, coeffs, dimS, u)
% This routine calculates for every u(:,i) the corresponding jacobian of the
% B-Spline function and stores it in dS(:,:,i).

% bspline:              holds the geometrical B-Spline structure
% coeffs(dimS,n1[,n2]): the coefficients we should use for evaluation, which are not
%                       necessarily the control points
% dimS:                 dimension of the space the curve/surface/DOF is embedded in,
%                       which implies that the first index of coeffs should have at
%                       least that many entries.
% u(dim,size_u):        the evaluation points with u(1,i) and u(2,i) being the $\xi$-
%                       and $\eta$-direction of the i-th evaluation point

dim = numel(bspline.order);
% Sanity checks
assert(dim == numel(bspline.number));
assert(size(u,1) == dim);
assert(dim == ndims(coeffs)-1);
assert(size(coeffs,1) >= dimS);

size_u = size(u,2);

% Preallocate dS

dS = zeros(dimS,dim,size_u);
if (dim == 1)
    % Fill a data structure dbspline with knots= $\tilde{U}$,number= $\tilde{n}$,..., as well as
    % an array dcoeffs with $\tilde{\bm{P}}_i$. Then evaluate this dbspline using bspline_eval.
    
    p = bspline.order - 1;
    n = bspline.number - 1;
    U = bspline.knots{1};
    U_tilda = U(2:size(U,2) - 1);
    dcoeffs = zeros(dimS,n);
    for d = 1:dimS
        for k = 1:size(dcoeffs,2)
            dcoeffs(d,k) = (p/(U(k+p+1) - U(k+1)))*(coeffs(d,k+1) - coeffs(d,k));
        end
    end
    dbspline = bspline;
    dbspline.number = n;
    dbspline.order = p;
    dbspline.knots{1} = U_tilda;
    dS(:,:,:) = bspline_eval(dbspline,dcoeffs,dimS,u);
    
elseif (dim == 2)
    % Compute derivative in $\xi$-direction:
    % Fill the dbspline structure as in the dim == 1 case, leaving the $\eta$-related
    % information unchanged. Store the result of bspline_eval in the first column
    % of dS.
    %     YOUR_CODE
    
    p = bspline.order(1) - 1;
    q = bspline.order(2) - 1;
    n1 = bspline.number(1);
    n2 = bspline.number(2);
    
    U = bspline.knots{1};
    U_tilda = U(2:size(U,2) - 1);
    dcoeffs = zeros(dimS,n1 - 1,n2);
    for d = 1:dimS
        for k1 = 1:(n1-1)
            
                dcoeffs(d,k1,:) = (p/(U(k1+p+1) - U(k1+1)))*(coeffs(d,k1+1,:) - coeffs(d,k1,:));
            
        end
    end
%     dbspline = bspline;
    dbspline.number = [n1-1,n2];
    dbspline.order = [p,bspline.order(2)];
    dbspline.knots{1} = U_tilda;
    dbspline.knots{2} = bspline.knots{2};
    dS(:,1,:) = bspline_eval(dbspline,dcoeffs,dimS,u);
    
    U = bspline.knots{2};
    U_tilda = U(2:size(U,2) - 1);
    dcoeffs = zeros(dimS,n1,n2 - 1);
    for d = 1:dimS
        
            for k2 = 1:n2-1
                dcoeffs(d,:,k2) = (q/(U(k2+q+1) - U(k2+1)))*(coeffs(d,:,k2+1) - coeffs(d,:,k2));
            end
        
    end
   % dbspline = bspline;
    dbspline.number = [n1 n2-1];
    dbspline.order = [bspline.order(1) q];
    dbspline.knots{2} = U_tilda;
    dbspline.knots{1} = bspline.knots{1};
    dS(:,2,:) = bspline_eval(dbspline,dcoeffs,dimS,u);
    
    
    % Compute derivative in $\eta$-direction:
    
end
end