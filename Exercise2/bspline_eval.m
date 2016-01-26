function S = bspline_eval(bspline, coeffs, dimS, u)

% bspline:              holds the geometrical B-Spline structure
% coeffs(dimS,n1[,n2]): the coefficients we should use for evaluation, which are NOT
%                       necessarily the control points
% dimS:                 dimension of the space the curve/surface/DOF is embedded in,
%                       which implies that the first index of coeffs should have at
%                       least that many entries.
% u(dim,size_u):        the evaluation points with u(1,i) = $\xi_i$ and u(2,i) = $\eta_i$

dim = numel(bspline.order);
% Sanity checks of the bspline structure
assert(dim == numel(bspline.number));
assert(size(u,1) == dim);
assert(dim == ndims(coeffs)-1);
assert(size(coeffs,1) >= dimS);

size_u = size(u,2);

%Preallocate S
if (dim == 1)
    S = zeros(dimS,size_u);
    % degree of bspline
    p = bspline.order-1;
    % number of control points
    n = bspline.number;
    U = bspline.knots{1};
%     sz = size(coeffs)
    for j=1:size_u
        i = bspline_findspan(n,p,u(1,j),U);
%         u(1,j)
        N = bspline_basisfuns(i,u(j),p,U);
        for d=1:dimS
            for k=1:p+1
                S(d,j) = S(d,j) + N(k)*coeffs(d,i-p+k-1);            
            end
        end
    end
elseif (dim == 2)
    S = zeros(dimS,size_u);
    p  = bspline.order(1)-1;
    q  = bspline.order(2)-1;
    n1 = bspline.number(1);
    n2 = bspline.number(2);
    U1 = bspline.knots{1};
    U2 = bspline.knots{2};

    for j=1:size_u
        i1 = bspline_findspan(n1,p,u(1,j),U1);
        i2 = bspline_findspan(n2,q,u(2,j),U2);
        N = bspline_basisfuns(i1,u(1,j),p,U1);
        M = bspline_basisfuns(i2,u(2,j),q,U2);
        L = N*M';
        for d=1:dimS
            for k1=1:p+1
                for k2=1:q+1
                    S(d,j) = S(d,j) + L(k1,k2)*coeffs(d,i1-p+k1-1,i2-q+k2-1);            
                end
            end
        end
    end
end
end