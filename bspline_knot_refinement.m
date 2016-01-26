function res = bspline_knot_refinement (bspline, coeffs, dimS, X, dir)

dim = numel(bspline.order);
assert(dim == numel(bspline.number));

if (dim == 1)
    assert(max(bspline.knots{1})>=max(X));
    assert(min(bspline.knots{1})<=min(X));

    r = size(X,2);

    res.order = bspline.order;
    res.number = bspline.number + r;
    res.knots{1} = sort([bspline.knots{1},X]);

    p = bspline.order -1;

    a = bspline_findspan(bspline.number, p, X(1), bspline.knots{1});
    b = bspline_findspan(bspline.number, p, X(r), bspline.knots{1});
    b=b+1;

    res.coeffs = zeros(dimS,res.number);

    res.coeffs(1:dimS,1:(a-p-1)) = coeffs(1:dimS,1:(a-p-1));
    res.coeffs(1:dimS,(b+r):res.number) = coeffs(1:dimS,b:bspline.number);

    i = b + p;
    k = b + p + r;
    for j=r:-1:1
        while X(j) <= bspline.knots{1}(i) && i > a
            res.coeffs(1:dimS,k-p-1) = coeffs(1:dimS,i-p-1);
            k = k-1;
            i=i-1;
        end
        res.coeffs(1:dimS,k-p-1) = res.coeffs(1:dimS,k-p);
        for l=1:p
            ind=k-p+l;
            alpha = res.knots{1}(k+l) - X(j);
            if(abs(alpha) == 0.0)
                res.coeffs(1:dimS,ind-1) = res.coeffs(1:dimS,ind);
            else
                alpha = alpha/(res.knots{1}(k+l)-bspline.knots{1}(i-p+l));
                res.coeffs(1:dimS,ind-1) = res.coeffs(1:dimS,ind-1).*alpha ...
                    + (1-alpha)*res.coeffs(1:dimS,ind);
            end
        end
        k = k-1;
    end
else
    error('Not supported yet');
end

end