function i = bspline_findspan(n,p,xi,U)
% This routine searches for a given xi the corresponding knot span.
% n:        no. of control points
% p:        degree of the B-Spline/NURBS
% xi:       the xi we search the knot span for
% U(n+p+1): knot vector

	% Check for xi = 1.
% 	YOUR_CODE
    if(xi == 1)
        i = n;
    else
        for j = 1:size(U,2) - 1
            if ((xi >= U(1,j)) && (xi < U(1,j + 1)))
                i = j;
            end
        end
    end

	% Here you have to implement a linear (or even binary?)
	% search for the knotspan i.
% 	YOUR_CODE
end