function N = bspline_basisfuns(i,xi,p,U)
% This routine determines for a given xi and its corresponding knot span i
% the non-vanishing basis functions N(1:p+1), which are evaluated at xi.
% i:        to xi corresponding knot span
% xi:       the xi we search the knot span for
% p:        degree of the B-Spline/NURBS
% U(n+p+1): knot vector

	% Preallocate the N array
% % 	YOUR_CODE
    U = U';
    N = zeros(p + 1,1); 
%     i = bspline_findspan(4,2,xi,U);
	% Unroll the first iteration of the outer loop, since
	% this is the special branch of the Cox-deBoor formula.
% 	YOUR_CODE
    N(1) = 1;
    

	for j=2:p+1
		% For k=1 there is no dependence on N(k-1) of the previous run.
		saved = 0.0;
		for k=1:j-1
			% Compute $N_{i-j+k,j-1}$ according to the Cox-deBoor formula
           
            num2 = U(i + k ) - xi;
            den2 = U(i + k ) - U(i -j + k + 1);
%             if (den1 < 1e-6)
%                 A = 0;
%             else
%                 A = (num1/den1)*saved;
%             end
            if (den2 < 1e-6)
                B = 0;
            else
                B = (num2/den2)*N(k);
            end
            tmp = saved + B;
% 			tmp = YOUR_CODE;
			% Precompute the first summand of $N_{i-j+k+1,j-1}$ for the next
			% iteration of the inner loop.
            num1 = (xi - U(i - j + k + 1 ));
            den1 = (U(i + k ) - U(i - j + k + 1));
            if (den1 < 1e-6)
                A = 0;
            else
                A = (num1/den1)*N(k);
            end
			saved = A;
			N(k) = tmp;
		end
		% Unroll the last iteration of the inner loop
% 		YOUR_CODE
        N(k + 1) = saved;
	end
end