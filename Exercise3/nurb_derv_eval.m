function dS = nurb_derv_eval(nurb, coeffs, dimS, u)

dim = numel(nurb.order);
assert(dim == numel(nurb.number));
assert(size(u,1) == dim);
assert(dim == ndims(coeffs)-1);
assert(size(coeffs,1) >= dimS);

size_u = size(u,2);

% Preallocate dS
dS = zeros(dimS,dim,size_u);

% Compute A,W and the corresponding jacobians dA,dW
A = bspline_eval(nurb,coeffs,dimS,u);
W = bspline_eval(nurb,nurb.coeffs(4,:,:),1,u);

dA = bspline_derv_eval(nurb,coeffs,dimS,u);
dW = bspline_derv_eval(nurb,nurb.coeffs(4,:,:),1,u);


for k=1:size_u
    % Compute dS(:,:,k)
    
    dS(:,:,k) = (dA(:,:,k) - (A(:,k)*dW(:,:,k))/W(1,k))/W(1,k);
end
end