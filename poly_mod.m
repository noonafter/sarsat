function [out_mod] = poly_mod(in_pre,mod_ploy)
% divide A(x) by g(x),output residue polynomial r(x), 
% that is, A(x) = q(x)g(x) + r(x), or, r(x) = A(x) mod g(x)
% in_pre: A(x), in_pre is longer than mod_poly
% mod_poly: g(x), the first and last bit of mod_poly is 1
% out_mod: the vector of r(x),
n = length(in_pre);
n_k_one = length(mod_ploy); % n-k+1
if n < n_k_one
    out_mod = [zeros(1,n_k_one-1-n) in_pre];
    return;
end
k = n - n_k_one + 1;
for idx_shift = 0:k-1
    if (mod_ploy(1) * in_pre(1+idx_shift) == 1) % requrie mod_ploy(1) is nonzero
        in_pre(1+idx_shift:n_k_one+idx_shift) = xor(in_pre(1+idx_shift:n_k_one+idx_shift), mod_ploy);
    end
end
out_mod = in_pre(end-n_k_one+2:end);
end

