function integral = b4_eta_trunc(x, y, x_data,y_data,H_data,k_llim,k_ulim,p)
    n = size(H_data,2);
    integral = 0;
    trunc_indices = p.K_vec >= k_llim & p.K_vec <= k_ulim;
    for impact = 1:n

        integral = integral +  p.b4_prefactor*p.dk * sum(p.K3_vec(trunc_indices) .* H_data(trunc_indices,impact) .* ...
                    besselj(0, p.K_vec(trunc_indices) .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
        
    end
end