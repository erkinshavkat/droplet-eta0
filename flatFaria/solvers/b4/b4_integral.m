function integral = b4_integral(x, y, impact, H, p, k_correction)
    if nargin < 6
        k_correction = 1;
    end
    K_vec_scaled = p.K_vec * k_correction;
    K3_vec_scaled = K_vec_scaled.^3;
    x_data = p.x_data;
    y_data = p.y_data;
    integral = p.b4_prefactor * sum(K3_vec_scaled .* H(:,impact) .* ...
        besselj(0, K_vec_scaled .* sqrt((x - x_data(impact)).^2 + (y - y_data(impact)).^2 )));
end