function integral = b4_eta(x, y, x_data,y_data,H_data,p)
    integral =  p.b4_prefactor*p.dk * sum(p.K3_vec .* H_data.* ...
                besselj(0, p.K_vec .* sqrt((x - x_data').^2 + (y - y_data').^2 )));
    integral=sum(integral);
end