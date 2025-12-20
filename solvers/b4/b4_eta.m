function integral = b4_eta(x, y, x_data,y_data,H_data,p)
    %% Formerly B4, now 3.10, computes eta at (x,y)
    % x, y (1x1), position to valuate eta
    % xdata, y_data (Nimpacts,1), position of past impacts
    % H_data (Nk, Nimpacts), H associated with each past impact
    integral =  p.b4_prefactor*p.dk * sum(p.K3_vec .* H_data.* ...
                besselj(0, p.K_vec .* sqrt((x - x_data').^2 + (y - y_data').^2 )));
    integral=sum(integral);
end