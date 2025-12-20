function integral = eta_varDtN(x, y, x_data,y_data,H_data,DtN,p)
    integral =  p.b4_prefactor/p.d0.*p.dk * sum(p.K_vec .*DtN.* H_data.* ...
                besselj(0, p.K_vec .* sqrt((x - x_data').^2 + (y - y_data').^2 )));
    integral=sum(integral);
end