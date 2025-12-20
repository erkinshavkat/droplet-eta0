function [rhs1, rhs2] = compute_rhs_full_IF(phi_hat, eta_hat, t, p)
% Computes the RHS of the hyperbolic part
% IF = integrating factor (combine with other IF files) -- M.Durey
% Including highest frequency in 2nd-derivative symbol (for surface
% tension term)
rhs1 = -p.g(t).*eta_hat + p.Bo*p.K2_deriv.*eta_hat;
rhs2 = DtN(phi_hat,p); % DtN computes phi_z_hat    
end