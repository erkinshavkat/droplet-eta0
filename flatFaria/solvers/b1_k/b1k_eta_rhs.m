function [d_eta, d_etaprime] = b1k_eta_rhs(eta_hat, etaprime_hat, t, p)

    
    d_eta=etaprime_hat;
    d_etaprime = -(p.d0*p.Bo* p.K2_deriv.*p.K2_deriv + p.d0*p.g(t)*p.K2_deriv).*eta_hat;
end