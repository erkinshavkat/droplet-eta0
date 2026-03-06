function H_corr = H_correction(t,k,p)
            a= sqrt(p.d0*p.Bo);b = sqrt(p.d0/p.Bo)*p.G/2;
            alpha = a.*k.^2+b;
            alpha_recip = 1/b .* exp(-a.*k.^2/b);
            H_corr=exp(-2*p.nu0*k.^2.*t).*sin(alpha.*t).*alpha_recip./k.^2 *b/a;