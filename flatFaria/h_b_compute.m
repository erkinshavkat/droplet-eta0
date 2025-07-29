N=100;

test_h = linspace(0.0001,0.006,N);
test_b= zeros(N,1);
for i=1:length(test_h)
    p=setup_IF_matt(0,test_h(i),zeros(2),2,1,1);
    test_b(i)=tanh(p.kf_mean * test_h(i))/p.kf_mean;
end
bofh = @(h)interp1(test_h,test_b,h,'cubic');
hofb = @(b)interp1(test_b,test_h,b,'cubic');
h1=0.005; h2=0.001;
b1=bofh(h1);
b2=bofh(h2);
meanb=b1*4/5 + b2*1/5;
meanh=hofb(meanb)
h0=0.0004;b0=bofh(h0);
meanchb=b0 *2/3 + meanb*1/3;
meanchh=hofb(meanchb)