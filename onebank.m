
%syms fGamma(gamma) ft(eta)
alpha=0.05;
thetaMin=0.08;
b=0.1;
s=1;
p=2-0.1/20;
x=1;
step=1/3000;
len=3000;

ft=@(eta) exp(-b*eta);
fGamma=@(gamma) exp(-b*gamma);
dft=@(eta) -b*ft(eta);
dfGamma=@(gamma) -b*fGamma(gamma);


q=[1,zeros(1,len-1)];
psi=[0,zeros(1,len-1)];
gamma=[0,zeros(1,len-1)];
theta0=(x+s*q(1)-p)/(alpha*s*q(1));
theta=[theta0,zeros(1,len-1)];



for i =2:len
    k=(i-1)/len;
    GammaDot=0;
    psiDot=0;
    if theta(i-1)<thetaMin
        z=vpa(((1-alpha*thetaMin)*(s-gamma(i-1)))/(alpha*thetaMin*ft(k)*fGamma(gamma(i-1))));
        GammaDot=vpa(-(z*dft(k))*fGamma(gamma(i-1))/(1+z*ft(k)*dfGamma(gamma(i-1))));
        psiDot=GammaDot*q(i-1); 
    end
    psi(i)=psi(i-1)+psiDot*step;
    gamma(i)=gamma(i-1)+GammaDot*step;
    qDot=dft(k)*fGamma(gamma(i-1))+GammaDot*ft(k)*dfGamma(gamma(i-1));
    q(i)=q(i-1)+qDot*step;
    num=(qDot*(s-gamma(i-1))*(p-x-psi(i-1))+GammaDot*q(i-1)*((s-gamma(i-1))*q(i-1)-(p-x-psi(i-1))));
    denom=(alpha*((s-gamma(i-1))*q(i-1))^2);
    thetaDot=num/denom;
    theta(i)=theta(i-1)+thetaDot*step;
    if theta(i)<0
        keyboard;
    end
end

figure 
plot(theta)
xlabel('Time')
ylabel('\theta')
