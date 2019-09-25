
%syms fGamma(gamma) ft(eta)
alpha=[0.05,0.03];
numAssets=2;
z=zeros(1,numAssets);

thetaMin=0.08;
b=[0.1,0.2];
s=[1,1.5];
p=3.15;
x=1;
step=1/3000;
len=3000;
% eta=1;
ft=@(eta,b) exp(-b*eta);
fGamma=@(gamma,b) exp(-b*gamma);
dft=@(eta,b) -b*ft(eta,b);
dfGamma=@(gamma,b) -b*fGamma(gamma,b);


q=ones(numAssets,len);
psi=zeros(numAssets,len);
gamma=zeros(numAssets,len);
theta0=1;
% (x+s.*q(:,1)-p)/(alpha.*s.*q(:,1));
theta=[theta0,zeros(1,len-1)];
eta=zeros(1,len);
for i =2:len
    k=(i-1)/len;
    GammaDot=zeros(numAssets,1);
    qDot=zeros(numAssets,1);
    psiDot=zeros(numAssets,1);
    for j=1:numAssets
        if theta(i-1)< thetaMin
    %             z(j)=vpa(((1-alpha(j)*thetaMin)*(s(j))-gamma(j,i-1)))/(alpha(j)*thetaMin*ft(k,b(j))*fGamma(gamma(j,i-1),b(j)));
    %             GammaDot(j)=vpa(-(z(j)*dft(k,b(j)))*fGamma(gamma(j,i-1))/(1+z(j)*ft(k,b(j))*dfGamma(gamma(j,i-1),b(j))));
                GammaDot(j)=etaDot*s(j);
                psiDot(j)=etaDot*s(j)*q(j,i-1);

        end      
        qDot(j)=dft(k,b(j))*fGamma(gamma(j,i-1),b(j))+GammaDot(j)*ft(k,b(j))*dfGamma(gamma(j,i-1),b(j))
        numerator=sum((qDot(:).*s(:)-eta(i-1)*s(:).*qDot(:)).*(alpha(:).*(s(:)-eta(i-1).*s(:)).*qDot(:)));
        denom=(x+sum(psiDot(:)+(s(:)-eta(i-1)*s(:)).*q(:,i-1)-p)).*sum(s(:).*q(:,i-1));
%       etaDot=sum((qDot(:).*s(:)+eta(i-1)*s(:).*qDot(:)).*(alpha(:).*(s:(:)-eta(i-1).*s(:)).*qDot(:)))./(x+sum(psiDot(:))+(s(:)-eta(i-1)*s(:)).*q(:,i)-p).*sum(s.*q(:,i))+sum(-alpha(:).*(qDot(:,i).*s(:)+eta(i-1).*s.*qDot(:))./(s(:).*q(:,i)));
        etaDot=numerator./denom +sum(-alpha(:).*(qDot(:).*s(:)-eta(i-1)*qDot(:).*s(:))./(s(:).*q(:,i-1)));
        psi(j,i)=psi(j,i-1)+psiDot(j)*step;
        gamma(j,i)=gamma(j,i-1)+GammaDot(j)*step;
        q(j,i)=q(j,i-1)+qDot(j)*step;
        eta(i)=eta(i-1)+etaDot*step;
    end

    
%     thetaDot=(qDot)*(s-gamma(i-1))*(p-x-psi(i-1))+GammaDot*q(i-1)*((s-gamma(i-1))*q(i-1)-(p-x-psi(i-1))))/(alpha*((s-gamma(i-1))*q(i-1))^2);
%     theta(i)=theta(i-1)+thetaDot*step;
    theta(i)=(x+sum(psi(:,i))+sum((s(:)-gamma(:,i)).*q(:,i))-p)./sum((alpha(:).*((s(:)-gamma(:,i)).*q(:,i))));
    
    if theta(i)<0
        keyboard;
    end
end

figure 
plot(gamma)
