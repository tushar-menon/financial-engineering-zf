
%syms fgammas(gammas) ft(eta)
alpha=0.05;
nBanks=5;
thetaMin=0.08;
b=0.1;
s=ones(nBanks);
p=ones(nBanks)*(2-0.9/20);
x=1*ones(nBanks);
step=1/3000;
len=3000;

ft=@(eta) exp(-b*eta);
fgammas=@(gammas) exp(-b*gammas);
dft=@(eta) -b*ft(eta);
dfgammas=@(gammas) -b*fgammas(gammas);


q=[1,zeros(1,len-1)];
psis=zeros(nBanks,len);
gammas=zeros(nBanks,len);
theta0=ones(nBanks,1);


for i=1:nBanks
  theta0(i,:)=(x(i)+s(i)*q(1)-p(i))/(alpha*s(i)*q(1)); 
end
thetas=[theta0,zeros(nBanks,len-1)];
for i =2:len
    for banki=1:nBanks
    k=(i-1)/len;
    gammaDot=0;
    psiDot=0; 
    
    if thetas(banki,i-1)< thetaMin
        z=vpa(((1-alpha*thetaMin)*(s(banki)-gammas(banki,i-1)))/(alpha*thetaMin*ft(k)*fgammas(gammas(banki,i-1))));
        %gammaDot=vpa(-(z*dft(k))*fgammas(gammas(banki,i-1))/(1+z*ft(k)*dfgammas(gammas(banki,i-1))));
        gammaDot=-qDot*(s(banki) - gammas(banki, i-1))*(p(banki)-x(banki)-psis(banki, i-1))/(q(i-1)*((s(banki)-gammas(banki,i-1))*q(i-1) - (p(banki)-x(banki)-psis(banki, i-1))));
        psiDot=gammaDot*q(i-1);
       
    end
    psis(banki,i)=psis(banki,i-1)+psiDot*step;
    gammas(banki,i)=gammas(banki,i-1)+gammaDot*step;
    
    qDot=dft(k)*fgammas(sum(gammas(:,i-1)))+gammaDot*ft(k)*dfgammas(gammas(banki,i-1));
    q(i)=q(i-1)+qDot*step;
    thetaDot=(qDot*(s(banki)-gammas(banki,i-1))*(p(banki)-x(banki)-psis(banki,i-1))+gammaDot*q(i-1)*((s(banki)-gammas(banki,i-1))*q(i-1)-(p(banki)-x(banki)-psis(banki,i-1))))/(alpha*((s(banki)-gammas(i-1))*q(i-1))^2);
    thetas(banki,i)=thetas(banki,i-1)+thetaDot*step;
    
%     if thetas(banki,i)<0
%         keyboard;
%     end
    end
end

figure 
plot(gammas)
