function Tan_22
global  Pr M K R B Ec Cr Nt Nb Sc A1 A2 A3 A4  phi1 phi2 F0
global phi3 A5 A r  lambda  S Fr Q
Pr=7;
Sc=2;
Nb=0.5;
Nt=0.1;
W=0.3;
A=0.1;
eta=2;
r=1;
lambda=0.3;
n=0.1;
S=0.1;
Fr=0.3;
M=0.2;
K=0.1;
Q=0.01;
R=1.5;
B=1.2;
Ec=0.2;
Cr=0.1;
phi1=0.05;          %  Single - Al2O3/EG
phi2=0.05;             %  Hy  -    Al2O3+Cu/EG
phi3=0.05;             %  Tri -    Al2O3+Cu+TiO2/EG
rhos1=3970;
rhos2=8933;
rhos3=4250;
ks1=40;
ks2=400;
ks3=8.9538;
rhocps1=3037050;
rhocps2=3439205;
rhocps3=2916350;
rhof=1115;
kf=0.253;
rhocpf=2709450;
sigmmaf=0.107;
sigmmas1=3.5*10^7;
sigmmas2=5.96*10^7;
sigmmas3=1.0*10^-12;
A1=(1/(((1-phi1)^2.5).*((1-phi2)^2.5).*((1-phi3)^2.5)));
A2=(1-phi1).*((1-phi2).*((1-phi3)+phi3.*(rhos3/rhof))+phi2.*(rhos2/rhof))+phi1.*(rhos1/rhof);
A5=(1-(phi1+phi2+phi3))+phi1.*(rhocps1/rhocpf)+phi2.*(rhocps2/rhocpf)+phi3.*(rhocps3/rhocpf);
D1=ks3+2.*kf-2.*phi3.*(kf-ks3);
D2=ks3+2.*kf+phi3.*(kf-ks3);
T1=(D1/D2).*kf;
D3=ks2+2.*T1-2.*phi2.*(T1-ks2);
D4=ks2+2.*T1+phi2.*(T1-ks2);
T2=(D3/D4).*T1;
D5=ks1+2.*T2-2.*phi1.*(T2-ks1);
D6=ks1+2.*T2+phi1.*(T2-ks1);
T3=(D5/D6).*T2;
A4=T3/kf;
C1=sigmmas3+2.*sigmmaf-2.*phi3.*(sigmmaf-sigmmas3);
C2=sigmmas3+2.*sigmmaf+phi3.*(sigmmaf-sigmmas3);
F1=(C1/C2).*sigmmaf;
C3=sigmmas2+2.*F1-2.*phi2.*(F1-sigmmas2);
C4=sigmmas2+2.*F1+phi2.*(F1-sigmmas2);
F2=(C3/C4).*F1;
C5=sigmmas1+2.*F2-2.*phi1.*(F2-sigmmas1);
C6=sigmmas1+2.*F1+phi1.*(F2-sigmmas1);
F3=(C5/C6).*F2;
A3=F3/sigmmaf;
%P1=A1/A2;
%P2=A3/A2;
% Initial values
for Q = [0      0.05      0.1       0.14      0.18] 
solinit = bvpinit(linspace(0, 10, 101),[0 0 0 0 0 0 0]); % Initial guess

%  bvp4c takes two functions defined below and give solution in structure
%  for2
  options = bvpset('RelTol',1e-7);
  sol = bvp5c(@odefun, @bcfun, solinit,options); % Solve BVP
  x = sol.x;
  y = sol.y; % y values
  if Q==0
   plot(x,y(4,:)    ,'r', 'linewidth',1.5)
  elseif Q==0.05
     plot(x,y(4,:)  ,'b','linewidth',1.5)
  elseif Q==0.1
     plot(x,y(4,:)  ,'k','linewidth',1.5)
   elseif Q==0.14
    plot(x,y(4,:)   ,'m','linewidth',1.5)
    elseif Q==0.18
   plot(x,y(4,:)   ,'g','linewidth',1.5)
%elseif F0==0
 %plot(x,y(2,:)   ,'m','linewidth',1.5)
 %elseif F0==0
% plot(x,y(2,:)   ,'p','linewidth',1.5)
 
 else 
      plot(x,y(4,:)   ,'m','linewidth',1.5)
   
   end
   
   
 xlabel('\eta')
  ylabel('\theta(\eta)')
     hold on
     grid on
 legend('Pr = 7','Pr = 9','Pr = 12','Pr = 15 ','Pr = 18');
%  hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%,%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+(3.98)*(10.53)*((1-n1)*(1-n2)) + 
    
end
%  Here I defined first order ODEs
function dydx = odefun(~,y)
dydx = [y(2);
        y(3);
        ((A+(A3/A2)*M+(A1/A2)*K)*y(2)+(1+Fr)*(y(2))^2-(y(1)-(A*eta)/2)*y(3))/((1-n)+n*W*y(3));
        y(5);
        (-((R/(A5*Pr))*(3*(1+(B-1)*y(4))^2*(B-1)*(y(5))^2)+(y(1)-(A*eta)/2)*y(5)+(1/A5)*Q*y(4)+(A1/A5)*n/2*W*Ec*(y(3))^3+(A1/A5)*(1-n)*Ec*(y(3))^2+Ec*((A1/A5)*K+(A3/A5)*M)*(y(2))^2+Fr*Ec*(y(2))^3+(1/A5)*(Nb*y(7)*y(5)+Nt*(y(5))^2)))/((A4/(A5*Pr))*(1+lambda*y(4))+(R/(A5*Pr))*(1+(B-1)*y(4))^3);
        y(7);
        ((Sc*A*eta)/2)*y(7)+Sc*Cr*y(6)-Sc*y(1)*y(7)-(Nt/Nb)*((-((R/(A5*Pr))*(3*(1+(B-1)*y(4))^2*(B-1)*(y(5))^2)+(y(1)-(A*eta)/2)*y(5)+(1/A5)*Q*y(4)+(A1/A5)*n/2*W*Ec*(y(3))^3+(A1/A5)*(1-n)*Ec*(y(3))^2+Ec*((A1/A5)*K+(A3/A5)*M)*(y(2))^2+Fr*Ec*(y(2))^3+(1/A5)*(Nb*y(7)*y(5)+Nt*(y(5))^2)))/((A4/(A5*Pr))*(1+lambda*y(4))+(R/(A5*Pr))*(1+(B-1)*y(4))^3))];
    end
% Here I define residual of boundary conditions
function res = bcfun(y, yinf)
res = [y(1)-S-F0;
       y(2)-r;
       y(4)-1;
       y(6)-1;  
       yinf(2);
       yinf(4);
       yinf(6);
       ]; 
end
%fprintf('%10.6f\t %10.6f\t %10.6f\n',y(end,2),y(end,4), y(end,6));
%fprintf('\nFirst solution:\n');
fprintf('f"(0) = %7.9f\n',             y(3));   %reduced skin friction
fprintf('-theta^prime(0)= %7.9f\n',   -y(5));  %reduced local Nusselt number
fprintf('-Phi^prime(0)=%7.9f\n',      -y(7));  %reduced local Sherwood number
%fprintf('Cfx=%7.9f\n',4*A2*eta^0.5*(1+(1/beta))*(y(3)));   %skin friction
%fprintf('Nux=%7.9f\n',-2*A4*eta^0.5*y(5));  %local Nusselt number
%fprintf('Shx=%7.9f\n',-2*eta^0.5*y(7));  %local Sherwood number
%fprintf('\n');
%fprintf('%3.2f   %7.6f    %7.6f\n',(4*eta^0.5)/(P1)*(1+(1/beta))*(y(3))),-2*eta^0.5*P4*y(5)),-2*eta^0.5*y(7))
 %fprintf('%0.6f \n',y(3,1));
 % fprintf('%0.6f \n',y(5,1));
 % fprintf('%0.6f \n',y(7,1));
end           