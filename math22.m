function math22
global  Pr  S  M  Rd K lamda lamda1 alpha A1 A2 A3 A4  A5 A6  phi1 phi2
Pr=20;
S=0;
lamda = 0;    
lamda1=0;
K=0;
Rd=0;
M=0;
alpha=pi/2;
ks1=40;                %Al2O3 = s1
ks2=400;               %Cu   =s2
rhos1=3970;
rhos2=8933;
phi1=0;  %0.1;
phi2=0;  %0.1;
rhocps1 = 3037050;
rhocps2 = 3439205;
rhof=1114;
kf=0.2520;
rhocpf=2690867;
betaTs1=85000;
betaTs2=167000;
betaTf=650000;
BB1=rhos1.*betaTs1;
BB2=rhof.*betaTf;
BB3=rhos2.*betaTs2;
FF1=BB1/BB2;
FF2=BB3/BB2;

sigmmaf=5.5.*10^-6;      %EG
sigmmas1=59.6.*10^6 ;
sigmmas2=35.83.*10^6;
A1=(1-phi2).*(1-phi1+phi1.*(rhos1/rhof))+phi2.*(rhos2/rhof);
A2=(1/(((1-phi1)^2.5).*((1-phi2)^2.5)));
A3=(1-phi2).*(1-phi1+phi1.*(rhocps1/rhocpf))+phi2.*(rhocps2/rhocpf);
A6=(1-phi1).*(1-phi2)+phi1.*FF1+phi2.*FF2;

EE1=5.*kf+ks1-5.*(kf-ks1).*phi1;      
EE2=5.*kf+ks1+(kf-ks1).*phi1;               
CC1=(EE1/EE2).*(kf);                              
EE3=5.*CC1+ks2-5.*(CC1-ks2).*phi2;                  
EE4=5.*CC1+ks2+(CC1-ks2).*phi2;             
A4=EE3/EE4;
DD1=2.*sigmmaf+sigmmas1-2.*(sigmmaf-sigmmas1).*phi1;       
DD2=2.*sigmmaf+sigmmas1+(sigmmaf-sigmmas1).*phi1;               
HH1=(DD1/DD2).*(sigmmaf);                              
DD3=2.*HH1+sigmmas2-2.*(HH1-sigmmas2).*phi2;                  
DD4=2.*HH1+sigmmas2+(HH1-sigmmas2).*phi2;             
A5=DD3/DD4;
% Initial values
for lamda1 = [0.1    0.5     0.9     1.3     1.7 ] 
sol1 = bvpinit(linspace(0,5,100),zeros(1,5));

%  bvp4c takes two functions defined below and give solution in structure
%  for2
  options = bvpset('RelTol',1e-3);
  sol = bvp4c(@bvpex1, @bcex1, sol1,options);
  x = sol.x;
  y = sol.y;
  if lamda1==0.1
   plot(x,y(3,:) ,'r--','linewidth',1.5)
  elseif lamda1==0.5
     plot(x,y(3,:)  ,'b--','linewidth',1.5)
  elseif lamda1==0.9
     plot(x,y(3,:)  ,'k--','linewidth',1.5)
   elseif lamda1==1.3
    plot(x,y(3,:)   ,'g--','linewidth',1.5)
    elseif lamda1==1.7
   plot(x,y(3,:)   ,'y--','linewidth',1.5)
%elseif d==0
% plot(x,y(3,:)   ,'m--','linewidth',1.5)
 else 
      plot(x,y(3,:)   ,'m--','linewidth',1.5)
   
   end
   
   
 xlabel('\eta')
  ylabel('theta(\eta)')
     hold on
     grid on
 %legend('d = 0.2','d= 0.5','d= 1','d=1.5','d=1.8');
%  hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%,%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%+(3.98)*(10.53)*((1-n1)*(1-n2)) + 
    
end
%  Here I defined first order ODEs
function ysol = bvpex1(~,y)
   ysol = [y(2);
          y(3);
          -(1/((A2/A1)-(lamda1/A1).*y(2)^2)).*(-y(2)^2-y(1).*y(3)-2*(lamda1/A1).*y(1).*y(2).*y(3)-(A6/A1).*lamda.*y(4).*cos(alpha)-(A5/A1).*M.*y(2)-A2.*K.*y(2));
          y(5);
          (Pr.*y(1).*y(5))/((A4/A3)+3.*A3.*Rd)];
end
% Here I define residual of boundary conditions
function res = bcex1(y, yinf)
res = [y(1)-S;
       y(2)-1;
       y(4)-1;  
       yinf(2);
       yinf(4)]; 
end
fprintf('%10.6f\t %10.6f\t %10.6f\n',y(3,1),y(5,1));
 
 %fprintf('%0.6f \n',y(3,1));
 % fprintf('%0.6f \n',y(5,1));
 % fprintf('%0.6f \n',y(7,1));
end          