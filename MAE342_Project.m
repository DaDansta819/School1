% Daniel Andrews & Darrin Neiman
clc, clear
format compact 
r2=2.25;r3=1.5; r4=1.875; r1=3.64; r04b=3.375; 
w2=1.0; a2=.75;
T2=linspace(155.78,329.0,232);
T3=linspace(303.87,29.44,232);
T4=linspace(349.79,264.09,232);
 for p=1:232
     t2=T2(p); t4=T4(p); t3=T3(p); 
    %Velocity
    a=[-r3.*sind(t3), r4.*sind(t4); r3.*cosd(t3), -r4.*cosd(t4)];
    b=[r2.*sind(t2).*w2; -r2*cosd(t2).*w2];
    c=a\b;
    w3=c(1);
    w4=c(2);
    %Acceleration
    d=[-r3.*sind(t3), r4.*sind(t4); r3.*cosd(t3), -r4.*cosd(t4)];
    e=[r2.*sind(t2).*a2+cosd(t2)*w2^2+r3*cosd(t3).*w3.^2-r4.*cosd(t4).*w4.^2 ;
        -r2.*(cosd(t2)*a2-sind(t2).*w2.^2)+r3.*sin(t3).*w3.^2-r4*sind(t4).*w4.^2];
    f=d\e;
    a3=f(1);
    a4=f(2);
        % Coupler Anaylsis
        % Rpo2=RBo2+RPB=RBo2P
        % Position
        Xp(p)=r04b*cosd(t4);
        Yp(p)=r04b*sind(t4);
        %Velocity
        Xpd(p)=-r04b*sind(t4)*w4;
        Ypd(p)=r04b*cosd(t4)*w4;
        %Acceleration
        Xpdd(p)=-r04b*(sind(t4)*a4+cosd(t4*w4^2));
        Ypdd(p)=r04b*(cosd(t4)*a4+sind(t4*w4^2));
         W3(p)=w3; W4(p)=w4; A3(p)=a3; A4(p)=a4;
 end
 W3(1)
 W4(1)
 A3(1)
 A4(1)
fprintf('Omega3 at full deployment is %4.6f rad/s \n',W3(232))
fprintf('Omega4 at full deployment is %4.6f rad/s \n\n',W4(232))
fprintf('Alpha3 at full deployment is %4.6f rad/s^2 \n',A3(232))
fprintf('Alpha4 at full deployment is %4.6f rad/s^2 \n\n',A4(232))
fprintf('Position of wheel in x direction at full deployment is %4.6f feet \n',Xp(232))
fprintf('Position of wheel in y direction at full deployment is %4.6f feet \n\n',Yp(232))
fprintf('Velocity of wheel in x direction at full deployment is %4.6f ft/s \n',Xpd(232))
fprintf('Velocity of wheel in y direction at full deployment is %4.6f ft/s \n\n',Ypd(232))
fprintf('Acceleration of wheel in x direction at full deployment is %4.6f ft/s^2 \n',Xpdd(232))
fprintf('Acceleration of wheel in y direction at full deployment is %4.6f ft/s^2 \n',Ypdd(232))
%Graphs
P=1:232;
X=subplot(2,1,1);plot(P,W3,P,W4); xlabel('steps');ylabel('rad/sec');
title('Omega shifts during deployment');legend('w3','w4')
Alpha=subplot(2,1,2);plot(P,A3,P,A4);xlabel('steps');ylabel('rad/sec^2');
title('Alpha shifts during deployment');legend('a3','a4')
figure
position=subplot(2,2,1);plot(Xp,Yp);title('Position');
xlabel('meters (x)'); ylabel('meters (y)');
Velocity=subplot(2,2,2);plot(Xpd,Ypd);title('Velocity');
xlabel('m/s (x)');ylabel('m/s (y)');
Acceleration=subplot(2,2,3);plot(Xpdd,Ypdd);title('Acceleration');
xlabel('m/s^2 (x)');ylabel('m/s^2 (y)');


%  T2=linspace(155.78,329.93,232);
%  T4=linspace(349.79,264.09,232);
%  T3=linspace(303.87,29.44,232);
%T=0; W2=1.5; A2=.75;
% A=zeros(1,232); %t3=t2; t4=t2; % B=t2; w3=t2; w4=t2
%tsq=linspace(124.15,29,232);
% for a=1:5
%     for b=1:a;
%         c=a*b;
%     end
%     c;
% end
% t2=zeros(232);
% for a=1:232
%     t2=linspace(158.16,390,a);
%     for t4=linspace(349.79,264.09,232)
%     D=t2*t4;
%      for t2=linspace(158.16,390,232)
%         D=C*A;
%     end
%     D
%     end
%     D;
% end
%  A=1:length(t2)
%      B=1:length(t4)
%          C=1:length(t3)
%t2=linspace(158.16,390,a);
%  tsq=linspace(124.15,29,232);
%     t2=p; t3=T+p; t4=T+129-p; B=T+10+p; w2=W2; a2=A2; %Made up have to change
%  for p=1:232
%      for q=1:232
%          for r=232
%                 p=t2; q=t3; r=t4;
%                 
%          end
%      end
%  end

%     end
%         end
%     end
% end

%Rpo4=RAo4+RPA
% Postion
% Xp=r4*cosd(t4)+rpa*cosd(t3+B);
% Yp=r2*sind(t2)+rpa*sind(t3+B);
% %Velocity
% Xpd=-r3*sind(t2)*w2-rpa*sind(t3+B)*w3;
% Ypd=r2*cosd(t3)*w2+rpa*cosd(t3+B)*w3;
% %Acceleration
% Xpdd=-r2*(sind(t2)*a2+cosd(t2*w2^2))-rpa*(sind(t3+B)*a3+cosd(t3+B)*w3^2);
% Ypdd=r2*(cosd(t2)*a2+sind(t2*w2^2))-rpa*(cosd(t3+B)*a3+sind(t3+B)*w3^2);
%W3(t2,t4,t3)=w3; W4(t2,t4,t3)=w4; A3(t2,t3,t4)=a3; A4(t2,t3,t4)=a4; 
% %         end
% %     end
% %    
% % end
% % %W3,W4,A3,A4,Xp,Yp,Xpd,Ypd, Xpdd,Ypdd
% % 
% % 
% % %Graphs
% % P=1:45;
% % X=subplot(2,1,1);plot(P,W3,P,W4); xlabel('degrees');ylabel('rad/sec');
% % Alpha=subplot(2,1,2);plot(P,A3,P,A4);xlabel('degrees');ylabel('rad/sec^2');
% % figure
% % position=subplot(2,2,1);plot(Xp,Yp);title('Position'); xlabel('meters'); ylabel('meters');
% % Velocity=subplot(2,2,2);plot(Xpd,Ypd);title('Velocity'); xlabel('m/s');ylabel('m/s');
% % ACC=subplot(2,2,3);plot(Xpdd,Ypdd);title('Acceleration');xlabel('m/s^2');ylabel('m/s^2');



%     t2=t2*(pi/180)
%     % 4 Bar 
% % Loop Closure r2+r3-r4-r1=0
% % Position
% syms r2 t2 r4 r1 r3 t3 t4
% eqn1 = r2*cos(t2)+r3*cos(t3)-r4*cos(t4)-r1 == 0;
% eqn2 = r2*sin(t2)+r3*sin(t3)-r4*sin(t4) == 0;
% sol=solve([eqn1,eqn2],[t3,t4]);
% T3Sol= sol.t3;
% T4Sol= sol.t4;
% t3
% t4

%    for u=1:1:2
%                 if u==1
%              t2=158.16; t3=-56; t4=349.79;
%                 else 
%                    t2=390; t3=29; t4=264.09; 
% 
% 
% %Acceleration
% %d=[-r3*sind(t3), r4*sind(t4); r3*cosd(t3), -r4*cosd(t4)];
% %e=[r2*sind(t2(p))*a2(p)+cosd(t2(p)).*w2(p).^2+r3*cosd(t3(p)).*w3(p).^2-r4*cosd(t4(p)).*w4(p).^2 ; -r2*(cosd(t2(p))*a2(p)-sind(t2(p)).*w2(p).^2)+r3*sin(t3(p)).*w3(p).^2-r4*sind(t4(p)).*w4(p).^2];
% %f=inv(d)*e;
% %a3(p)=f(1);
% %a4(p)=f(2);
% % 6-bar coupler
% % Rpo2=RAo2+RPA
% % Postion
% %Xp=r2*cosd(t2(p))+rpa*cosd(t3(p)+B(p))
% %Yp=r2*sind(t2(p))+rpa*sind(t3(p)+B(p))
% %Velocity
% %Xpd=-r3*sind(t2)*w2-rpa*sind(t3+B)*w3
% %Ypd=r2*cosd(t3)*w2+rpa*cosd(t3+B)*w3
% %Acceleration
% %Xpdd=-r2*(sind(t2)*a2+cosd(t2*w2^2))-rpa*(sind(t3+b)*a3+cosd(t3+b)*w3^2)
% %Ypdd=r2*(cosd(t2).*a2+sind(t2.*w2.^2))-rpa*(cosd(t3+b)*a3+sind(t3+b).*w3.^2)
% 
% %w3
% % Force Anaylsis
% 
% 
% 
% % Daniel Andrews
% % angles=[deployment, retraction, landing]
% %t2=[30,0,30]; t3=[30,0,30]; t4=[110,0,110]; B=[40,0,40]; 
% %t2=30; t3=30; t4=110; B=40; w2=1; a2=1; %t2=theta2 B=beta
% % 4 Bar 
% % Loop Closure r2+r3-r4-r1=0
% % Position
% %syms t3 t4
% %eqn1 = r2*cosd(t2)+r3*cosd(t3)-r4*cosd(t4)-r1 == 0;
% %eqn2 = r2*sind(t2)+r3*sind(t3)-r4*sind(t4) == 0;
% %sol=solve([eqn1,eqn2],[t3,t4]);
% %T3Sol= sol.t3
% %T4Sol= sol.t4
% %Velocity
% %r2*e^(i*t2)
% %-r2*w2*sind(t2)-r3*sind(t3)*w3+r4*w4*sind(t4)=0
% %r2*w2*cosd(t2)+r3*cosd(t3)*w3-r4*cosd(t4)*w4=0
% %Acceleration
% %-r3*sind(t3)*a3-r3*cosd(t3)*w3^2-r4*sind(t4)*a4+r4*cosd(t4)*w4^2=r2*sind(t2)*a2+r2*cosd(t2)*w2^2
% %r3*cos(t3)*a3-r3*sind(t3)*w3^2-r4*cosd(t4)*a4+r4*sind(t4)*w4^2=-r2*cosd(t2)*a2+r2*sind(t2)*w2^2
% Daniel Andrews


