% %Daniel Andrews
%Force Anaylsis
clc, clear
disp('What the forces are without considering Fi and Moment of Inertia')
names={'F1x','F1y','F2x','F2y','F3x','F3y','F4x','F4y','T'};
a=zeros(9,9);
r1=3.0625; r2=2.25; r3=1.5; r4op=3.375; r4=1.875; %Lengths
r2m=8.09; r3m=3.74; r4m=5.46; % Masses
a(1)=-1; a(2,2)=-1; a(4,3)=-1; a(5,4)=-1; a(7,5)=-1; a(8,6)=-1; 
a(1,3)=1;a(2,4)=1; a(4,5)=1;a(5,6)=1; a(7,7)=1; a(8,8)=1;
a(9,9)=1;
tf=68.05;fgx=125.16; fgy=894.88;fa=990;
for p=1:3
       if p==1 %Deployment 
       disp('Deployment')
       t2=263.9; t3=29.44; t4=329.93; fgy=0;fgx=0; 
        elseif p==2 %retraction
           disp('Retraction')
           w3=-0.41591; w2=-0.94692; a3=2.3066; a4=-2.4067;
           t2=349.85; t3=303.87; t4=155.78; fa=-fa;
        elseif p==3 %landing
           disp('Landing')
          t2=263.9; t3=29.44; t4=329.93;   fgy=894.88; fgx=125.16; 
       end 
   
   a(3,3)=-r2*sind(t2); a(3,4)=r2*cosd(t2); a(6,5)=-r3*sind(t3);a(6,6)=r2*cosd(t3);a(9,7)=-r4*sind(t4);a(9,8)=r4*cosd(t4); %These values calculated
   
  A=array2table(a,'VariableNames',names,'RowNames',names)
   b=zeros(9,1);
   b(3)=fgy*r4op*cosd(t2)-r4op*sind(t2)*fgx;
   b(7)=-fa*cosd(tf);
   b(8)=fa*sind(tf);
   B=b
   c=a\b;
   C=array2table(c,'RowNames',names)
 
end

disp('What the forces are with considering Fi and Moment of Inertia')
% %Daniel Andrews
%Force Anaylsis
 clear
names={'F1x','F1y','F2x','F2y','F3x','F3y','F4x','F4y','T'};
a=zeros(9,9);
r1=3.0625; r4=2.25; r3=1.5; r2op=3.375; r2=1.875; %Lengths
%Center of gravity given that it is in the middle of the bar.
cg=[ r4/2, r3/2, r2op/2, r2/2]; 
m=[5.46, 3.74, 8.09, 8.09*(r2./r2op)];% Masses lb
u=cg.*m;
w=1/6; %width ft
I=(1/12)*w^2; % Inertia
MI=[r2,r3,r2op].^2.*I; % Inertia 
a(1)=-1; a(2,2)=-1; a(4,3)=-1; a(5,4)=-1; a(7,5)=-1; a(8,6)=-1; 
a(1,3)=1;a(2,4)=1; a(4,5)=1;a(5,6)=1; a(7,7)=1; a(8,8)=1;
a(9,9)=1; % Torque value
tf=68.05;fgx=125.16; fgy=894.88;fa=990;
w4=1; a4=.75;
M4=MI(3)*-a4;    
for p=1:3
      if p==1 %Deployment 
        disp('Deployment')
      t4=329.93; t2=263.9; t3=29.44;  fgy=0;fgx=0; 
       w3=-1.665534; w2=1.279751; a3=-2.413351; a2=-1.152065;
           Ax2=-r2*(sind(t2)*a2+cosd(t2*w2^2));
           Ay2=r2*(cosd(t2)*a2+sind(t2*w2^2));
           Ax3=-r3*(sind(t3)*a3+cosd(t3*w3^2));
           Ay3=r3*(cosd(t3)*a3+sind(t3*w3^2));
           Ax4=-r4*(sind(t4)*a4+cosd(t4*w4^2));
           Ay4=r4*(cosd(t4)*a4+sind(t4*w4^2));
           M3=MI(2)*-a3; M2=MI(1)*-a2;
   elseif  p==2 %retraction
           disp('Retraction')
           w3=-0.41591; w2=-0.94692; a3=2.3066; a2=-2.4067;
           t2=349.85; t3=303.87; t4=155.78; fa=-fa;
           Ax2=-r2*(sind(t2)*a2+cosd(t2*w2^2));
           Ay2=r2*(cosd(t2)*a2+sind(t2*w2^2));
           Ax3=-r3*(sind(t3)*a3+cosd(t3*w3^2));
           Ay3=r3*(cosd(t3)*a3+sind(t3*w3^2));
           Ax4=-r4*(sind(t4)*a4+cosd(t4*w4^2));
           Ay4=r4*(cosd(t4)*a4+sind(t4*w4^2));
           M3=MI(2)*-a3; M2=MI(1)*-a2;
   elseif p==3 %landing
           disp('Landing')
           % Acceleration of Mechanism is 0
          t2=263.9; t3=29.44; t4=329.93; fgy=894.88; fgx=125.16;
          Ax2=0;Ay2=0;Ax3=0;Ay3=0;Ax4=0;Ay4=0;
       end 
   
   r2x=r2op*cosd(t2);
   r2y=r2op*sind(t2);
   fix2=u(1)*Ax2;fiy2=u(1)*Ay2;
   fix3=u(2)*Ax3;fiy3=u(2)*Ay3;
   fix4=u(3)*Ax3;
   fiy4=u(3)*Ay3;
   a(3,3)=-r2*sind(t2); a(3,4)=r2*cosd(t2); a(6,5)=-r3*sind(t3);a(6,6)=r2*cosd(t3);a(9,7)=-r2*sind(t4);a(9,8)=r2*cosd(t4); %These values calculated
  AA=array2table(a,'VariableNames',names,'RowNames',names)
   b=zeros(9,1);
   b(1)=-fgx-fix2;
   b(2)=-fgy-fiy2;
   b(3)=fgy*r2x-r2y*fgx+fiy2*.5*r2x-fix2*.5*r2y-M2;
   b(4)=-fix3;
   b(5)=-fiy3;
   b(6)=-M3+fiy3*.5*r3*cosd(t3)-fix3*.5*r3*sind(t3);
   b(7)=-fa*cosd(tf);
   b(8)=fa*sind(tf);
   b(9)=-M4+fa*sin(tf)*r4*cos(t4)-fa*cos(tf)*r4*sind(t4)+fix4*.5*r4*sind(t4)+fiy4*.5*r4*cosd(t4);
 BB=b
   answers=a\b;
   CC=array2table(answers,'RowNames',names)
end

    
%     Xpdd(p)=-r*(sind(t)*a+cosd(t*w^2));
%         Ypdd(p)=r04b*(cosd(t)*a4+sind(t*w4^2));


 %a(3,3)=40; a(3,4)=45; a(6,5)=4;a(6,6)=4;a(9,5)=4; a(9,6)=4; %These values calculated
%a(3,3)=3; a(3,4)=3; a(6,5)=3;a(6,6)=3;a(9,5)=3; a(9,6)=3; %These values calculated

  
