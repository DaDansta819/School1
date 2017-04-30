
%Force Anaylsis
 clc,clear
names={'F1x','F1y','F2x','F2y','F3x','F3y','F4x','F4y','T'};
a=zeros(9,9);
r1=3.0625; r2=2.25; r3=1.5; r4op=3.375; r4=1.875; %Lengths
%Center of gravity given that it is in the middle of the bar.
cg=[ 24/2, r3/2, r4op/2, r4/2]; 
m=[5.46, 3.74, 8.09, 8.09*(r4./r4op)];% Masses lb
u=cg.*m;
w=1/6; %width ft
I=(1/12)*w^2; % Inertia
MI=[r2,r3,r4op].^2.*I; % Inertia 
M3=MI(2); M4=MI(3);
a(1)=-1; a(2,2)=-1; a(4,3)=-1; a(5,4)=-1; a(7,5)=-1; a(8,6)=-1; 
a(1,3)=1;a(2,4)=1; a(4,5)=1;a(5,6)=1; a(7,7)=1; a(8,8)=1;
a(3,9)=1; % Torque value
fgx=125.16; fgy=894.88; %Force from Ground
fa=-990;ta=68.05; %Force and angle of Acutuator
w2=1; a2=.75;
M2=MI(1)*-a2; 
% Mechanism Matrix Values [1,2]   
% 1=is Closed 2=Open 
W3=[ -0.41591,-1.665534];
W4=[-0.94692,-0.94692];
A3=[2.3066,-2.413351];
A4=[-2.4067,-1.152065];
T2=[155.78,329.93]; 
T3=[303.87,29.44];
T4=[349.85,263.9];
disp('What the forces are without considering Fi and Moment of Inertia')
for p=1:3
       if p==1 %Deployment 
       disp('Deployment')
       t2=T2(1); t3=T3(1); t4=T4(1); fgy=0;fgx=0; 
        elseif p==2 %retraction
           disp('Retraction')
           t2=T2(2); t3=T3(2); t4=(2); fa=-fa;
        elseif p==3 %landing 
           disp('Landing')
          t2=T2(2); t3=T3(2); t4=T4(2);  
          fgy=894.88; fgx=125.16; fa=0;
       end 
   r2x=r4op*cosd(t2);r2y=r4op*sind(t2);
   r3x=r3*cosd(t3);r3y=r3*sind(t3);
   r4x=r4*cos(t2); r4y=r4*cos(t2);
   %fix2=u(1)*Ax2;fiy2=u(1)*Ay2;
   %fix3=u(2)*Ax3;fiy3=u(2)*Ay3;
   %fix4=u(3)*Ax3;
   %fiy4=u(3)*Ay3;
   a(3,3)=-r2*sind(t2); a(3,4)=r2*cosd(t2); a(6,5)=-r3*sind(t3);a(6,6)=r3*cosd(t3);a(9,7)=-r4*sind(t4);a(9,8)=r4*cosd(t4);
  AA=array2table(a,'VariableNames',names,'RowNames',names)
   b=zeros(9,1);
   b(1)=-fa*cosd(ta);
   b(2)=fa*sind(ta);
   b(3)=-M2+fa*sin(ta)*r4x-fa*cos(ta)*r4y;%+fix4*.5*r4y+fiy4*.5*r4x;
   %b(4)=-fix3;
   %b(5)=-fiy3;
   %b(6)=-M3+fiy3*.5*r3x-fix3*.5*r3x;
   %b(7)=-fgx-fix2;
   b(8)=-fgy;%-fiy2;
   b(9)=fgy*r2x-r2y*fgx-M4;%+fiy2*.5*r2x-fix2*.5*r2y-M4;
   B=b
   c=a\b;
   C=array2table(c,'RowNames',names)
 
end

disp('What the forces are with considering Fi and Moment of Inertia')

for p=1:3
      if p==1 %Deployment 
        AAA={'Deployment'}; %Mechanism is closed
     disp(AAA)
      t2=155.78; t4=349.85; t3=303.87;  fa=-fa;
       w3=W3(1); w4=W3(1); a3=-2.413351; a4=-1.152065;
           Ax2=-r2*(sind(t2)*a2+cosd(t2*w2^2));
           Ay2=r2*(cosd(t2)*a2+sind(t2*w2^2));
           Ax3=-r3*(sind(t3)*a3+cosd(t3*w3^2));
           Ay3=r3*(cosd(t3)*a3+sind(t3*w3^2));
           Ax4=-r4*(sind(t2)*a4+cosd(t2*w4^2));
           Ay4=r4*(cosd(t2)*a4+sind(t2*w4^2));
           M3=MI(2)*-a3; M2=MI(1)*-a2;
   elseif  p==2 %retraction
           AAA={'Retraction'};
           disp(AAA)
           w3=-0.41591; w4=-0.94692; a3=2.3066; a4=-2.4067;  
           t2=329.93; t3=29.44; t4=263.9;  fgy=0;fgx=0; 
           Ax2=-r2*(sind(t2)*a2+cosd(t2*w2^2));
           Ay2=r2*(cosd(t2)*a2+sind(t2*w2^2));
           Ax3=-r3*(sind(t3)*a3+cosd(t3*w3^2));
           Ay3=r3*(cosd(t3)*a3+sind(t3*w3^2));
           Ax4=-r4*(sind(t2)*a4+cosd(t2*w4^2));
           Ay4=r4*(cosd(t2)*a4+sind(t2*w4^2));
           M3=MI(2)*-a3; M2=MI(1)*-a2;
   elseif p==3 %landing
           AAA={'Landing'};
           disp(AAA)
           % Acceleration of Mechanism is 0
          t4=263.9; t3=29.44; t2=329.93; fgy=894.88; fgx=125.16;
          Ax2=0;Ay2=0;Ax3=0;Ay3=0;Ax4=0;Ay4=0;
       end 
   
   r2x=r4op*cosd(t2);r2y=r4op*sind(t2);
   r3x=r3*cosd(t3);r3y=r3*sind(t3);
   r4x=r4*cos(t2); r4y=r4*cos(t2);
   fix2=u(1)*Ax2;fiy2=u(1)*Ay2;
   fix3=u(2)*Ax3;fiy3=u(2)*Ay3;
   fix4=u(3)*Ax3;
   fiy4=u(3)*Ay3;
   a(3,3)=-r2*sind(t2); a(3,4)=r2*cosd(t2); a(6,5)=-r3*sind(t3);a(6,6)=r3*cosd(t3);a(9,7)=-r4*sind(t4);a(9,8)=r4*cosd(t4);
  AA=array2table(a,'VariableNames',names,'RowNames',names)
   b=zeros(9,1);
   b(1)=-fa*cosd(ta);
   b(2)=fa*sind(ta);
   b(3)=-M2+fa*sin(ta)*r4x-fa*cos(ta)*r4y+fix4*.5*r4y+fiy4*.5*r4x;
   b(4)=-fix3;
   b(5)=-fiy3;
   b(6)=-M3+fiy3*.5*r3x-fix3*.5*r3x;
   b(7)=-fgx-fix2;
   b(8)=-fgy-fiy2;
   b(9)=fgy*r2x-r2y*fgx+fiy2*.5*r2x-fix2*.5*r2y-M4;
   disp(AAA)
 BB=b
 disp(AAA)
   answers=a\b;
   CC=array2table(answers,'RowNames',names)
end

    




  
