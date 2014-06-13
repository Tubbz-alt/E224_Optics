D=3.718; %Distance Axicon-DS Gold Mirror
alpha=0.006021;%Laser angle of convergence (theoretical)

figure(1)
clf();
hold on;

rds_max=D*alpha-0.004765;
x_max=rds_max*cos(linspace(0,2*pi,1000));
y_max=rds_max*sin(linspace(0,2*pi,1000));
plot(x_max,y_max,'r');

rds_min=D*alpha-0.019;
x_min=rds_min*cos(linspace(0,2*pi,1000));
y_min=rds_min*sin(linspace(0,2*pi,1000));
plot(x_min,y_min,'r');

%Cutting from P1
for a=linspace(0.005,0.018,20)
x=ones(1,1000)*a;
ymax=sqrt(0.019^2-a^2);
y=linspace(-ymax,ymax,1000);
r_cut=x+i*y;
transfo_cut=r_cut.*(1-D*alpha./(abs(r_cut)));
plot(transfo_cut);
end

%Cutting from P3

D_P3=3.690; %Distance P3-DS Gold Mirror
R_P3=0.0125; %Rayon de P3
%Distance P3-Laser
%delta_x=-0.013;
delta_y=0;

for delta_x=linspace(-0.027,-0.015,20)
t=linspace(-pi,pi,1000);
P3_cut=delta_x+i*delta_y + R_P3*(1/sqrt(2)*cos(t)+i*sin(t));
P3_cut_f=P3_cut(abs(P3_cut)<=0.0188);
%plot(P3_cut_f,'m');
transfo_P3cut=P3_cut_f.*(1-D_P3*alpha./(abs(P3_cut_f)));
plot(transfo_P3cut,'g');
end