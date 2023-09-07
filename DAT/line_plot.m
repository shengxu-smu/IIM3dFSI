clear all

load xc.dat
load yc.dat
load zc.dat

load xu.dat
load yv.dat
load zw.dat

load cxcycz.dat
ks=1;
frc=[cxcycz(:,1),cxcycz(:,2+(ks-1)*9:1+9*ks)];

figure(1)
plot(xc,xu,'bo-',yc,yv,'ro-',zc,zw,'go-')
legend('u(x)','v(y)','w(z)')

figure(2)
plot(frc(:,1),frc(:,2),'b-',frc(:,1),frc(:,3),'r-',frc(:,1),frc(:,4),'g-')
legend('cx','cy','cz')

figure(3)
plot(frc(:,1),frc(:,6),'b-',frc(:,1),frc(:,7),'r-',frc(:,1),frc(:,8),'g-',frc(:,1),frc(:,5),'k-')
legend('ra','rb','rc','vol')

figure(4)
plot(frc(:,1),frc(:,9),'b-')
legend('zsc')

figure(5)
plot(frc(:,1),frc(:,10),'r-')
legend('zsct')


clear all

