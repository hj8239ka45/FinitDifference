clc;clear;
L=30;W=50;
x = 1:L;
y = 1:W;
[xx, yy] = meshgrid(x, y);
yy=51-yy;
theta = 0;
for n=1:100
    x1 = sin(n*pi*xx/L);
    y1 = sinh(n*pi*yy/L)/sinh(n*pi*W/L);
    theta = theta + 2/pi*((-1)^(n+1)+1)/n.*x1.*y1;
end
[c,h]=contour(xx, yy, theta);
clabel(c,h);
title('theta'); 
xlabel('x/L') ;
ylabel('y/L');
axis tight