clear;clc;close all;

%使用者自訂
x = 0.00125;%x取樣-->0.02m(or 0.01)
t=0.025;  %設定時間間隔為0.1
Time = 10;%取10秒內結果

%基本參數
y = x;   %y取樣-->0.02m
W = 0.04;%寬
L = 0.08;%長
sample_y = 1 + W/y;
sample_x = 1 + L/2/x; %把銅和鋼分開取


%建立 mesh indep.test
for i=1:sample_y*sample_x
    T(i)=700;
end
k_x1=sample_y*sample_y+1;
k_x2=sample_y*(sample_y+1);
k_y1=sample_y*(2*sample_y-1)+1;
k_y2=sample_y*sample_y*2;
T(k_x1:k_x2)=1000;
T(k_y1:k_y2)=800;
den=2;
n=(0.02/x);
pow=0;
while (n)>=1
    n=n/2;
    pow=pow+1;
end
n=pow;
for num=1:sample_y-2
    k_z1=1/den^n*(num*(k_x1)+(sample_y-1-num)*k_y1)
    k_z2=1/den^n*(num*(k_x2)+(sample_y-1-num)*k_y2)
    T(k_z1:k_z2)=1/den^n*(num*T(k_x1:k_x2)+(sample_y-1-num)*T(k_y1:k_y2));
end


R = 2*10^-5;%contact resistance m^2*K/W
%AISI_1010
k_st=31.3; %熱傳導係數 w/(m*K)
c_st=1168; %熱容係數   j/(kg*K)
rho_st=7832; %密度       kg/m^3
alpha_st=k_st/(rho_st*c_st);%     m^2/s
Fo_st=t*alpha_st/x^2;%傅立葉數
Bi_st=y/(R*k_st);   %biot number

if (1/(2*Bi_st+4)-Fo_st)>0    %判別Fourier number 合理性
    "Fo_st為合理值"
else
    "Fo_st不合理，重新設定題目"  %因為Fo為題目給定值
end

%COPPER 設定幾乎同上，改為銅之參數
k_c=(379+366)/2;
c_c=(417+433)/2;
rho_c=8960;
alpha_c=k_c/(rho_c*c_c);
Bi_c=y/(R*k_c);
Fo_c=t*alpha_c/x^2;

if (1/(2*Bi_c+4)-Fo_c)>0    %判別Fourier number 合理性
    "Fo_c為合理值"
else
    "Fo_c不合理，重新設定t間隔"
end


%有限差分法矩正運算
%setting C
C=-4*eye(sample_y,sample_x);
for i=1:sample_y
    for j=1:sample_x
        if i==j-1||i==j+1
            C(i,j)=1;
        end
    end
end
C(1,2)=2;C(sample_y,sample_x-1)=2;

%setting D
D = zeros(sample_y,sample_x);
for i=1:sample_y
    if i~=1&&i~=sample_y
        k_y1 = ((i-1)*sample_y+1);  k_y2 = i*sample_y;      %上下界運算
        k_x1 = (i-2)*sample_y+1;    k_x2 = (i+1)*sample_y;  %上下界運算
        D(k_y1:k_y2,k_x1:k_x2) = horzcat(eye(sample_y,sample_x),C,eye(sample_y,sample_x));
    end
end
D(1:sample_y,1:2*sample_x) = horzcat(C,2*eye(sample_y,sample_x));
k_y1=(sample_y-1)*sample_y+1;   k_y2=sample_y*sample_y; %上下界運算
k_x1=(sample_x-2)*sample_x+1;   k_x2=sample_x*sample_x; %上下界運算
D(k_y1:k_y2,k_x1:k_x2) =horzcat(2*eye(sample_y,sample_x),C);

Fo_R=0;     %Fo分布設定
for i=1:sample_y^2*2
    for j=1:sample_y^2*2
        if i<=sample_y^2
            Fo_R(i,j)=Fo_c;
        else
            Fo_R(i,j)=Fo_st;
        end
    end
end
%setting A
A = Fo_R.*[D zeros(sample_y^2,sample_x^2);zeros(sample_y^2,sample_x^2) D];     %整合有限差分法運算
%setting E
E_c =2*Fo_c*Bi_c*eye(sample_y,sample_x);
E_st=2*Fo_st*Bi_st*eye(sample_y,sample_x);
%setting B
B = zeros(sample_y^2*2,sample_x^2*2);
k_y1=(sample_y-1)*sample_y+1;   k_y2=sample_y*sample_y; %上下界運算
k_x1=(sample_x-1)*sample_x+1;   k_x2=(sample_x+1)*sample_x; %上下界運算
B(k_y1:k_y2,k_x1:k_x2) = horzcat(-E_c,E_c);
k_y1=sample_y*sample_y+1;   k_y2=(sample_y+1)*sample_y; %上下界運算
B(k_y1:k_y2,k_x1:k_x2) = horzcat(E_st,-E_st);
I=eye(sample_y^2*2,sample_x^2*2);%單位矩陣-->用來考量上一時間點的數據

%% explicit
T1=I+A+B; %整合有限差分法(外顯)的邊界矩正以及上一時間點資料
Temp1=T';        %溫度時間資料(溫度隨時間疊加)
T_time1 = zeros(sample_y*sample_x*2,Time/t);%(18筆)溫度資料*(100筆)時間資料
T_time1(:,1) = Temp1(:,1);%第一筆不隨時間改變
for n=2:1:Time/t
    Temp1=T1*Temp1;
    T_time1(:,n) = Temp1(:);
end
%% implicit
T2=I-A-B; %整合有限差分法(內顯)的邊界矩正以及上一時間點資料
Temp2=T';    %溫度時間資料(溫度隨時間疊加)
T_time2 = zeros(sample_y*sample_x*2,Time/t);%(18筆)溫度資料*(100筆)時間資料
T_time2(:,1) = Temp2(:,1);%第一筆不隨時間改變
for n=2:1:Time/t
    Temp2=inv(T2)*Temp2;
    T_time2(:,n) = Temp2(:);
end

T_data1=zeros(sample_y,sample_x*2,Time/t);
T_data2=zeros(sample_y,sample_x*2,Time/t);
k=1;
for i=1:sample_y
    k=i;
    for j=1:sample_x*2
        T_data1(i,j,:)=T_time1(k,:);
        T_data2(i,j,:)=T_time2(k,:);
        k=k+sample_x;
    end
end

tt=t:t:Time/t*t;
figure(1);
subplot(3,2,1);
T_plot(:,:)=T_data1((sample_y+1)/2,:,:);
T_plot2(:,:)=T_data2((sample_y+1)/2,:,:);
plot(tt,T_plot(1,:));hold on;plot(tt,T_plot2(1,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T1,2");
xlabel("Time(s)");
ylabel("Temp(K)");
legend("Explit.","Implit.",'Location','northeast');
subplot(3,2,2);%1+(sample_x-1)/2
plot(tt,T_plot(1+(sample_x-1)/2,:));hold on;plot(tt,T_plot2(1+(sample_x-1)/2,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T2,2");
xlabel("Time(s)");
ylabel("Temp(K)");
subplot(3,2,3);%1+(sample_x-1)/2*2
plot(tt,T_plot(1+(sample_x-1)/2*2,:));hold on;plot(tt,T_plot2(1+(sample_x-1)/2*2,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T3,2");
xlabel("Time(s)");
ylabel("Temp(K)");
subplot(3,2,4);%2+(sample_x-1)/2*2
plot(tt,T_plot(2+(sample_x-1)/2*2,:));hold on;plot(tt,T_plot2(2+(sample_x-1)/2*2,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T4,2");
xlabel("Time(s)");
ylabel("Temp(K)");
subplot(3,2,5);%2+(sample_x-1)/2*3
plot(tt,T_plot(2+(sample_x-1)/2*3,:));hold on;plot(tt,T_plot2(2+(sample_x-1)/2*3,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T5,2");
xlabel("Time(s)");
ylabel("Temp(K)");
subplot(3,2,6);%2+(sample_x-1)/2*4
plot(tt,T_plot(2+(sample_x-1)/2*4,:));hold on;plot(tt,T_plot2(2+(sample_x-1)/2*4,:));
axis([ 0, Time,700 ,1000]);
title("Theta scheme for T6,2");
xlabel("Time(s)");
ylabel("Temp(K)");


% 動畫製作
% figure(2);
% i=1:sample_x*2;
% j=1:sample_y;
% [XX,YY]=meshgrid(i*x,j*y);
% for tt=1:Time/t
%     subplot(2,1,1);
%     [c,h]=contour(XX,YY,T_data1(:,:,tt),[700:20:1000]);%等溫線圖-->[700:20:1000]畫線區間
%     
%     hold on;
%     %streamline
%     u = zeros(sample_y,sample_x*2);
%     for a=1:sample_x*2
%         if a~=1
%             u(:,a) = (-T_data1(:,a,tt)+T_data1(:,a-1,tt))/sample_x;
%         end
%     end
%     v = zeros(sample_y,sample_x*2);
%     quiver(XX,YY,u,v)
%     hold off;
%     
%     clabel(c,h);%標示線的溫度數值
%     title("Theta distribution (Explit.)");
%     xlabel("Time(s)");
%     ylabel("Temp(K)");
%     subplot(2,1,2);
%     [c,h]=contour(XX,YY,T_data2(:,:,tt),[700:20:1000]);%等溫線圖-->[700:20:1000]畫線區間
%     
%     hold on;
%     %streamline
%     u = zeros(sample_y,sample_x*2);
%     for a=1:sample_x*2
%         if a~=1
%             u2(:,a) = (-T_data2(:,a,tt)+T_data2(:,a-1,tt))/sample_x;
%         end
%     end
%     v = zeros(sample_y,sample_x*2);
%     quiver(XX,YY,u2,v)
%     hold off;
%     
%     clabel(c,h);%標示線的溫度數值
%     title("Theta distribution (Implit.)");
%     xlabel("X(m)");
%     ylabel("Y(m)");
%     % <動畫抓取
%     frames(1)=getframe(gcf);
%     [image,map]=frame2im(frames(1));
%     [im,map2]=rgb2ind(image,128);
%     % <gif製作 -->在同一資料夾中建立heattransfer_project.gif的gif檔案
%     if tt==1
%         imwrite(im,map2,'heattransfer_project.gif','gif','writeMode','overwrite','delaytime',0.1,'loopcount',inf);
%     else
%         imwrite(im,map2,'heattransfer_project.gif','gif','writeMode','append','delaytime',0.1);
%     end
%     % /gif製作>
%     % /動畫抓取>
% end

% %% 比教semi-infinite
% if x==0.02
%     T_semi = [700,700,700,700,700,700,700,700,700,800,800,800,800,800,800,800,800,800];
% elseif x==0.01
%     T_semi = [700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,700,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800,800];
% end
% T_sem(:,1) = T_semi;
% T_s = ((k_c*rho_c*c_c)^0.5*700+(k_st*rho_st*c_st)^0.5*800)/((k_c*rho_c*c_c)^0.5+(k_st*rho_st*c_st)^0.5)
% for n=2:1:Time/t
%     for i=1:sample_x^2
%         T_sem(i,n) = erf(fix((sample_x^2-i)/sample_x)*x/(4*alpha_c*n))*(700-T_s)+T_s;
%     end
%     for i=sample_x^2+1:sample_x^2*2
%         T_sem(i,n) = erf(fix((i-sample_x^2-1)/sample_x)*x/(4*alpha_c*n))*(800-T_s)+T_s;
%     end
%     %T_semi_1(10:18,n) = erf(x/(4*alpha_st*n))*(800-700)+700;
% end
% 
% for i=1:sample_y
%     k=i;
%     for j=1:sample_x*2
%         T_sem_data(i,j,:)=T_sem(k,:);
%         k=k+sample_x;
%     end
% end
% T_sem_plot(:,:) = T_sem_data((sample_y+1)/2,:,:);
% xx = x:x:2*sample_x*x;
% figure
% for i=1:Time/t
%     plot(xx,T_sem_plot(:,i));
%     hold on;
% end
k=0;
for j= 1:sample_x*2
    k=k+1;
    if j==sample_x
        k=k-1;
    else
        T_data_semi(k,:) = T_data1(3,j,:);
    end
end
xx = x:x:(2*sample_x-1)*x;
figure
for i=1:Time/t
    plot(xx,T_data_semi(:,i));
    hold on;
end
