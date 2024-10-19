clc;
clear;
n = -1: 0.01: 1;
a = 1.5;
f = 4;
x_1 = a * sin(2*pi*f*n); % 正弦序列
n = -1: 0.02: 1;
a = 0.1;
x_2 = a.^n; % 实指数序列

x_conv = conv(x_1,x_2); % 卷积

zx = n;
zy = n;
for n=1:length(zx) 
    for m=1:length(zy) 
        z(n,m)=zx(n)+1i*zy(m); 
        Zt_tmp_1=0; 
        for ii=1:length(x_1) 
            Zt_tmp_1=Zt_tmp_1+x_1(ii)*((z(n,m))^(-(ii-1))); 
        end 
        Zt_tmp_2=0; 
        for ii=1:length(x_2) 
            Zt_tmp_2=Zt_tmp_2+x_2(ii)*((z(n,m))^(-(ii-1))); 
        end 
        Zt_conv=0;
        for ii=1:length(x_conv) 
            Zt_conv=Zt_conv+x_conv(ii)*((z(n,m))^(-(ii-1))); 
        end 
        Zt_1(n,m)=Zt_tmp_1;
        Zt_2(n,m)=Zt_tmp_2;
        Zt_c(n,m)=Zt_conv;
        z_product=Zt_tmp_1*Zt_tmp_2;
        Zt(n, m) = z_product;
    end 
end

figure(1);
view(3);
subplot(2,1,1);
mesh(zx,zy,abs(Zt));
title('Z变换的乘积')
xlabel('Re(z)');
ylabel('Im(z)');

subplot(2,1,2);
mesh(zx,zy,abs(Zt_c));
title('卷积的Z变换')
xlabel('Re(z)');
ylabel('Im(z)');

% subplot(2,1,2);
% mesh(zx,zy,abs(Zt_2));
% xlabel('Re(z)');
% ylabel('Im(z)');

% w = (-1: 0.005: 1) * pi;
% for n=1:length(w) 
%     Zw_tmp_1=0; 
%     Zw_tmp_2=0; 
%     for ii=1:length(x_1) 
%         Zw_tmp_1=Zw_tmp_1+x_1(ii)*(exp(-1i*w(n)*(-(ii-1)))); 
%     end 
%     for ii=1:length(x_2) 
%         Zw_tmp_2=Zw_tmp_2+x_2(ii)*(exp(-1i*w(n)*(-(ii-1)))); 
%     end 
%     zw_1(n)=Zw_tmp_1; 
%     zw_2(n)=Zw_tmp_2; 
% end

% figure(2);
% title('频率响应');
% subplot(2,1,1);
% plot(w/pi,abs(zw_1));
% xlabel('f/pi');
% ylabel('Magnitude');

% subplot(2,1,2);
% plot(w/pi,abs(zw_2));
% xlabel('f/pi');
% ylabel('Magnitude');
