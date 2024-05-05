clc;
close all;
clear;

%% Q1.1

n0 = 0;
n1 = 1;
n = -10:10;
delta = @(n) ((n-n0)==0);
figure

myheaviside = @(n) ((n>=n0)==1);

subplot(3,2,1)
x1 = 4*delta(n+5) + 2*delta(n)-delta(n-3);
stem(n,x1,'filled')
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_1[n] = 4\delta[n+5] + 2\delta[n]-\delta[n-3]$','Interpreter','latex','Color','black','FontSize',10)
ylim([-0.5 4.5])
grid on

subplot(3,2,2)
x2 = (-0.8).^n .* myheaviside(n);
stem(n,x2,'filled')
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_2[n] = (-0.8)^n u[n]$','Interpreter','latex','Color','black','FontSize',10)
grid on
ylim([-1.2 1.2])

subplot(3,2,3)
x3 = myheaviside(n-5) - myheaviside(n+5);
stem(n,x3,'filled')
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_3[n] = u[n-5]- u[n+5]$','Interpreter','latex','Color','black','FontSize',10)
ylim([-1.2 1.2])
grid on

subplot(3,2,4)
x4 = cos(0.5*pi*n) + 3*sin(pi*n);
stem(n,x4,'filled')
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_4=cos(0.5\pi n) +3sin(\pi n)$','Interpreter','latex','Color','black','FontSize',10)
ylim([-1.2 1.2])
grid on

subplot(3,2,[5,6])
x5 = 1.*n;
x5(1) = x1(1)+ x2(1);
for i = 2:21
    x5(i) = x5(i-1) + x1(i)+x2(i);
end
stem(n,x5,"filled")
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_5=\sum_{i=-\infty}^n(x_1+x_2)[i]$','Interpreter','latex','Color','black','FontSize',10)
ylim([-1 8])
grid on
%%

%% Q1.2
clc;
close all;
clear;

t =-5:0.01:5;
r = t.*heaviside(t);
figure('Renderer', 'painters', 'Position', [500 60 700 800]);
subplot(6,1,1)
x1 = (1-exp(-t)).*heaviside(t);
plot(t,x1)
xlabel('t','Interpreter','latex','Color','black','FontSize',10)
title('$x_1(t) = 1-e^{-t}u(t)$','Interpreter','latex','Color','black','FontSize',10)
grid on

subplot(6,1,[2,3])
x2 = 3*(t+3).*heaviside(t+3) +heaviside(t) - heaviside(t-2) ;
plot(t,x2)
xlabel('t','Interpreter','latex','Color','black','FontSize',10)
title('$x_2(t) = 3r(t+3)+u(t)-u(t-2)$','Interpreter','latex','Color','black','FontSize',10)
grid on


subplot(6,1,4)
x3 = sinc(t).*(heaviside(t+3) - heaviside(t-3));
plot(t,x3)
xlabel('t','Interpreter','latex','Color','black','FontSize',10)
title('$x_3(t) = sinc(t)(u(t+3)-u(t-3))$','Interpreter','latex','Color','black','FontSize',10)
grid on

subplot(6,1,[5,6])
syms t
x4 = 0.4*cos(pi*t)+1;
x5 = x4*cos(30*pi*t);
fplot(x4)
hold on
fplot(x5)
legend("X_{4} = 0.4cos(\pi t)+1","X_{5} = x_4 cos(30\pi t)")
xlabel('t','Interpreter','latex','Color','black','FontSize',10)
title('$x_1(t) \; \; \& \; \;x_5(t)$','Interpreter','latex','Color','black','FontSize',10)
grid on

%% Q2
clc
n = -20:20;
n2 = -40:40;
n0=0;
delta = @(n) ((n-n0)==0);
x1 = heaviside(n+7) - heaviside(n-8) +0.5*delta(n) - 4* delta(n-3) +2*delta(n+4);
x2 = (3/5).^n .* (heaviside(n+12) - heaviside(n-12));
x3 = n.*sin(0.5*n).*( heaviside(n + 15) - heaviside(n - 15));

tic
y11 = myConv(x1,x2);
toc

tic
y12 = conv(x1,x2);
toc

figure 

subplot(2,1,1)
stem(n2,y11)
xlim([-20 20])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_1 * x_2 \; by \; myConv$','Interpreter','latex','Color','black','FontSize',10)
grid on

subplot(2,1,2)
stem(n2,y12)
xlim([-20 20])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_1 * x_2 \; by \; conv$','Interpreter','latex','Color','black','FontSize',10)
grid on


% 
tic
y21 = myConv(x2,x3);
toc

tic
y22 = conv(x2,x3);
toc

figure 

subplot(2,1,1)
stem(n2,y21)
xlim([-28 28])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_2 * x_3 \; by \; myConv$','Interpreter','latex','Color','black','FontSize',10)
grid on

subplot(2,1,2)
stem(n2,y22)
xlim([-28 28])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_2 * x_3 \; by \; conv$','Interpreter','latex','Color','black','FontSize',10)
grid on
%
tic
y31 = myConv(x1,x3);
toc

tic
y32 = conv(x1,x3);
toc

figure 

subplot(2,1,1)
stem(n2,y31)
xlim([-23 23])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_1 * x_3 \; by \; myConv$','Interpreter','latex','Color','black','FontSize',10)
grid on

subplot(2,1,2)
stem(n2,y32)
xlim([-23 23])
xlabel('n','Interpreter','latex','Color','black','FontSize',10)
title('$x_2 * x_3 \; by \; conv$','Interpreter','latex','Color','black','FontSize',10)
grid on
%%
A = [1 -1 1;-1 1 -1;1 -1 1];
B = [1 2 3 4 5;6 7 8 9 10;11 12 13 14 15;16 17 18 19 20;21 22 23 24 25];

ConvFull = conv2(A,B);
ConvSame = conv2(A,B,'same');

%% Q3
clc
clear
img = imread('pic.jpg');
imshow(img);

R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);

figure;
subplot(2,2,1);
imshow(img)
title('Orginal image')

subplot(2,2,2);
Rimg = img;
Rimg(:,:,2:3) =0;
imshow(Rimg)
title('Red Channel');

subplot(2,2,3);
Gimg = img;
Gimg(:,:,1) =0;
Gimg(:,:,3) =0;
imshow(Gimg)
title('Green Channel');

subplot(2,2,4);
Bimg = img;
Bimg(:,:,1:2) =0;
imshow(Bimg)
title('Blue Channel');

figure

subplot(1,2,1)
RGB_sum = uint8(mean(cat(3, R, G, B),3));
imshow(RGB_sum)
title('Mean of R,G,B matrices')

subplot(1,2,2)
I = rgb2gray(img);
imshow(I)
title('Grayscale of image done by rgb2gray')


figure
subplot(1,2,1)
img_double = im2double(RGB_sum);
imshow(img_double);
title("convert wises to double via im2double")
subplot(1,2,2)
imshow(I)
title('Grayscale of image done by rgb2gray')

EDKx = [-1 0 1;-2 0 2;-1 0 1];
EDKy = [-1 -2 -1;0 0 0;1 2 1];

IMGedgesx = conv2(EDKx,img_double);
IMGedgesy = conv2(EDKy,img_double);
IMGedge = sqrt((IMGedgesx.^2) + (IMGedgesy.^2));

figure

subplot(2,2,1)
imshow(RGB_sum)
title('Mean of R,G,B matrices')
subplot(2,2,2)
imshow(IMGedgesx)
title('horizontal edges of image')
subplot(2,2,3)
imshow(IMGedgesy)
title('vertical edges of image')
subplot(2,2,4)
imshow(IMGedge)
title('magnitude of edges')

img_double = im2double(img);

Rdouble = img_double(:,:,1);
Gdouble = img_double(:,:,2);
Bdouble = img_double(:,:,3);

R_x = conv2(Rdouble, EDKx, 'same');
R_y = conv2(Rdouble, EDKy, 'same');
G_x = conv2(Gdouble, EDKx, 'same');
G_y = conv2(Gdouble, EDKy, 'same');
B_x = conv2(Bdouble, EDKx, 'same');
B_y = conv2(Bdouble, EDKy, 'same');

R_mag = sqrt(R_x.^2 + R_y.^2);
G_mag = sqrt(G_x.^2 + G_y.^2);
B_mag = sqrt(B_x.^2 + B_y.^2);

img_double(:,:,1) = R_mag;
img_double(:,:,2) = G_mag;
img_double(:,:,3) = B_mag;

figure
imshow(img_double)
title('Color edge detected image')

%%
clc

BBK = 1/9*ones(3);
Bimg1 = imconv(img , BBK);

figure
subplot(1,2,1)
imshow(img)
title('Original image')
subplot(1,2,2)
imshow(Bimg1)
title('Slightly blurred image')

figure

kernel_size = 15;
kernel = ones(kernel_size) / kernel_size^2;
Bimg2 = imconv(img , kernel);
subplot(3,1,1)
imshow(img)
title('Original image')
subplot(3,1,2)
imshow(Bimg1)
title('Slightly blurred image')
subplot(3,1,3)
imshow(Bimg2)
title('Blurred image')

ISK = [0 -1 0; -1 5 -1; 0 -1 0];
SharpenImg = imconv(img,ISK);
figure
subplot(1,2,1)
imshow(img)
title('Original image')
subplot(1,2,2)
imshow(SharpenImg)
title('Sharpen image')
%% Q4.1
clc
clear
y = audioread('Music.wav');
fs = 48000;
figure
stem(y)
mean_y = mean(y,2);

%soundsc(y,fs)
%pause(22)
%soundsc(mean_y,fs)
%% Q4.2
alpha=0.8;
nd=0.2*fs;
echoMusix = filter([1 ,zeros(1,nd) , alpha] , 1 , mean_y);
soundsc(echoMusix,fs)
%% Q4.3
RecoveredMusic = filter(1,[1 ,zeros(1,nd) , alpha],echoMusix);
soundsc(RecoveredMusic,fs);
%% Q4.4
alpha=0.8;
nd=0.1*fs;
f = zeros(1,3*nd+1);
f(1)=1;
f(nd)=alpha;
f(2*nd)=alpha^2;
f(3*nd)=alpha^3;

echoMusix = filter(f , 1 , mean_y);
soundsc(echoMusix,fs)
RecoveredMusic = filter(1,f,echoMusix);
pause(22)
soundsc(RecoveredMusic,fs)

figure
subplot(3,1,1)
stem(mean_y)
title('Original music')
grid on
subplot(3,1,2)
stem(echoMusix)
title('Echo music')
grid on
subplot(3,1,3)
stem(RecoveredMusic)
title('Recoverd music')
geid on
%% Q4.5
Mean=0;
Var=0.01;
z1 = imnoise(mean_y,'Gaussian',Mean,Var);
%soundsc(z1,fs)

u = -0.1 + rand(length(mean_y),1)*(0.2);
z2 = u+mean_y;
figure
subplot(2,1,1)
stem(z1)
title('Normal distribution noise')
subplot(2,1,2)
stem(z2)
title('Uniform distribution noise')
soundsc(z2,fs)
%% Q4.6
clc
t = 0:1/fs:(length(mean_y)-1)/fs; 
f = linspace(1000,2000,length(t)); 
y = 0.1*sin(2*pi.*f.*t);
y = transpose(y);
MusicSine = mean_y + y;
soundsc(MusicSine,fs)
audiowrite('MusicSine.wav',MusicSine,fs)

%% Functions
function y = myConv(x, h)
y = zeros(1, length(x) + length(h) - 1);
for i = 1:length(y)
    for j = 1:length(h)
        if i-j+1 > 0 && i-j+1 <= length(x)
            y(i) = y(i) + x(i-j+1) * h(j);
        end
    end 
end 
end

function Bimg = imconv(img , kernel)
BimgR = conv2(kernel,im2double(img(:,:,1)));
BimgG = conv2(kernel,im2double(img(:,:,2)));
BimgB = conv2(kernel,im2double(img(:,:,3)));

Bimg = cat(3,BimgR,BimgG,BimgB);
end

