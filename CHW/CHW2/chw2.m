clc 
close 
%% p2
clc
clear
ramp = @(t) 2*(mod(t,2)-1.5).*(mod(t,2)>=3/2 & mod(t,2)<=2) + 2*(0.5-mod(t,2)).*(mod(t,2)>=0 & mod(t,2)<=1/2);
ak_ramp = FSC_calculator(ramp, 1000, 2, -10, 10);

Sawtooth_wave = @(t) (t - floor(t));
ak_Sawtooth_wave = FSC_calculator(Sawtooth_wave, 1000, 1, -10, 10);

my_fun = @(t) exp(2*mod(t,3)) + 2*(mod(t,3))^3;
ak_my_fun = FSC_calculator(my_fun, 1000, 3, -10, 10);
aks = {ak_ramp,ak_Sawtooth_wave,ak_my_fun};
names = {'Triangle wave','Sawtooth wave ',' exponential function' };

for i=1:length(aks)
    subplot(2,length(aks),i);
    ak_plotter(real(aks{i}));
    title(strcat('real part of ak coefiisents for ',names{i}),'Interpreter','latex','FontSize',10)
    xlabel('k','Interpreter','latex','FontSize',13)
    hold on;
    grid on;
    ylim([1.2*min(real(aks{i})) 1.2*max(real(aks{i}))]);

    subplot(2,length(aks),3+i);
    ak_plotter(imag(aks{i}));
    title(strcat('imaginary part of ak coefiisents for ',names{i}),'Interpreter','latex','FontSize',10)
    xlabel('k','Interpreter','latex','FontSize',13)
    hold on;
    grid on;
    ylim([1.2*min(imag(aks{i})) 1.2*max(imag(aks{i}))]);
 end
%% p3
clc
t = -5:0.01:5;

y = FS_calculator(ak_ramp , -10 , 10 ,2 ,t);
figure
subplot(3,1,1)
plot(t,y)
hold on
fplot(ramp)
title(strcat(names{1},'function'),'Interpreter','latex','FontSize',10)
xlabel('t','Interpreter','latex','FontSize',13)

y = FS_calculator(ak_Sawtooth_wave , -10 , 10 ,1 ,t);
subplot(3,1,2)
plot(t,y)
hold on 
fplot(Sawtooth_wave)
title(strcat(names{2},'function'),'Interpreter','latex','FontSize',10)
xlabel('t','Interpreter','latex','FontSize',13)

y = FS_calculator(ak_my_fun , -10 , 10 ,3 ,t);
subplot(3,1,3)
plot(t,y)
hold on
fplot(my_fun)
title(strcat(names{3},'function'),'Interpreter','latex','FontSize',10)
xlabel('t','Interpreter','latex','FontSize',13)

%% p4
Square_wave = @(t) (1+square(pi*(t/4+3/4),50))/2;
save_Square_Wave_gif(Square_wave , 30 , t)

%% Q2.1
clc
clear

num = 400102222;
digits = int32(cell2mat(cellfun(@str2num, split(num2str(num),""),UniformOutput=false)));
call = number_to_sound(digits,8192,100);
filename = 'StudentID.wav';
audiowrite(filename, call, 8192);
%% Q2.2
clc
clear

fs=8192;
N=1000;
plot_every_digit_fourier_trnasform(fs,N)

%% Q2.3
clc
clear
x1 = number_dector('dialing1');
x2 = number_dector('dialing2');
x3 = number_dector('dialing3');


%% Q2.4
clc
images = load_numbers_image();
show_number(images,x1,'dialing1')
show_number(images,x2,'dialing2')
show_number(images,x3,'dialing3')

%% Q2.5
clc
clear 
x1 = real_number_dector('realDialing1');
x2 = real_number_dector('realDialing2');
x3 = real_number_dector('realDialing3');
images = load_numbers_image();
show_number(images,x1,'Real Dialing1')
show_number(images,x2,'Real Dialing2')
show_number(images,x3,'Real Dialing3')




%% Q3.1
clc
clear 
[y,Fs]=audioread("Audio/sound.wav");

y=(y(:,1)+y(:,2))/2;
N=length(y);
x=Fs*(-N/2:N/2-1)/N;
dfft = fourier_transform(y , N);

figure
plot(x,dfft)
title('dfft for given sound','Interpreter','latex','FontSize',15)
xlabel('f','Interpreter','latex','FontSize',13)
%% Q3.2
mu = 0;
sigma = sqrt(0.0025);
r = normrnd(mu,sigma,[N,1]);
noisySound = y + r;
audiowrite("noisySound.wav",y,Fs);

%% Q3.3

my_filter=Filt();
recoverd_music=filter(my_filter,noisySound);

audiowrite("recoverdmusic.wav",recoverd_music,Fs);

%% Q3.4
figure
subplot(3,1,1)
dfft = fourier_transform(noisySound , N);
plot(x,dfft)
title('Discrete fourier transform of of noisy music')

subplot(3,1,2)
dfft = fourier_transform(y , N);
plot(x,dfft)
title('Discrete fourier transform of of music')


subplot(3,1,3)
dfft = fourier_transform(recoverd_music, N);
plot(x,dfft)
title('fft Discrete fourier transform of of recoverd music')

%% Q3.5
band_width=2300;

func_IdealFilter=@(w) (abs(w)<band_width);

arr_IdealFilter = arrayfun(func_IdealFilter,x);
arr_IdealFilter = arr_IdealFilter';

filteredMusic=arr_IdealFilter.*fftshift(fft(noisySound));
time_domain_filteredMusic=(ifft(ifftshift(filteredMusic)));
audiowrite("filteredNoisyMusic1.wav",time_domain_filteredMusic,Fs);


%%

filteredMusic2=conv(ifft(arr_IdealFilter),noisySound,'same');
audiowrite("filteredNoisyMusic2.wav",filteredMusic2,Fs);

figure
subplot(2,1,1)
N=length(time_domain_filteredMusic);
x=Fs*(-N/2:N/2-1)/N;
dfft = fourier_transform(time_domain_filteredMusic , N);
plot(x,dfft)
title('Discrete fourier transform of of recoverd music via dfft')

subplot(2,1,2)
N=length(filteredMusic2);
x=Fs*(-N/2:N/2-1)/N;
dfft = fourier_transform(filteredMusic2 , N);
plot(x,dfft)
title('Discrete fourier transform of of recoverd music via convolution')


%% Functions

function Hd = Filt()
    %FILT Returns a discrete-time filter object.

    % Butterworth Bandpass filter

    % All frequency values are in Hz.
    Fs = 48000;  % Sampling Frequency

    N   = 10;
    Fc1 = 50;
    Fc2 = 1200;

    h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
    Hd = design(h, 'butter');
end

function images = load_numbers_image()
    number0 = imread("Numbers\0.png");
    number1 = imread("Numbers\1.png");
    number2 = imread("Numbers\2.png");
    number3 = imread("Numbers\3.png");
    number4 = imread("Numbers\4.png");
    number5 = imread("Numbers\5.png");
    number6 = imread("Numbers\6.png");
    number7 = imread("Numbers\7.png");
    number8 = imread("Numbers\8.png");
    number9 = imread("Numbers\9.png");
    images = {number0 ,number1 ,number2, number3, number4, number5, number6, number7, number8, number9};
end

function show_number(images,x1,name)
    h = figure;
    set(h,'Position',[450 360 800 420]);
    for i=1:length(x1)
        h = subplot(1,length(x1),i);
        imshow(images{x1(i)+1})
        pos = get(h,'OuterPosition');
        pos(3) = 0.103; % set width to 100% of figure
        set(h,'OuterPosition',pos);
    end
    pos = get(h,'OuterPosition');
    pos(3) = 0.113; % set width to 100% of figure
    set(h,'OuterPosition',pos);
    sgtitle(name)

end



function numbers = real_number_dector(str)
    [y,fs]=audioread("Audio/"+str+".wav");
    [runstarts,runends] = real_zero_remover(y); 
    y = y';
    runends = runends';
    runstarts = runstarts';
    runends= [1  runends];
    runstarts = [runstarts length(y)];
    numbers = [];
    n = length(runstarts);
    for i = 1 : n-1
        my_fft=fftshift(fft(y(runends(i):runstarts(i))));
        N = runstarts(i) - runends(i);
        fft_oneside=my_fft(1:N);
        dfft=abs(fft_oneside)/(N/2);
        [maxVals, maxIdxs] = maxk(dfft, 4);
        number = match_number_frequency(maxIdxs,N);
        numbers = [numbers number];
    end
    numbers = int8(numbers);
end






function numbers = number_dector(str)
    [y,fs]=audioread("Audio/"+str+".wav");
    N = 1000;
    y = zero_remover(y);
    y = y';
    n = size(y);
    n = n(2)/1000;
    numbers = [];
    for i = 1 :n
        my_fft=fftshift(fft(y((i-1)*1000+1:1000*i)));
        fft_oneside=my_fft(1:N);
        dfft=abs(fft_oneside)/(N/2);
        [maxVals, maxIdxs] = maxk(dfft, 4);
        number = match_number_frequency(maxIdxs,N);
        numbers = [numbers number];
    end
    numbers = int8(numbers);
end

function number = match_number_frequency(maxIdxs,N)
    maxIdxs = maxIdxs - N/2;
    cindition = find(maxIdxs < 0);
    c = setdiff(1:length(maxIdxs), cindition);
    positiveMaxIdxs = maxIdxs(c);
    frequencies = positiveMaxIdxs*(8192/N);
    number = match_frequency_to_number(frequencies);

end

function number = match_frequency_to_number(frequencies)
    row_frequencies = [697 770 852 941];
    column_frequencies = [1209 1336 1477];
    for i=1:2
        row_frequencies1 = row_frequencies - frequencies(1);
        row_frequencies2 = row_frequencies - frequencies(2);
    end
    [~, idx] = min(abs([row_frequencies1 row_frequencies2]));
    match_row_frequency = mod(idx,4);
    for i=1:2
        column_frequencies1 = column_frequencies - frequencies(1);
        column_frequencies2 = column_frequencies - frequencies(2);
    end
    [~, idx2] = min(abs([column_frequencies1 column_frequencies2]));
    match_column_frequency = mod(idx2,3);
    if(match_column_frequency==0)
        match_column_frequency=3;
    end
    number = (match_row_frequency-1)*3 + match_column_frequency;



    if(match_column_frequency == 2 && match_row_frequency ==0)
        number = 0;
    end
end


function y = zero_remover(x)
    transitions = diff([0 ;x == 0; 0]); %find where the array goes from non-zero to zero and vice versa
    runstarts = find(transitions == 1);
    runends = find(transitions == -1); %one past the end
    runlengths = runends - runstarts;
    % keep only those runs of length 100 or less:
    runstarts(runlengths < 100) = [];
    runends(runlengths < 100) = [];
    %expand each run into a list indices:
    indices = arrayfun(@(s, e) s:e-1, runstarts, runends, 'UniformOutput', false);
    indices = [indices{:}];  %concatenate the list of indices into one vector
    % Remove those zeros which are consecuitve 3 in number
    x(indices) = [] ;
    y = x;
    %credit : https://www.mathworks.com/matlabcentral/answers/399882-removing-consecutive-zeros-of-a-certain-length

end

function [runstarts,runends] = real_zero_remover(x)
    transitions = diff([0 ;x == 0; 0]);
    runstarts = find(transitions == 1);
    runends = find(transitions == -1); 
    runlengths = runends - runstarts;
    runstarts(runlengths < 100) = [];
    runends(runlengths < 100) = [];
end


function plot_every_digit_fourier_trnasform(fs,N)
    x=fs*(-N/2:N/2-1)/N;
    figure
    for i=0:9
        subplot(4,3,i+1);
        if(i==9)
            subplot(4,3,[10,11,12])
        end
        plot_fourier_digits(x,i,N,fs)
        title('dfft for'+string(i),'Interpreter','latex','FontSize',10)
        xlabel('f','Interpreter','latex','FontSize',13)
    end
end

function plot_fourier_digits(x,number,N,fs)
    func = digit_to_sound(number);
    n=0:1/fs:(N-1)/fs;
    xd = arrayfun(func,n);
    dfft = fourier_transform(xd , N);
    plot(x,dfft)
end

function dfft = fourier_transform(signal , N)
    my_fft=fftshift(fft(signal));
    fft_oneside=my_fft(1:N);
    dfft=abs(fft_oneside)/(N/2);
end




function call = number_to_sound(x,Fs,space)
    t1 = 0:1/Fs:999/Fs;
    spaces = zeros(1, space);
    call = zeros(1,0);
    for i=1:length(x)
        y = digit_to_sound(x(i));
        y = arrayfun(y,t1);
        call = cat(2 , call , y);
        if(i<length(x))
            call = cat(2 , call , spaces);
        end
    end
end



function y = digit_to_sound(x)
    y = @(t) 0;
    a = @(t) 0;
    b = @(t) 0;
    if(mod(x,3) ==1)
        a = @(t) cos(2*pi*1209*t);
    end

    if(mod(x,3) == 2)
        a = @(t) cos(2*pi*1336*t);
    end
    if(mod(x,3) == 0)
        a = @(t) cos(2*pi*1477*t);
    end

    if(x<=3)
        b = @(t) cos(2*pi*697*t);
    end

    if(x<=6 && x>3)
        b = @(t) cos(2*pi*770*t);
    end

    if(x<=9 && x>6 )
        b = @(t) cos(2*pi*852*t);
    end
    y = @(t) a(t) + b(t);
    if(x==0)
        y = @(t) cos(2*pi*1336*t) + cos(2*pi*941*t);
    end
end

function ak_plotter(ak)
    stem(-10:10,ak);
end





function save_Square_Wave_gif(func , K , t)
filename = 'FS_Simulation.gif';
for k = 1:K
    clf
    fplot(func)
    ylim([-0.3 , 1.3])
    grid on
    hold on
    title('Square Wave function for k =' + string(k),'Interpreter','latex','FontSize',10)
    xlabel('t','Interpreter','latex','FontSize',13)
    ak = FSC_calculator(func, 1000, 8, -k, k);
    plot(t,FS_calculator(ak , -k , k ,8 ,t))
    pause(0.4)
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
end


function ak = FSC_calculator(xc, Fs, T, k_min, k_max)
    N = floor(Fs * T);
    n = linspace(0,N-1,N)/Fs;
    xd = arrayfun(xc,n);
    ak0 = fftshift(fft(xd,N)/N);
    ak = ak0(N/2+k_min+1:N/2+k_max+1);
end

function y = FS_calculator(ak, k_min, k_max, T, t)
    y=0;
    for i=1:k_max-k_min+1
        y = y + ak(i)*exp(((2*pi/T)*(i+k_min-1)*1i).*t);
    end
end
