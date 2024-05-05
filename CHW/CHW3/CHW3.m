clear
clc
close all
%%
clc
close all
ecg = importdata("CHW3\ecg.mat");
ecg1 = ecg(:,1)';
ecg2 = ecg(:,2)';
plot_ecg('ECG1',ecg1,360)
plot_ecg('ECG2',ecg2,360)

%%
clc
Fs = 360;
Hd1 = Build_BandPass_Filter(55,65,Fs);
fvtool(Hd1);
%%
filteredECG1 = filter(Hd1, ecg1);
filteredECG2 = filter(Hd1, ecg2);

plot_ecg('Filtered ECG1',filteredECG1,360)
plot_ecg('Filtered ECG2',filteredECG2,360)

%%
Hd2 = HighPassFilter(1 , 3,Fs);
fvtool(Hd2);
%%
filteredECG1 = filter(Hd2, ecg1);
filteredECG2 = filter(Hd2, ecg2);

plot_ecg('Filtered Baseline Wander ECG1',filteredECG1,360)
plot_ecg('Filtered Baseline Wander ECG2',filteredECG2,360)
%%
Hd3 = dfilt.cascade(Hd1, Hd2);
fvtool(Hd3);
%%
filteredECG1 = filter(Hd3, ecg1);
filteredECG2 = filter(Hd3, ecg2);

plot_ecg('Filtered Baseline Wander and urban electricity noise',filteredECG1,360)
plot_ecg('Filtered Baseline Wander and urban electricity noise',filteredECG2,360)

%%%%%%






%%
clear
clc 
close all
img = imread("sample_19201280.bmp");
imshow(img)

% Convert the RGB image to YCbCr color space
img_ycbcr = rgb2ycbcr(img);

% Extract the Y, Cb, and Cr channels
Y = img_ycbcr(:,:,1);
Cb = img_ycbcr(:,:,2);
Cr = img_ycbcr(:,:,3);

y_downsampled = imresize(Y, 0.5, 'nearest');
cb_downsampled = imresize(Cb, 0.5, 'nearest');
cr_downsampled = imresize(Cr, 0.5, 'nearest');

% replace the original Cb component with the downsampled version
img_ycc(:,:,1) = y_downsampled;
img_ycc(:,:,2) = cb_downsampled;
img_ycc(:,:,3) = cr_downsampled;

% convert back to RGB color space
img_downsampled = ycbcr2rgb(img_ycc);
figure
imshow(img_downsampled)
%%
double_Cb = double(cb_downsampled) - 128;
double_Cr = double(cr_downsampled) - 128;
double_y = double(y_downsampled) - 128;



dct2_Cb = block_dct(double_Cb);
dct2_Cr = block_dct(double_Cr);
dct2_Y = block_dct(double_y);
%%
clc
K = [1 5 10 20 50 80];
Quantisized_Cr = Quantisize(false,dct2_Cr,K);
Quantisized_Cb = Quantisize(false,dct2_Cb,K);
Quantisized_Y = Quantisize(true,dct2_Y,K);
%%
clc
new_y = 256*dct2_inverse(Quantisized_Y)+128;
new_Cr = 256*dct2_inverse(Quantisized_Cr)+128;
new_Cb = 256*dct2_inverse(Quantisized_Cb)+128;

%%
figure
for i =1:6
    subplot(2,3,i)
    img_ycc(:,:,1) = reshape(new_y(i,:,:),[640 960]);
    img_ycc(:,:,2) = reshape(new_Cb(i,:,:),[640 960]);
    img_ycc(:,:,3) = reshape(new_Cr(i,:,:),[640 960]);
    img_new = ycbcr2rgb(img_ycc);
    imshow(img_new)
    title(strcat("Quality =" , int2str(K(i))));
    sgtitle("JPG Image With Diffrent Qualities")
    imwrite(img_new, strcat(strcat("IMG with Quality =" , int2str(K(i))),".jpg"),"jpg")
end

%%
function plot_ecg(name,ecg,sampling_rate)
    n = length(ecg);
    t = n/sampling_rate;
    x = 0:1/sampling_rate:t-1/sampling_rate;
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(x,ecg);
    xlim([0 t]);
    xlabel('S',Interpreter='latex',Color=[0.25, 0.25, 0.25],FontSize=17);
    ylabel('mV',Interpreter='latex',Color=[0.25, 0.25, 0.25],FontSize=17)
    title(append(name,' Time Domain'),Interpreter='latex',Color=[0.15, 0.15, 0.15],FontSize=28)
    grid on
    max_y = max(abs(min(ecg)),max(ecg));
    ylim([-1.5*max_y 1.5*max_y])

    [f,func_fft] = plot_furier_transform(ecg);

    figure('units','normalized','outerposition',[0 0 1 1])
    plot(f ,func_fft)
    xlabel('f',Interpreter='latex',Color=[0.25, 0.25, 0.25],FontSize=17);
    ylabel('abs',Interpreter='latex',Color=[0.25, 0.25, 0.25],FontSize=17)
    title(append(name,' Frequency Domain'),Interpreter='latex',Color=[0.15, 0.15, 0.15],FontSize=28)
    grid on
    ylim([1.5*min(func_fft) 1.5*max(func_fft)])
end

function [f,func_fft] = plot_furier_transform(func)
    func_fft = fft(func);
    f_m = length(func_fft)/2;
    f = -f_m:f_m-1;
    func_fft = abs(fftshift(func_fft));
end

function Hd = Build_BandPass_Filter(fstop1,fstop2,Fs)
    Fpass1 = fstop1-7;        % First Passband Frequency
    Fstop1 = fstop1;          % First Stopband Frequency
    Fstop2 = fstop2;          % Second Stopband Frequency
    Fpass2 = fstop2+7;        % Second Passband Frequency
    Dpass1 = 0.025;           % First Passband Ripple
    Dstop  = 0.001;           % Stopband Attenuation
    Dpass2 = 0.055;           % Second Passband Ripple
    dens   = 20;              % Density Factor

    % Calculate the order from the parameters using FIRPMORD.
    [N, Fo, Ao, W] = firpmord([Fpass1 Fstop1 Fstop2 Fpass2]/(Fs/2), [1 0 ...
                          1], [Dpass1 Dstop Dpass2]);

    % Calculate the coefficients using the FIRPM function.
    b  = firpm(N, Fo, Ao, W, {dens});
    Hd = dfilt.dffir(b);
end

function Hd = HighPassFilter(Fstop , Fpass,Fs)
    Dstop = 0.0001;          % Stopband Attenuation
    Dpass = 0.055;           % Passband Ripple
    dens  = 20;              % Density Factor

    % Calculate the order from the parameters using FIRPMORD.
    [N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);

    % Calculate the coefficients using the FIRPM function.
    b  = firpm(N, Fo, Ao, W, {dens});
    Hd = dfilt.dffir(b);

end

function arr = block_dct(double_arr)
    [h, w] = size(double_arr);
    num_blocks_h = floor(h / 8);
    num_blocks_w = floor(w / 8);
    dct_blocks = zeros(8, 8, num_blocks_h, num_blocks_w);
    for i = 1:num_blocks_h
        for j = 1:num_blocks_w
            block = double_arr((i-1)*8+1:i*8, (j-1)*8+1:j*8);
            dct_block = dct2(block);
            dct_blocks(:, :, i, j) = dct_block;
        end
    end
    arr = dct_blocks;
end

function Quantisized = Quantisize(lum,dct2_arr,K)
    [a, b, h, w] = size(dct2_arr);
    Quantisized = zeros(length(K),a,b,h,w);
    alpha = 1;
    for quality = K
        disp(quality)
        if(lum)
                    Quantisize_mat = [16 11 10 16 24 40 51 61 ;
                          12 12 14 19 26 58 60 55 ;
                          14 13 16 24 40 57 69 56 ;
                          14 17 22 29 51 87 80 62 ;
                          18 22 37 56 68 109 103 77 ;
                          24 35 55 64 81 104 113 92 ;
                          49 64 78 87 103 121 120 101 ;
                          72 92 95 98 112 100 103 99];

        else
                        Quantisize_mat = [17 18 24 47 99 99 99 99 ;
                              18 21 26 66 99 99 99 99 ;
                              24 26 56 99 99 99 99 99 ;
                              47 66 99 99 99 99 99 99 ;
                              99 99 99 99 99 99 99 99 ;
                              99 99 99 99 99 99 99 99 ;
                              99 99 99 99 99 99 99 99 ;
                              99 99 99 99 99 99 99 99 ];

        end
    
        if(quality > 50)
                Quantisize_mat = Quantisize_mat * (100 - quality) / 50;
        else
                Quantisize_mat = Quantisize_mat * 50 / quality;
        end
         dct2_arr(:,:,1:h,1:w)= dct2_arr(:,:,1:h,1:w) ./ Quantisize_mat;
         Quantisized(alpha,:,:,:,:)=dct2_arr;
         alpha=alpha+1;
    end
    
end

function temp = dct2_inverse(Quantisized)
    [a,b,c,d,e] = size(Quantisized);
    temp=zeros(a,b*d,c*e);
    for i=1:a
        for j =1:d
            for l=1:e
                temp(i,8*(j-1)+1:8*j,8*(l-1)+1:8*l)=reshape(idct2(Quantisized(i,:,:,j,l)),[8 8]);
            end
        end
    end
    
end
