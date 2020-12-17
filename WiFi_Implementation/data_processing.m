%% Main part
clear all; close all; clc

Fs = 20; % Sampling frequency

ret = read_log_file('out_0610_2012');   % out_0508_1159 % out_0514_1046   
                                        % out_0518_1033 has heart rate from 4450 to 4900
[n_packets,~] = size(ret);

% choose which channel to use
nr1 = 1;
nc1 = 2;

nr2 = 1;
nc2 = 1;

% start time 
base =  ret{1, 1}.timestamp;

figure();
% find best subcarrier in chosen channel
for subcarrier = 1:4:56
    subplot(3,5,subcarrier / 4 + 0.75);
    for packet = 1:n_packets
        csi_data = ret{packet}.csi;
        if(csi_data ~= 0)
            data1(packet) = csi_data(nr1,nc1,subcarrier);
            data2(packet) = csi_data(nr2,nc2,subcarrier);
        end
    end

    ang1 = angle(data1);
    ang2 = angle(data2);  
    plot(1:n_packets, mod(ang1 - ang2, 2*pi), 'r');
end

%% Find best subcarrier

varSub = zeros(1,56);
data1 = zeros(56,n_packets);
data2 = zeros(56,n_packets);

for i = 1:56
    for packet = 1:n_packets
        csi_data = ret{packet}.csi;
        if(csi_data ~= 0)
            data1(i,packet) = csi_data(nr1,nc1,i);
            data2(i,packet) = csi_data(nr2,nc2,i);
        end
    end
    
    ang1Temp = angle(data1(i,:));
    
    ang2Temp = angle(data2(i,:));
    
    angTemp = mod(ang1Temp - ang2Temp, 2*pi);
    angTemp = hampel(angTemp, Fs*4/5,1.5);
    varSub(i) = sum(abs(angTemp-mean(angTemp)))./length(angTemp);
end
figure()
stem(varSub,'filled','Linewidth',1);
%title('Mean absolute deviation of each subcarrier');
xlabel('Subcarrier Index');
ylabel('Mean absolute deviation');
[value,index] = sort(varSub,'DESCEND');
best_sc = index(2);

%% One plot: fetch data with best subcarrier and best channel
data1 = zeros(1,n_packets);
data2 = zeros(1,n_packets);
%best_sc = 13;
for packet = 1:n_packets
    csi_data = ret{packet}.csi;
    if(csi_data ~= 0)
        data1(packet) = csi_data(nr1,nc1,best_sc);
        data2(packet) = csi_data(nr2,nc2,best_sc);
    end
end

ang1 = angle(data1);
ang2 = angle(data2); 

raw_data = mod(ang1 - ang2, 2*pi);

figure()
plot(1:n_packets, ang1,'r');
hold on;
plot(1:n_packets, ang2,'b');
hold on;
plot(1:n_packets, mod(ang1 - ang2, 2*pi),'y','Linewidth',2);
legend('Phase 1','Phase 2','Phase Difference');
xlabel('Number of Packets')
ylabel('Phase(rad)')
%title('Phase Difference Signal');

h_data = hampel(raw_data, Fs*4/5,1.5);
h_m_data = movmean(h_data,Fs);

figure();
plot(1:n_packets, raw_data, 'r');
hold on
plot(1:n_packets, h_data, 'b','Linewidth',3);
hold on
plot(1:n_packets, h_m_data, 'g','Linewidth',3);
legend('Raw Data','After Outlier Cancellation','After Noise Cancellation')
xlabel('Number of Packets')
ylabel('Phase Shift')
%title('CSI Phase shift vs Number of packets');
%hold on
%plot(h_m_data, 'ro');
%y_3 = xcorr(y_2,length(y_2)-1);
%[autocor,lags] = xcorr(h_m_data,'coeff');
%y3 = autocorr(h_m_data,length(h_m_data)-1);
%figure();
%plot(y3);
%xlabel('Lag (days)')
%ylabel('Autocorrelation')


%% analyze in frequency domain                    
T = 1/Fs;             % Sampling period       
X1 = floor(n_packets*0.2);
X2 = floor(n_packets*0.8);
L = X2-X1+1;        % Length of signal
t = (0:L-1)*T;        % Time vector

nfft = 2^nextpow2(L); % length of signal in power of 2
% remove DC component
ex_data = h_m_data(X1:X2);
data = ex_data-mean(ex_data);
data = detrend(data);
dcblker = dsp.DCBlocker;
data = dcblker(data');
%%%%%%%%%%%%%%%%%%%%%%%
fft_data = fft(data', nfft)/nfft; %FFT
fft_data1 = fft_data(1:nfft/2);  %take first half
x_fft = Fs*(0:(nfft/2)-1)/nfft;
figure();
plot(x_fft, abs(fft_data1),'LineWidth',2);
%title('Single-Sided Amplitude Spectrum of data');
xlabel('f(Hz)');
ylabel('Amplitude');
w_band = [1/5,1/3];
w = w_band ./ (Fs / 2);
[b, a] = butter(2,w,'bandpass');
H_bp = freqz(b,a,nfft/2);
figure()
plot(x_fft, abs(H_bp),'LineWidth',2);
xlabel('f(Hz)');
ylabel('Amplitude');
figure()
%filtered_data = filter(b,a,data');
filtered_data = filter(b,a,data);

plot(t,filtered_data,'LineWidth',3);
hold on;
plot(t,ex_data, 'r', 'LineWidth',3);
%title('Filtered data after Buttferworth bandpass filter');
xlabel('Time(s)');
ylabel('Phase difference');
legend('filtered data','unfiltered data')
time = L * T;
breath_butt = length(findpeaks(filtered_data))/time;

%% emd decomposition
[imf,residual] = emd(ex_data);
for i = 1:size(imf,2)
    rate = length(findpeaks(imf(:,i)))./time;
    if (rate<0.3) && (rate>0.15)
        breath_emd = rate;
    end
    if (rate<1.7) && (rate>0.5)
        heart_emd = rate;
    end
end


figure()
subplot(4,1,1);
plot(t,imf(:,1),'r','Linewidth',2);
title('Level 4')
hold on;
subplot(4,1,2);
plot(t,imf(:,2),'r','Linewidth',2);
title('Level 3')
hold on;
subplot(4,1,3);
plot(t,imf(:,3),'r','Linewidth',2);
title('Level 2')
hold on;
subplot(4,1,4);
plot(t,imf(:,4),'r','Linewidth',2);
title('Level 1')
hold on;
xlabel('Time(s)');


%% DWT
%temp2 = hampel(temp,2, 0.1);
temp2 = ex_data;

%load sumsin 
plot(temp2)
title('Signal')

for i = 1:floor(log2(Fs))
    [c,l] = wavedec(temp2,1,'db2');
    temp2 = appcoef(c,l,'db2');
    if (i == 1)
        cd1 = detcoef(c,l,1);
    end
    if (i == 2)
        cd2 = detcoef(c,l,1);
    end
    if (i == 3)
        cd3 = detcoef(c,l,1);
    end
    if (i == 4)
        cd4 = detcoef(c,l,1);
    end
    if (i == 5)
        cd5 = detcoef(c,l,1);
    end
    if (i == 6)
        cd6 = detcoef(c,l,1);
    end
    if (i == 7)
        cd7 = detcoef(c,l,1);
    end
end

%[c,l] = wavedec(ex_data,4,'db2');
%approx = appcoef(c,l,'db2');
%length(findpeaks(approx))./(length(approx) ./ Fs * 16)


subplot(5,1,1)
x = [0:length(temp2)-1]*T;
values = spcrv([[x(1) x x(end)];[temp2(1) temp2 temp2(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',2)
title('Approximation Coefficients')

subplot(5,1,2)
x = [0:length(cd4)-1]*T;
values = spcrv([[x(1) x x(end)];[cd4(1) cd4 cd4(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',2)
title('Level 4 Detail Coefficients')
subplot(5,1,3)
x = [0:length(cd3)-1]*T;
values = spcrv([[x(1) x x(end)];[cd3(1) cd3 cd3(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',2)
title('Level 3 Detail Coefficients')
subplot(5,1,4)
x = [0:length(cd2)-1]*T;
values = spcrv([[x(1) x x(end)];[cd2(1) cd2 cd2(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',2)
title('Level 2 Detail Coefficients')
subplot(5,1,5)
x = [0:length(cd1)-1]*T;
values = spcrv([[x(1) x x(end)];[cd1(1) cd1 cd1(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',2)
title('Level 1 Detail Coefficients')
xlabel('Time(s)');


breath_dwt = length(findpeaks(temp2))./(length(temp2) ./ Fs * (2^floor(log2(Fs))));

cd = ifft(fft(cd3,56)+fft(cd4,56));
x = [0:length(cd)-1]*T;
figure()
values = spcrv([[x(1) x x(end)];[cd(1) cd cd(end)]],3);
plot(values(1,:),values(2,:),'Linewidth',3);
xlabel('Time(s)');
ylabel('Phase difference');
%plot([0:length(cd)-1]*T,cd);
heart_dwt = length(findpeaks(cd))./(length(cd) ./ Fs * (Fs/2/1.875));

%% compute average breathing rate and heart rate

breath_rate = (breath_butt + breath_emd + breath_dwt)./3
heart_rate = (heart_emd + heart_dwt)./2