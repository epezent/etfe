clear all
close all
clc

fs = 1000;
n    = 10000;
w    = 2000;
o    = w/2;
win  = hamming(w);
k    = (o - w + n) / o;
nfft = 2000;
pcf = 1 / (fs*sum(win.^2));

t = 0:1/fs:10-1/fs;
x = 0.5 * sin(2*pi*25*t) + 0.5 * sin(2*pi*50*t) + randn(size(t));
[b,a] = butter(2,0.2);
y = filter(b,a,x);

% preallocate the space for the psd results
mypsd = zeros(nfft/2+1, k);
mycsd = zeros(nfft/2+1, k);
% step 1: loop through the data, 512 points at a time, with 256 points overlap
for i = 0:k-1
    % step 2: apply a hamming window
    temp1 = x(1+o*i:w+i*o)'.*win;
    temp2 = y(1+o*i:w+i*o)'.*win;
    % step 3: calculate FFT and take just the first half
    temp1 = fft(temp1,nfft);
    temp2 = fft(temp2,nfft);
    % step 4: calculate the "periodogram" by taking the absolute value squared
    temp2 = temp1 .* conj(temp2);
    temp1 = abs(temp1).^2;
    % save the results in the storage variable
    mypsd(:, i+1) = temp1(1:(nfft/2+1));   
    mycsd(:, i+1) = temp2(1:(nfft/2+1));   
end
% step 5: take the average of all the periodograms
mypsd = mean(mypsd,2);
mycsd = mean(mycsd,2);
% normalizing factor
mypsd = mypsd * pcf;
mycsd = mycsd * pcf;
% ignore the DC and Nyquist value
mypsd(2:end-1) = mypsd(2:end-1) * 2;
mycsd(2:end-1) = mycsd(2:end-1) * 2;
% find tfe
mytfe = conj(mycsd) ./ mypsd;
% find freqs
f = 0:fs/nfft:fs/2;

figure
subplot(5,1,1)
hold on
plot(t,x)
plot(t,y)
legend('x','y')
grid on

subplot(5,1,2)
pwelch(x, win, [], nfft, fs);
hold on
plot(f, 10*log10(mypsd),'r--');
legend('welch','clone')
grid on

subplot(5,1,3)
cpsd(x, y, win, [], nfft, fs);
hold on
plot(f, 10*log10(abs(mycsd)),'r--');
legend('cpsd','clone')
grid on

subplot(5,1,4)
tfestimate(x, y, win, [], nfft, fs)
hold on
plot(f, 20*log10(abs(mytfe)),'r--')

subplot(5,1,5)
hold on
txy = tfestimate(x, y, win, [], nfft, fs);
plot(f, 180/pi*angle(txy))
plot(f, 180/pi*angle(mytfe),'r--')
legend('tfestimate','clone')


