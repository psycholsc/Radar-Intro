% clc
close all
clear all
%% para
f_c = 10e9;
dur = 10e-6;
bw  = 10e6;
f_s = 100e6;
tar = 3000;
R_max=6000;
c=3e8;
%% waveform
t=0:1/f_s:(R_max*2/c-1/f_s);
k=bw/dur;
s_ref=rectpuls(t/dur/2).*exp(1j*pi*k*t.*t);%.*exp(-1j*2*pi*f_c*t);
%% 
figure;
subplot(211)
plot(real(s_ref))
xlabel('sample (1/f_s) ')
ylabel('amplitude')
title('real part of LFM signal (time domain)')
grid on
subplot(212)
plot(imag(s_ref))
xlabel('sample (1/f_s) ')
ylabel('amplitude')
title('imag part of LFM signal (time domain)')
grid on
figure;
plot((t-max(t)/2)/2e-5*50,abs(abs(fftshift(fft(s_ref)))))
xlabel('Frequency (MHz)')
ylabel('amplitude')
title('amplitude of LFM signal (freq domain)')
grid on
%%
delay = tar*2/c;
% s_recv = rectpuls((t-delay*1.125)/dur*2).*exp(1j*2*pi*k*(t-delay).*(t-delay)).*exp(-1j*2*pi*f_c*delay);
s_recv = rectpuls((t-delay*1.25)/dur).*exp(1j*pi*k*(t-delay).*(t-delay)).*exp(-1j*2*pi*f_c*delay);
% cor=xcorr(s_recv,s_ref);
% figure;plot(10*log(abs(cor)))
figure;
subplot(211)
plot(real(s_recv))
xlabel('sample (1/f_s) ')
ylabel('amplitude')
title('real part of received signal (time domain)')
grid on
subplot(212)
plot(imag(s_recv))
xlabel('sample (1/f_s) ')
ylabel('amplitude')
title('imag part of received signal (time domain)')
grid on
figure;
plot((t-max(t)/2)/2e-5*50,abs(abs(fftshift(fft(s_recv)))))
xlabel('Frequency (MHz)')
ylabel('amplitude')
title('amplitude of received signal (freq domain)')
grid on

%%
s_ref_f=fft(s_ref);
s_recv_f=fft(s_recv);
corr=conj(s_ref_f).*s_recv_f;
corr_t=ifft(corr);
corr_t=corr_t/max(corr_t);
% figure;plot(abs(corr))
figure;
plot(20*log10(abs(corr_t)))
xlabel('sample (1/f_s) ')
ylabel('magnitude (dB)')
title('magnitude of filtered signal')
grid on
[~,index]=max(abs(corr_t));
axis([index-100 index+100 -50 0])
hold on;
plot(20*log10(abs(sqrt(0.5)*ones(1,numel(t)))))
plot(20*log10(abs(2/3/pi*ones(1,numel(t)))))