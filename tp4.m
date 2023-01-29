clear all
close all
clc

Te = 5*1e-4;
f1 = 500;
f2 = 400;
f3 = 50;
t = 0:Te:5-Te;
fe = 1/Te;
N = length(t);

fshift = (-N/2:N/2-1)*(fe/N);
f = (0:N-1)*(fe/N);

x = sin(2*pi*f1*t)+sin(2*f2*pi*t)+sin(2*pi*f3*t);
y = fft(x);
%{
subplot(2,1,1)
plot(t,x)

subplot(2,1,2)
plot(fshift, fftshift(abs(y)));
%}

k = 1;
wc = 2*pi*250;
wc1 = 2*pi*500;
wc2 = 2*pi*1000;

h = (k*1j*((2*pi*f)/wc))./(1+1j*((2*pi*f)/wc));
h1 = (k*1j*((2*pi*f)/wc1))./(1+1j*((2*pi*f)/wc1));
h2 = (k*1j*((2*pi*f)/wc2))./(1+1j*((2*pi*f)/wc2));

G = 20*log(abs(h));
G1 = 20*log(abs(h1));
G2 = 20*log(abs(h2));

P =  angle(h);
P1 = angle(h1);
P2 = angle(h2);

subplot(3,1,1)
semilogx(f,G)
legend("Module de h(t)")
grid on

% subplot(2,1,1)
% semilogx(f,G,f,G1,f,G2,"Linewidth",1.5);
% title("Diagramme de Bode")
% xlabel("rad/s")
% ylabel("decibel")
% legend(" fc=50","fc2=500","fc3=1000")
% grid on
% 
% subplot(2,1,2)
% semilogx(f,P,f,P1,f,P2,"Linewidth",1.5)
% grid on
% %legend("P","P1","P2")
% %filtrage

filter =  [h(1:floor(N/2)),flip(h(1:floor(N/2)))];
  
filtred_signal = y.*filter; 

temp_signal = ifft(filtred_signal,"symmetric");

plot(fshift,fftshift(2*abs(fft(temp_signal))/N))


[music,fs] = audioread('test.wav');
music = music';
N = length(music);

spect_music = fft(music);


f = (0:N-1)*(fs/N);
fshift = (-N/2:N/2-1)*(fs/N);

plot(fshift, fftshift(abs(spect_music)));

fc = 4500; 
H = 1./(1+1j*(f/fc).^1); 

filter = [H(1:floor(N/2)),flip(H(1:floor(N/2)))]; 

filterd_signal_freq = spect_music(1:end-1).*filter; 

filterd_signal_tmp = ifft(filterd_signal_freq,"symmetric");

plot(fshift(1:end-1), fftshift(abs(fft(filterd_signal_tmp))));


semilogx(f(1:floor(N/2)),abs(H(1:floor(N/2))),"linewidth",1.5