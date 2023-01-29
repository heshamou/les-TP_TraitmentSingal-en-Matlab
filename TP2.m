clear all
close all
clc
 
%Question 1:

[x,fe]=audioread("C:\Users\zaina\Desktop\MyFile\S'5\Traitement de signal\tp2\Phrase.m4a");

%Question 2:

sound(x,fe)
plot(x)

%Question 3:

sound(x,fe*2)
sound(x,fe/2)
sound(x,fe)



%Question 4:

stem(x)
%Après avoir tracé la fonction x, nous pouvons voir que la pièce
%Rien ne sert est contenu dans la partie de l'indice 1025 jusqu'à 2.05e05

%Question 5:


riennesertde=x(1025:(2.05e05));

%sound(riennesertde,fe)

%Question 6:

courir=x((2.239e+05):(2.642e+05));
%sound(courir,fe)

ilfaut=x((2.987e+05):(3.81e+05));


partirapoint=x((4e+05):(5.28e+05));




%Question 7:
parole=[riennesertde ,partirapoint ,ilfaut,courir];



% Synthèse et analyse spectrale d’une gamme de musique

%Synthèse d’une gamme de musique

%Question 1:

fe=8192;
te=1/fe;
N=5000;
t=(0:N-1)*te;
do=10*cos(2*pi*262*t);
%sound(do,fe)
re=10*cos(2*pi*294*t);
%sound(re,fe)
mi=10*cos(2*pi*330*t);
%sound(mi,fe)
fa=10*cos(2*pi*349*t);
%sound(re,fa)
sol=10*cos(2*pi*392*t);
%sound(sol,fe)
la=10*cos(2*pi*440*t);
%sound(la,fe)
si=10*cos(2*pi*494*t);
%sound(si,fe)
do2=10*cos(2*pi*523*t);
%sound(do2,fe)
musique=[do,re,mi,fa,sol,la,si,do2];
sound(musique,fe)
 

%Spectre de la gamme de musique

%Question 2:

f=(0:N-1)*(fe/N);
spectre_musique=fft(musique);
signalAnalyzer(abs(fftshift(spectre_musique)));
