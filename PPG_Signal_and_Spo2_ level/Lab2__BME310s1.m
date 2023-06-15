% This code reads the PPG and ECG data from the .CSV files and estimate
% Average SpO2 level in the blood
clear
clc

%% File and Data Acquisition
% All 200 Hz, Fingertip 
f1 = "Lab2_BME310_dataset.csv";

%% Choose File and Convert Spreadsheet Data into Table
chosenFile = f1;
data = readtable(chosenFile); 

%Number of samples to remove from start and end to clean the time-data sample
sC = 300;                             %Start cutoff (#samples)
eC = 0;                               %End cutoff (#samples)
%Apply cutoffs
N = height(data);
data = data(sC+1:N-eC,:);
N = N-sC-eC;
fN =1; %File specific number used later for reference values

%% Wave Variable Setups
%Sampling frequency, time intervals
Fs = 200;                               %Sample rate, Hz
tInt = 1/Fs;                            %Time interval between data points (s)
tElapsed = tInt*(N-1);                  %Total time elapsed for recording (s)
t = linspace(0,tElapsed,N)';            %Time vector (s)

%Data obtained from the table MAX86150 board
%PPGs are inverted to positively correlate with varying blood volume*
PPG_IR = -data.Var3;                    %IR Light Data
PPG_Red = -data.Var4;                   %Red Light Data

%% Filtering the PPG waveforms
[Fd,fPPG_IR,fPPG_Red,sPPG_IR,sPPG_Red] = PPG_filter(Fs,N,PPG_IR,PPG_Red);

%% Find PPG Peaks
T_s = 1/Fd;                             %Calculate dominant wave period (s)
T_n = T_s*Fs;                           %Calculate dominant wave period (#samples)
minPDist = 0.8*T_n;

%Perform peak detection using findpeaks
[oPPGpks, oPPG_pk_locs] = findpeaks(fPPG_IR,'MinPeakDistance',minPDist);
%Minimums found through peak detection of the negative wave
[oMins, oMin_locs] = findpeaks(-fPPG_IR,'MinPeakDistance',minPDist);
oMins = -oMins;                         %Use negative to get back to original value

%% 
[PPGpks,mins,PPG_pk_locs,min_locs,numPPGPeaks] = noise_removal(Fs,fPPG_IR,oPPGpks,oPPG_pk_locs,oMins,oMin_locs,N,minPDist);

%% Plot Raw and filtered PPG signals 
figure(1)
subplot(2,1,1)
plot(1:N,detrend(PPG_IR,1),'b',1:N,detrend(PPG_Red,1),'r')
title('Original PPGs')
legend('Detrended Raw IR PPG','Detrended Raw Red PPG')
ylabel('Amplitude')
xlabel('Sample #')
xlim([0 1000])

subplot(2,1,2)
plot(1:N,sPPG_IR,'b',1:N,sPPG_Red,'r')
title('Filtered PPGs')
legend('Filtered IR PPG','Filtered Red PPG')
ylabel('Amplitude')
xlabel('Sample #')
xlim([0 1000])
