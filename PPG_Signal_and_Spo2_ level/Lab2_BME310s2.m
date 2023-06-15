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

%% SpO2 Calculations
%% Calculating SpO2 Ratio from PPG Peaks and Mins (AC/DC values)
%Use original, unfiltered waves for the DC (offset) values
%Use filtered waves for the AC values (peak-min distance)
DC_IR = zeros(numPPGPeaks,1);           %Store DC values at every pulse
DC_Red = zeros(numPPGPeaks,1);
AC_IR = zeros(numPPGPeaks,1);           %Store AC values at every pulse
AC_Red = zeros(numPPGPeaks,1);
R = zeros(numPPGPeaks,1);               %Store Ratio of AC/DCs for SpO2 Calc.
for n = 1:numPPGPeaks                   %Loop through each pulse (i.e. peak-min pair)
    %Red Pulse
    DC_Red(n) = PPG_Red(min_locs(n));
    AC_Red(n) = sPPG_Red(PPG_pk_locs(n))-sPPG_Red(min_locs(n)); 

    %IR Pulse
    DC_IR(n) = PPG_IR(min_locs(n));
    AC_IR(n) = sPPG_IR(PPG_pk_locs(n))-sPPG_IR(min_locs(n));
    
    %Calculate R
    R(n) = [fill in the blank]; % use Equation 2
end

%% SpO2 Formulas from Ratio (calibration-less)
mean(R);
SpO2 = [fill in the blank]; %use Equation 1

%Remove overshot values above 100% SpO2 if nec.
SpO2 = SpO2(SpO2<100);
avgSpO2 = mean(SpO2);
SpO2str = ['Average SpO2 level = ',num2str(round(avgSpO2,2)),'%']






