function [PPGpks,mins,PPG_pk_locs,min_locs,numPPGPeaks] = noise_removal(Fs,fPPG_IR,oPPGpks,oPPG_pk_locs,oMins,oMin_locs,N,minPDist)
numPPGPeaks = length(oPPGpks);
numPPGMins = length(oMins);

%Vectors for verified peaks, mins (removed artifacts)
raPPGpks = zeros(numPPGPeaks,1);        
ramins = zeros(numPPGPeaks,1);

%Binary array for points to be ignored in the future (found as artifacts)
ignorePts = zeros(N,1);

%Remove outlier peaks and mins (noise/motion artifacts)
%Peaks
c=0;
for i=1:numPPGPeaks
    currPeak = oPPGpks(i);
    currPeakLoc = oPPG_pk_locs(i);
    startPeak = max(1,i-3);
    endPeak = min(numPPGPeaks,i+3);
    %Check the 3 previous and 3 following peaks (not including current peak)
    wPks = [oPPGpks(startPeak:endPeak)];
    avgPeakHeight = mean(wPks) - currPeak/length(wPks);
    
    %Verify that peak is not noise based on other peaks in the window
    if currPeak>0 && (currPeak<2*avgPeakHeight || currPeak<0.5)
        c=c+1;
        raPPGpks(c) = currPeak;
    else %Noise
        ignorePts(currPeakLoc) = 1;
    end
end
%Keep reference vectors, but copy them to new vectors (for comparison later on)
PPGpks = raPPGpks(1:c,1);
%Mins
c=0;
for i=1:numPPGMins
    currMin = oMins(i);
    currMinLoc = oMin_locs(i);
    startMin = max(1,i-3);
    endMin = min(numPPGMins,i+3);
    %Check the 3 previous and 3 following mins (not including current min)
    wMins = [oMins(startMin:endMin)];
    avgMinHeight = mean(wMins) - currMin/length(wMins);
    
    %Verify that min is not noise based on other mins in the window
    if currMin<0 && currMin>2*avgMinHeight
        c=c+1;
        ramins(c) = currMin;
    else %Noise
        ignorePts(currMinLoc) = 1;
    end
end
mins = ramins(1:c,1);

%Use findpeaks for a final set of peaks and mins based on peaks with noise
%removed
minPkHeight = 0.4*median(PPGpks);
minMinHeight = -0.4*median(mins);

[PPGpks, PPG_pk_locs] = findpeaks(fPPG_IR,'MinPeakDistance',minPDist,'MinPeakHeight',minPkHeight);
%Minimums found through peak detection of the negative wave
[mins, min_locs] = findpeaks(-fPPG_IR,'MinPeakDistance',minPDist,'MinPeakHeight',minMinHeight);
mins = -mins;                         %Use negative to get back to original value

%Finally, remove the points that were classified before as noise
%remove the points that were classified before as noise
i=1;
while i <= length(PPGpks)
    if ignorePts(PPG_pk_locs(i)) == 1
        PPGpks(i) = [];
        PPG_pk_locs(i) = [];
    end
    i=i+1;
end
i=1;
while i <= length(mins)
    if ignorePts(min_locs(i)) == 1
        mins(i) = [];
        min_locs(i) = [];
    end
    i=i+1;
end

%% Removing Non-Pair PPG Peaks and Mins for SpO2 Calculations
%Remove the first min if it is before a max
if min_locs(1) < PPG_pk_locs(1)
    min_locs = min_locs(2:end);
    mins = mins(2:end);
end
%Remove the last max if it has no following min
if min_locs(end) < PPG_pk_locs(end)
    PPG_pk_locs = PPG_pk_locs(1:end-1);
    PPGpks = PPGpks(1:end-1);
end
numPeaks = length(PPGpks);

tpks = zeros(numPeaks,1);
tmins = zeros(numPeaks,1);
tPPG_pk_locs = zeros(numPeaks,1);          
tmin_locs = zeros(numPeaks,1);          %Validated peaks and mins stored here

%Determine the threshold for PPG max-min distance (to validate pairs)
MaxToMinDist = zeros(numPeaks,1);
for i = 1:numPeaks
    smallestDist = N;
    for j = 1:length(mins)
        dist = min_locs(j) - PPG_pk_locs(i);
        if dist > 0 && dist < smallestDist 
            smallestDist = dist;
        end
    end
    %Only add to vector if a valid peak-min distance was found
    if smallestDist < N
        MaxToMinDist(i) = smallestDist;
    end
end
%Eliminate zeros from vector and find average to use as a threshold
MaxToMinDist = MaxToMinDist(MaxToMinDist>0);
MaxToMinThresh = min(1.5*median(MaxToMinDist),Fs*1);    %Maximum Threshold (#samples) between each peak and min                         
o=0;                                                    %Offset value for alignment
c=0;                                                    %Index counter for true peaks/mins

for i = 1:numPeaks
    %Ensure that scanning only occurs within bounds
    if i+o > length(min_locs)
        break;
    end
    
    %Missed peak(s), one or more minimums has come before
    %Loop until the min is caught back up to recent max
    while i+o <= length(min_locs) && min_locs(i+o) < PPG_pk_locs(i)
        o=o+1;                          %Increase offset for alignment
    end
    
    %Difference (#samples) between max and corresponding min must be within the thresh
    if min_locs(i+o) - PPG_pk_locs(i) <= MaxToMinThresh
        c=c+1;                          
        tPPG_pk_locs(c) = PPG_pk_locs(i);       
        tmin_locs(c) = min_locs(i+o);
        tpks(c) = PPGpks(i);
        tmins(c) = mins(i+o);           %Add all validated data to corresponding array
    else
        o=o-1;                          %Decrease offset for alignment
    end
end
%Only the true peak and min pairs remain
PPGpks = tpks(1:c,1);
mins = tmins(1:c,1);
PPG_pk_locs = tPPG_pk_locs(1:c,1);
min_locs = tmin_locs(1:c,1);
numPPGPeaks = c;
end

