function mseed_cross_corr_2()

clc;
close all;

%% =========================================================
% Parameters
% =========================================================

targetChannel = 'BHZ';
fLow  = 0.05;
fHigh = 10.0;
filterOrder = 4;
shiftSeconds = 600;
frequencyList = 0.2:0.02:3;       % 0.20 -1.0 Hz

fprintf('\n============================================================\n');
fprintf('MiniSEED Goertzel cross-spectrum calculation\n');
fprintf('============================================================\n');

%% =========================================================
% Select directory
% =========================================================

dataDir = uigetdir(pwd,'Select MiniSEED directory');
if isequal(dataDir,0)
    return;
end

files = dir(fullfile(dataDir,'*'));
files = files(~[files.isdir]);

if isempty(files)
    error('No files found in selected directory.');
end

fprintf('\nData directory:\n%s\n',dataDir);

%% =========================================================
% Select reference station
% =========================================================

[fileRef,pathRef] = uigetfile(fullfile(dataDir,'*'), ...
    'Select reference MiniSEED file');
if isequal(fileRef,0)
    return;
end

refFile = fullfile(pathRef,fileRef);
fprintf('\nReference station: %s\n',fileRef);

%% =========================================================
% Input window parameters
% =========================================================

answer = inputdlg({ ...
    'Number of repeated periods per window:', ...
    'Number of periods to slide each window:'}, ...
    'Goertzel window parameters',1,{'20','1'});

if isempty(answer)
    return;
end

repeatCycles = str2double(answer{1});
slideCycles  = str2double(answer{2});

if ~isfinite(repeatCycles) || repeatCycles <= 0 || ...
   ~isfinite(slideCycles)  || slideCycles <= 0
    error('Repeat cycles and slide cycles must be positive numbers.');
end

fprintf('\nWindow parameters:\n');
fprintf('  Repeat cycles = %.0f\n',repeatCycles);
fprintf('  Slide cycles  = %.0f\n',slideCycles);

%% =========================================================
% Read and filter reference station
% =========================================================

fprintf('\nReading reference station...\n');
ref = readBHZMiniSEED(refFile,targetChannel);

if isempty(ref.data)
    error('No valid BHZ data found in reference file.');
end

fprintf('  Sample rate = %.6f Hz\n',ref.fs);
fprintf('  Start       = %s\n',datestr(ref.startTime));
fprintf('  End         = %s\n',datestr(ref.endTime));
fprintf('  Samples     = %d\n',length(ref.data));

ref.data = bandpassFilter(ref.data,ref.fs,fLow,fHigh,filterOrder);

%% =========================================================
% Process every other station independently
% =========================================================

for k = 1:length(files)

    fileName = files(k).name;
    filePath = fullfile(dataDir,fileName);

    if strcmp(fileName,fileRef)
        continue;
    end

    fprintf('\n============================================================\n');
    fprintf('Target station: %s\n',fileName);
    fprintf('============================================================\n');

    %% ---------------------------------------------------------
    % Read target station
    % ----------------------------------------------------------

    try
        target = readBHZMiniSEED(filePath,targetChannel);
    catch ME
        fprintf('  ERROR reading file: %s\n',ME.message);
        continue;
    end

    if isempty(target.data)
        fprintf('  No valid BHZ data.\n');
        continue;
    end

    fprintf('  Sample rate = %.6f Hz\n',target.fs);
    fprintf('  Start       = %s\n',datestr(target.startTime));
    fprintf('  End         = %s\n',datestr(target.endTime));
    fprintf('  Samples     = %d\n',length(target.data));

    target.data = bandpassFilter(target.data,target.fs, ...
        fLow,fHigh,filterOrder);

    %% ---------------------------------------------------------
    % Determine common time for THIS pair only
    % ----------------------------------------------------------

    commonStart = max(ref.startTime,target.startTime);
    commonEnd   = min(ref.endTime,target.endTime);

    if commonEnd <= commonStart
        fprintf('  No common time interval.\n');
        continue;
    end

    fprintf('\n  Pair common interval:\n');
    fprintf('    Start = %s\n',datestr(commonStart));
    fprintf('    End   = %s\n',datestr(commonEnd));
    fprintf('    Duration = %.3f min\n',minutes(commonEnd-commonStart));

    % Only the common START is moved forward by 600 s.
    processStart = commonStart + seconds(shiftSeconds);
    processEnd   = commonEnd;

    if processEnd <= processStart
        fprintf('  Common interval is shorter than the 600 s stabilization shift.\n');
        continue;
    end

    fprintf('  Processing interval after +600 s:\n');
    fprintf('    Start = %s\n',datestr(processStart));
    fprintf('    End   = %s\n',datestr(processEnd));

    %% ---------------------------------------------------------
    % Align the two data sets to the pair start time
    % ----------------------------------------------------------

    [refUse,refStartActual] = extractTimeRange(ref,processStart,processEnd);
    [targetUse,targetStartActual] = extractTimeRange(target,processStart,processEnd);

    if isempty(refUse) || isempty(targetUse)
        fprintf('  Unable to extract pair data.\n');
        continue;
    end

    % The two sample rates should normally be identical. If not,
    % resample the target to the reference sample rate.
    if abs(target.fs-ref.fs) > 1e-6
        fprintf('  Resampling target: %.6f -> %.6f Hz\n', ...
            target.fs,ref.fs);

        targetUse = resample(targetUse,round(ref.fs),round(target.fs));
        fsUse = ref.fs;
    else
        fsUse = ref.fs;
    end

    N = min(length(refUse),length(targetUse));
    refUse = double(refUse(1:N));
    targetUse = double(targetUse(1:N));

    if N <= 0
        fprintf('  No usable paired samples.\n');
        continue;
    end

    fprintf('  Actual paired start: reference=%s, target=%s\n', ...
        datestr(refStartActual),datestr(targetStartActual));
    fprintf('  Paired samples = %d\n',N);

    %% ---------------------------------------------------------
    % Frequency loop
    % ----------------------------------------------------------

    for fi = 1:length(frequencyList)

        f = frequencyList(fi);
        period = 1/f;

        windowSeconds = repeatCycles * period;
        slideSeconds  = slideCycles  * period;

        windowSamples = round(windowSeconds*fsUse);
        slideSamples  = round(slideSeconds*fsUse);

        % Recalculate the actual window/slide lengths from integer
        % samples. The same sample indices are used for both stations.
        if windowSamples > N
            fprintf('\n  %.2f Hz: data shorter than one window.\n',f);
            continue;
        end

        if slideSamples < 1
            slideSamples = 1;
        end

        fprintf('\n------------------------------------------------------------\n');
        fprintf('Frequency = %.2f Hz\n',f);
        fprintf('Period    = %.6f s\n',period);
        fprintf('Window    = %.6f s (%d samples)\n', ...
            windowSamples/fsUse,windowSamples);
        fprintf('Slide     = %.6f s (%d samples)\n', ...
            slideSamples/fsUse,slideSamples);
        fprintf('------------------------------------------------------------\n');

        %% -----------------------------------------------------
        % Goertzel window loop
        % ------------------------------------------------------

        % First calculate and save ALL window complex cross-spectra.
        % Screening is deliberately performed only after the complete
        % frequency has been calculated, so isolated large values cannot
        % distort the statistics used for screening.
        GrefList = [];
        GtargetList = [];
        crossList = [];
        amplitudeList = [];
        windowNoList = [];

        startIndex = 1;
        windowNo = 0;

        while startIndex + windowSamples - 1 <= N

            i1 = startIndex;
            i2 = startIndex + windowSamples - 1;

            refWindow = refUse(i1:i2);
            targetWindow = targetUse(i1:i2);

            refWindow = refWindow - mean(refWindow);
            targetWindow = targetWindow - mean(targetWindow);

            Gref = goertzelComplex(refWindow,f,fsUse);
            Gtarget = goertzelComplex(targetWindow,f,fsUse);
            crossValue = Gtarget * conj(Gref);

            windowNo = windowNo + 1;
            GrefList(windowNo,1) = Gref;
            GtargetList(windowNo,1) = Gtarget;
            crossList(windowNo,1) = crossValue;
            amplitudeList(windowNo,1) = abs(crossValue);
            windowNoList(windowNo,1) = windowNo;

            startIndex = startIndex + slideSamples;
        end

        if isempty(crossList)
            fprintf('  No complete windows at %.2f Hz.\n',f);
            continue;
        end

        % Robust statistics for absolute cross-spectrum amplitude.
        amplitudeMedian = median(amplitudeList);
        amplitudeMAD = median(abs(amplitudeList - amplitudeMedian));
        amplitudeScale = 1.4826 * amplitudeMAD;
        KAmplitude = 5;

        % Robust statistics for window-to-window amplitude change.
        amplitudeDelta = [0; abs(diff(amplitudeList))];
        deltaMedian = median(amplitudeDelta(2:end));
        deltaMAD = median(abs(amplitudeDelta(2:end) - deltaMedian));
        deltaScale = 1.4826 * deltaMAD;
        KDeltaAmplitude = 5;

        if amplitudeScale > 0
            amplitudeThreshold = amplitudeMedian + KAmplitude*amplitudeScale;
        else
            amplitudeThreshold = KAmplitude*amplitudeMedian;
        end

        if deltaScale > 0
            deltaThreshold = deltaMedian + KDeltaAmplitude*deltaScale;
        else
            deltaThreshold = KDeltaAmplitude*deltaMedian;
        end

        fprintf('  Screening amplitude: median=%.4e, MAD=%.4e, threshold=%.4e\n', ...
            amplitudeMedian,amplitudeMAD,amplitudeThreshold);
        fprintf('  Screening delta:     median=%.4e, MAD=%.4e, threshold=%.4e\n', ...
            deltaMedian,deltaMAD,deltaThreshold);

        % First apply the two amplitude-based criteria only.
        % The surviving windows are used to obtain a preliminary phase.
        amplitudeValidList = false(windowNo,1);

        for i = 1:windowNo
            if i == 1
                deltaValue = 0;
            else
                deltaValue = abs(amplitudeList(i)-amplitudeList(i-1));
            end

            amplitudeValidList(i) = isValidCrossSpectrum( ...
                GrefList(i),GtargetList(i),crossList(i), ...
                amplitudeThreshold,deltaThreshold,deltaValue,0,inf);
        end

        % =============================================================
        % Split the complex cross-spectrum into two propagation branches.
        %
        % Phase = 0 deg is the boundary between the two propagation
        % directions.  Do NOT change the sign or conjugate the complex
        % values here.  The two branches are screened and accumulated
        % independently.
        % =============================================================

        crossPhase = angle(crossList) * 180/pi;

        positiveBranch = amplitudeValidList & (crossPhase >= 0);
        negativeBranch = amplitudeValidList & (crossPhase < 0);

        positiveCount = sum(positiveBranch);
        negativeCount = sum(negativeBranch);

        fprintf('  Phase branches after amplitude screening:\n');
        fprintf('    0 ... +180 deg : %d windows\n',positiveCount);
        fprintf('   -180 ... 0 deg  : %d windows\n',negativeCount);

        % Preliminary resultant is calculated independently for each
        % branch.  This prevents the two propagation directions from
        % cancelling each other before phase screening.
        preliminaryCrossPositive = complex(0,0);
        preliminaryCrossNegative = complex(0,0);

        if positiveCount > 0
            preliminaryCrossPositive = sum(crossList(positiveBranch));
        end

        if negativeCount > 0
            preliminaryCrossNegative = sum(crossList(negativeBranch));
        end

        phaseThreshold = 90;

        preliminaryPhasePositive = NaN;
        preliminaryPhaseNegative = NaN;

        if abs(preliminaryCrossPositive) > 0 && ...
                isfinite(real(preliminaryCrossPositive)) && ...
                isfinite(imag(preliminaryCrossPositive))
            preliminaryPhasePositive = angle(preliminaryCrossPositive);
        end

        if abs(preliminaryCrossNegative) > 0 && ...
                isfinite(real(preliminaryCrossNegative)) && ...
                isfinite(imag(preliminaryCrossNegative))
            preliminaryPhaseNegative = angle(preliminaryCrossNegative);
        end

        fprintf('  Positive branch preliminary phase = ');
        if isfinite(preliminaryPhasePositive)
            fprintf('%+.3f deg\n',preliminaryPhasePositive*180/pi);
        else
            fprintf('N/A\n');
        end

        fprintf('  Negative branch preliminary phase = ');
        if isfinite(preliminaryPhaseNegative)
            fprintf('%+.3f deg\n',preliminaryPhaseNegative*180/pi);
        else
            fprintf('N/A\n');
        end

        fprintf('  Phase screening threshold = %.1f deg\n',phaseThreshold);

        % Phase difference is calculated against the preliminary
        % resultant of the corresponding branch only.
        phaseDifference = NaN(windowNo,1);

        for i = 1:windowNo
            if ~amplitudeValidList(i)
                continue;
            end

            if crossPhase(i) >= 0
                if isfinite(preliminaryPhasePositive)
                    phaseDifference(i) = ...
                        angle(crossList(i) * exp(-1i*preliminaryPhasePositive)) ...
                        * 180/pi;
                end
            else
                if isfinite(preliminaryPhaseNegative)
                    phaseDifference(i) = ...
                        angle(crossList(i) * exp(-1i*preliminaryPhaseNegative)) ...
                        * 180/pi;
                end
            end
        end

        % =============================================================
        % Final screening
        % =============================================================

        deltaAmplitudeList = zeros(windowNo,1);

        if windowNo > 1
            deltaAmplitudeList(2:end) = abs( ...
                amplitudeList(2:end) - amplitudeList(1:end-1));
        end

        validList = false(windowNo,1);

        for i = 1:windowNo
            validList(i) = isValidCrossSpectrum( ...
                GrefList(i), ...
                GtargetList(i), ...
                crossList(i), ...
                amplitudeThreshold, ...
                deltaThreshold, ...
                deltaAmplitudeList(i), ...
                phaseDifference(i), ...
                phaseThreshold);
        end

        % The validity function treats NaN phase differences as invalid.
        % Accumulate the two propagation branches independently.
        totalCrossPositive = complex(0,0);
        totalCrossNegative = complex(0,0);

        accumulationCountPositive = 0;
        accumulationCountNegative = 0;
        rejectedCount = 0;

        countListPositive = [];
        energyListPositive = [];
        phaseListPositive = [];
        realListPositive = [];
        imagListPositive = [];

        countListNegative = [];
        energyListNegative = [];
        phaseListNegative = [];
        realListNegative = [];
        imagListNegative = [];

        for i = 1:windowNo

            deltaValue = deltaAmplitudeList(i);
            valid = validList(i);

            if valid
                if crossPhase(i) >= 0

                    totalCrossPositive = ...
                        totalCrossPositive + crossList(i);

                    accumulationCountPositive = ...
                        accumulationCountPositive + 1;

                    totalAmplitude = abs(totalCrossPositive);
                    totalPhase = angle(totalCrossPositive)*180/pi;

                    countListPositive(end+1,1) = ...
                        accumulationCountPositive;
                    energyListPositive(end+1,1) = totalAmplitude;
                    phaseListPositive(end+1,1) = totalPhase;
                    realListPositive(end+1,1) = real(totalCrossPositive);
                    imagListPositive(end+1,1) = imag(totalCrossPositive);

                    fprintf('  Window %5d  POSITIVE  valid  |Cross|=%.4e  ', ...
                        windowNoList(i),amplitudeList(i));
                    fprintf('TotalAmplitude=%.4e  Phase=%+.3f deg\n', ...
                        totalAmplitude,totalPhase);

                else

                    totalCrossNegative = ...
                        totalCrossNegative + crossList(i);

                    accumulationCountNegative = ...
                        accumulationCountNegative + 1;

                    totalAmplitude = abs(totalCrossNegative);
                    totalPhase = angle(totalCrossNegative)*180/pi;

                    countListNegative(end+1,1) = ...
                        accumulationCountNegative;
                    energyListNegative(end+1,1) = totalAmplitude;
                    phaseListNegative(end+1,1) = totalPhase;
                    realListNegative(end+1,1) = real(totalCrossNegative);
                    imagListNegative(end+1,1) = imag(totalCrossNegative);

                    fprintf('  Window %5d  NEGATIVE  valid  |Cross|=%.4e  ', ...
                        windowNoList(i),amplitudeList(i));
                    fprintf('TotalAmplitude=%.4e  Phase=%+.3f deg\n', ...
                        totalAmplitude,totalPhase);
                end

            else
                rejectedCount = rejectedCount + 1;

                if ~amplitudeValidList(i)
                    fprintf('  Window %5d  REJECTED  |Cross|=%.4e  Delta=%.4e\n', ...
                        windowNoList(i),amplitudeList(i),deltaValue);
                else
                    fprintf('  Window %5d  REJECTED  |Cross|=%.4e  Branch=%s  PhaseDiff=%+.3f deg\n', ...
                        windowNoList(i),amplitudeList(i), ...
                        ternaryBranch(crossPhase(i)),phaseDifference(i));
                end
            end
        end

        totalValidCount = accumulationCountPositive + accumulationCountNegative;

        %% ---------------------------------------------------------
        % Final result for this frequency / station pair
        % ----------------------------------------------------------

        if totalValidCount == 0
            fprintf('  No valid windows at %.2f Hz.\n',f);
            continue;
        end

        finalEnergyPositive = abs(totalCrossPositive);
        finalPhasePositive = angle(totalCrossPositive)*180/pi;

        finalEnergyNegative = abs(totalCrossNegative);
        finalPhaseNegative = angle(totalCrossNegative)*180/pi;

        fprintf('\nFinal result: %s relative to %s\n',fileName,fileRef);
        fprintf('  Frequency       = %.2f Hz\n',f);

        fprintf('  Positive branch (0 ... +180 deg):\n');
        fprintf('    Valid windows = %d\n',accumulationCountPositive);
        fprintf('    Total complex = %+.8e %+.8ei\n', ...
            real(totalCrossPositive),imag(totalCrossPositive));
        fprintf('    Total energy  = %.8e\n',finalEnergyPositive);
        fprintf('    Total phase   = %+.6f deg\n',finalPhasePositive);

        fprintf('  Negative branch (-180 ... 0 deg):\n');
        fprintf('    Valid windows = %d\n',accumulationCountNegative);
        fprintf('    Total complex = %+.8e %+.8ei\n', ...
            real(totalCrossNegative),imag(totalCrossNegative));
        fprintf('    Total energy  = %.8e\n',finalEnergyNegative);
        fprintf('    Total phase   = %+.6f deg\n',finalPhaseNegative);

        %% ---------------------------------------------------------
        % Plot accumulation curves
        % ----------------------------------------------------------

        figure('Color','w', ...
            'Name',sprintf('%s - %.2f Hz',fileName,f));

        subplot(2,1,1);
        hold on;

        if ~isempty(countListPositive)
            plot(countListPositive,energyListPositive,'LineWidth',1.5);
        end

        if ~isempty(countListNegative)
            plot(countListNegative,energyListNegative,'LineWidth',1.5);
        end

        grid on;
        xlabel('Accumulated valid windows');
        ylabel('|Total cross-spectrum|');
        title(sprintf('%s vs reference, %.2f Hz: accumulated energy', ...
            fileName,f),'Interpreter','none');

        legend('0 ... +180 deg','-180 ... 0 deg','Location','best');
        hold off;

        subplot(2,1,2);
        hold on;

        if ~isempty(countListPositive)
            plot(countListPositive,phaseListPositive,'LineWidth',1.5);
        end

        if ~isempty(countListNegative)
            plot(countListNegative,phaseListNegative,'LineWidth',1.5);
        end

        grid on;
        xlabel('Accumulated valid windows');
        ylabel('Phase (degree)');
        title(sprintf('Separate branch accumulation phases, %.2f Hz',f));

        legend('0 ... +180 deg','-180 ... 0 deg','Location','best');
        hold off;

        drawnow;

        %% ---------------------------------------------------------
        % Save result in memory for this pair/frequency
        % ----------------------------------------------------------

        result.frequency = f;
        result.repeatCycles = repeatCycles;
        result.slideCycles = slideCycles;
        result.windowSamples = windowSamples;
        result.slideSamples = slideSamples;

        result.validWindowsPositive = accumulationCountPositive;
        result.validWindowsNegative = accumulationCountNegative;
        result.validWindows = totalValidCount;

        result.totalCrossPositive = totalCrossPositive;
        result.totalCrossNegative = totalCrossNegative;

        result.finalEnergyPositive = finalEnergyPositive;
        result.finalEnergyNegative = finalEnergyNegative;

        result.finalPhasePositive = finalPhasePositive;
        result.finalPhaseNegative = finalPhaseNegative;

        result.cross = crossList;
        result.crossAmplitude = amplitudeList;
        result.crossPhase = crossPhase;

        result.branchPositive = positiveBranch;
        result.branchNegative = negativeBranch;

        result.valid = validList;
        result.rejectedWindows = rejectedCount;

        result.amplitudeMedian = amplitudeMedian;
        result.amplitudeMAD = amplitudeMAD;
        result.amplitudeThreshold = amplitudeThreshold;
        result.deltaMedian = deltaMedian;
        result.deltaMAD = deltaMAD;
        result.deltaThreshold = deltaThreshold;

        result.countPositive = countListPositive;
        result.energyPositive = energyListPositive;
        result.phasePositive = phaseListPositive;
        result.realPositive = realListPositive;
        result.imagPositive = imagListPositive;

        result.countNegative = countListNegative;
        result.energyNegative = energyListNegative;
        result.phaseNegative = phaseListNegative;
        result.realNegative = realListNegative;
        result.imagNegative = imagListNegative;

        result.preliminaryPhasePositive = preliminaryPhasePositive*180/pi;
        result.preliminaryPhaseNegative = preliminaryPhaseNegative*180/pi;

        result.phaseDifference = phaseDifference;
        result.phaseThreshold = phaseThreshold;

        assignin('base','lastGoertzelResult',result);

    end
end

fprintf('\n============================================================\n');
fprintf('Finished.\n');
fprintf('============================================================\n');

end


%% ========================================================================
% Read BHZ MiniSEED using the existing rdmseed structure
% ========================================================================
function out = readBHZMiniSEED(filePath,targetChannel)

mseed = rdmseed(filePath);

out = struct('data',[],'fs',[],'startTime',[], ...
             'endTime',[],'records',[]);

if isempty(mseed)
    return;
end

ch = {mseed.ChannelIdentifier};
idx = find(strcmp(ch,targetChannel));

if isempty(idx)
    return;
end

bhz = mseed(idx);

% Use the actual fields already used in the original program.
fs = bhz(1).SampleRate;

st = datetime(bhz(1).RecordStartTimeMATLAB,'ConvertFrom','datenum');
last = bhz(end);
et = datetime(last.RecordStartTimeMATLAB,'ConvertFrom','datenum') + ...
     seconds((last.NumberSamples-1)/last.SampleRate);

parts = {};
validRecords = [];

for r = 1:length(bhz)
    rec = bhz(r);

    recStart = datetime(rec.RecordStartTimeMATLAB,'ConvertFrom','datenum');
    recFs = rec.SampleRate;
    recData = double(rec.d(:));

    if isempty(recData)
        continue;
    end

    % Check that sample rate is consistent.
    if abs(recFs-fs) > 1e-6
        error('Sample rate changes inside %s.',filePath);
    end

    % Check continuity using record header time and number of samples.
    if ~isempty(validRecords)
        prev = bhz(validRecords(end));
        prevStart = datetime(prev.RecordStartTimeMATLAB,'ConvertFrom','datenum');
        prevEnd = prevStart + seconds((prev.NumberSamples-1)/prev.SampleRate);
        expectedStart = prevEnd + seconds(1/fs);
        gapSamples = seconds(recStart-expectedStart)*fs;

        if abs(gapSamples) > 1.0
            fprintf('  WARNING: gap/discontinuity in %s before record %d: %.3f samples\n', ...
                filePath,r,gapSamples);
        end
    end

    parts{end+1} = recData;
    validRecords(end+1) = r;
end

if isempty(parts)
    return;
end

out.data = vertcat(parts{:});
out.fs = fs;
out.startTime = st;
out.endTime = et;
out.records = bhz(validRecords);

end


%% ========================================================================
% Band-pass filter
% ========================================================================
function y = bandpassFilter(x,fs,fLow,fHigh,filterOrder)

x = double(x(:));
x = detrend(x);

if fHigh >= fs/2
    error('Upper filter frequency %.3f Hz is too high for fs=%.3f Hz.', ...
        fHigh,fs);
end

[b,a] = butter(filterOrder,[fLow fHigh]/(fs/2),'bandpass');
y = filtfilt(b,a,x);

end


%% ========================================================================
% Extract samples corresponding to a requested time range
% ========================================================================
function [x,tActual] = extractTimeRange(s,startTime,endTime)

fs = s.fs;

if endTime < s.startTime || startTime > s.endTime
    x = [];
    tActual = [];
    return;
end

startTime = max(startTime,s.startTime);
endTime   = min(endTime,s.endTime);

% Convert time difference to sample indices using the MiniSEED start time.
i1 = round(seconds(startTime-s.startTime)*fs)+1;
i2 = round(seconds(endTime-s.startTime)*fs)+1;

i1 = max(i1,1);
i2 = min(i2,length(s.data));

if i1 > i2
    x = [];
    tActual = [];
    return;
end

x = double(s.data(i1:i2));
tActual = s.startTime + seconds((i1-1)/fs);

end


%% ========================================================================
% Goertzel complex value at a specified frequency
% ========================================================================
function G = goertzelComplex(x,f,fs)

x = double(x(:));
N = length(x);

k = round(f*N/fs);
fk = k*fs/N;

% Direct Goertzel recurrence. The final complex value is equivalent
% to the DFT component at the selected Goertzel bin.
w = 2*pi*fk/fs;
coeff = 2*cos(w);

s1 = 0;
s2 = 0;

for n = 1:N
    s0 = x(n) + coeff*s1 - s2;
    s2 = s1;
    s1 = s0;
end

G = s1 - exp(-1i*w)*s2;

end



function name = ternaryBranch(phaseValue)

if phaseValue >= 0
    name = '0...+180';
else
    name = '-180...0';
end

end

%% ========================================================================
% Cross-spectrum validity judgment
% ========================================================================
function valid = isValidCrossSpectrum(Gref,Gtarget,crossValue, ...
    amplitudeThreshold,deltaThreshold,deltaValue,phaseDifference,phaseThreshold)

% Reject non-finite numerical results first.
if ~isfinite(real(Gref)) || ~isfinite(imag(Gref)) || ...
   ~isfinite(real(Gtarget)) || ~isfinite(imag(Gtarget)) || ...
   ~isfinite(real(crossValue)) || ~isfinite(imag(crossValue))
    valid = false;
    return;
end

% Criterion 1: absolute cross-spectrum amplitude is abnormally large.
if abs(crossValue) > amplitudeThreshold
    valid = false;
    return;
end

% Criterion 2: amplitude changes abnormally strongly from the previous
% window. The first window has deltaValue = 0 and is therefore retained.
if deltaValue > deltaThreshold
    valid = false;
    return;
end

% Criterion 3: the cross-spectrum must contribute in the same half-plane
% as the preliminary resultant. A phase difference larger than 90 degrees
% means that this window opposes the coherent accumulation.
if ~isfinite(phaseDifference)
    valid = false;
    return;
end

if abs(phaseDifference) > phaseThreshold
    valid = false;
    return;
end

valid = true;

end
