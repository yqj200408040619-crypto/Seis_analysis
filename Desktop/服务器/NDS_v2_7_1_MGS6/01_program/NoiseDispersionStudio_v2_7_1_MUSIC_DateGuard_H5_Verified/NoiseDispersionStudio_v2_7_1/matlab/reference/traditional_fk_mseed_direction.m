function traditional_fk_mseed_source_direction()
% Traditional frequency-wavenumber (F-K) beamforming.
% This program is independent of the pair-phase algorithm.
%
% Coordinate input: longitude,latitude
% Direction: 0=East, 90=North, 180=West, 270=South.
% Returned direction is propagation direction.

clc;
close all;

%% Parameters
targetChannel = 'HHZ';
frequency = 0.3;
segmentSeconds = 1000;
frequencyHalfWidth = 0.005;

vMin = 1500;
vMax = 3500;
vStep = 100;

directionMin = 0;
directionMax = 360;
directionStep = 1;

nIntermediatePlots = 4;

%% Select directory
dataDir = uigetdir(pwd,'Select MiniSEED directory');
if isequal(dataDir,0)
    return;
end

files = dir(fullfile(dataDir,'*'));
files = files(~[files.isdir]);
files = files(isMiniSEEDFile(files));

if isempty(files)
    error('No MiniSEED files found in selected directory.');
end

fprintf('\n============================================================\n');
fprintf('Traditional F-K source-direction estimation\n');
fprintf('============================================================\n');
fprintf('Directory: %s\n',dataDir);
fprintf('MiniSEED files: %d\n',length(files));

%% Enter coordinates
nStation = length(files);
prompt = cell(nStation,1);
defaults = cell(nStation,1);

for i = 1:nStation
    prompt{i} = sprintf('%s  Longitude,Latitude (deg):',files(i).name);
    defaults{i} = '0,0';
end

answer = inputdlg(prompt,'Station coordinates: longitude,latitude',1,defaults);
if isempty(answer)
    return;
end

station = repmat(struct( ...
    'name','', ...
    'fileName','', ...
    'filePath','', ...
    'longitude',NaN, ...
    'latitude',NaN, ...
    'data',[], ...
    'fs',NaN, ...
    'startTime',[], ...
    'endTime',[]),nStation,1);

for i = 1:nStation
    parts = split(strtrim(answer{i}),',');
    if numel(parts) ~= 2
        error('Invalid coordinate for %s. Use longitude,latitude.',files(i).name);
    end

    lon = str2double(strtrim(parts{1}));
    lat = str2double(strtrim(parts{2}));

    if ~isfinite(lon) || ~isfinite(lat) || lon < -180 || lon > 180 || lat < -90 || lat > 90
        error('Invalid coordinate for %s.',files(i).name);
    end

    [~,name,~] = fileparts(files(i).name);
    station(i).name = name;
    station(i).fileName = files(i).name;
    station(i).filePath = fullfile(files(i).folder,files(i).name);
    station(i).longitude = lon;
    station(i).latitude = lat;
end

%% Read MiniSEED
for i = 1:nStation
    fprintf('\nReading %s (%d/%d)\n',station(i).fileName,i,nStation);

    station(i).data = readChannelMiniSEED(station(i).filePath,targetChannel);

    if isempty(station(i).data.data)
        error('No %s data found in %s.',targetChannel,station(i).fileName);
    end

    station(i).fs = station(i).data.fs;
    station(i).startTime = station(i).data.startTime;
    station(i).endTime = station(i).data.endTime;

    fprintf('  fs    = %.6f Hz\n',station(i).fs);
    fprintf('  start = %s\n',datestr(station(i).startTime));
    fprintf('  end   = %s\n',datestr(station(i).endTime));
    fprintf('  N     = %d\n',length(station(i).data.data));
end

fs = station(1).fs;
for i = 2:nStation
    if abs(station(i).fs-fs) > 1e-6
        error('Sample rates are not identical.');
    end
end

%% Common interval
commonStart = max([station.startTime]);
commonEnd = min([station.endTime]);

if commonEnd <= commonStart
    error('No common time interval exists among the stations.');
end

fprintf('\nCommon interval:\n');
fprintf('  %s -> %s\n',datestr(commonStart),datestr(commonEnd));

if seconds(commonEnd-commonStart) < segmentSeconds
    error('Common interval is shorter than segmentSeconds.');
end

%% Local East/North coordinates
R = 6371000;
lat0 = deg2rad(mean([station.latitude]));
lon0 = deg2rad(mean([station.longitude]));

lat = deg2rad([station.latitude]);
lon = deg2rad([station.longitude]);

xEast = R*cos(lat0).*(lon-lon0);
yNorth = R*(lat-lat0);

fprintf('\nStation coordinates:\n');
for i = 1:nStation
    fprintf('  %-20s E = %10.1f m  N = %10.1f m\n', ...
        station(i).name,xEast(i),yNorth(i));
end

%% Extract common data segment
segmentSamples = round(segmentSeconds*fs);
data = zeros(nStation,segmentSamples);

for i = 1:nStation
    [x,tActual] = extractTimeRange(station(i).data, ...
        commonStart,commonStart+seconds(segmentSeconds));

    if length(x) < segmentSamples
        error('Insufficient data at station %s.',station(i).name);
    end

    data(i,:) = double(x(1:segmentSamples));

    if i == 1
        actualStart = tActual;
    end
end

fprintf('\nAnalysis segment:\n');
fprintf('  start = %s\n',datestr(actualStart));
fprintf('  duration = %.1f s\n',segmentSeconds);

%% Normalize station amplitudes
for i = 1:nStation
    data(i,:) = data(i,:) - mean(data(i,:));
    rmsValue = sqrt(mean(data(i,:).^2));
    if rmsValue > 0
        data(i,:) = data(i,:)/rmsValue;
    end
end

%% FFT
N = segmentSamples;
window = hann(N).';
X = zeros(nStation,N);

for i = 1:nStation
    X(i,:) = fft(data(i,:).*window);
end

freqAxis = (0:N-1)*fs/N;
frequencyMask = abs(freqAxis-frequency) <= frequencyHalfWidth;
frequencyMask = frequencyMask & freqAxis <= fs/2;
frequencyBins = find(frequencyMask);

if isempty(frequencyBins)
    error('No FFT bin exists near %.6f Hz.',frequency);
end

fprintf('\nFFT bins: %.6f -> %.6f Hz\n', ...
    freqAxis(frequencyBins(1)),freqAxis(frequencyBins(end)));

Xf = mean(X(:,frequencyBins),2);

%% Array geometry
figure('Color','w','Name','Array geometry');
plot(xEast/1000,yNorth/1000,'ko','MarkerFaceColor','k');
hold on;
for i = 1:nStation
    text(xEast(i)/1000,yNorth(i)/1000,['  ' station(i).name], ...
        'Interpreter','none');
end
axis equal;
grid on;
xlabel('East (km)');
ylabel('North (km)');
title('Station geometry');
quiver(0,0,5,0,0,'LineWidth',1);
text(5.2,0,'E');
quiver(0,0,0,5,0,'LineWidth',1);
text(0,5.2,'N');
hold off;

%% Traditional F-K scan
velocityList = vMin:vStep:vMax;
directionList = directionMin:directionStep:directionMax-directionStep;
powerMatrix = zeros(length(velocityList),length(directionList));

snapshotVelocityIndex = unique(round(linspace(1,length(velocityList), ...
    nIntermediatePlots)));

figure('Color','w','Name','F-K fitting process');
snapshotNumber = 0;

for iv = 1:length(velocityList)
    v = velocityList(iv);

    for itheta = 1:length(directionList)
        theta = deg2rad(directionList(itheta));
        projected = xEast*cos(theta) + yNorth*sin(theta);

        % Conventional plane-wave steering.
        steering = exp(-1i*2*pi*frequency*projected/v);
        beam = sum(Xf(:).*steering(:));
        powerMatrix(iv,itheta) = abs(beam)^2;
    end

    if any(iv == snapshotVelocityIndex)
        snapshotNumber = snapshotNumber + 1;

        subplot(2,2,snapshotNumber);
        imagesc(directionList,velocityList(1:iv),powerMatrix(1:iv,:));
        axis xy;
        xlabel('Propagation direction (degree)');
        ylabel('Phase velocity (m/s)');
        title(sprintf('F-K search up to %.0f m/s',v));
        colorbar;
        hold on;

        [~,idx] = max(powerMatrix(1:iv,:),[],'all','linear');
        [bestIv,bestIt] = ind2sub(size(powerMatrix(1:iv,:)),idx);

        plot(directionList(bestIt),velocityList(bestIv), ...
            'kx','MarkerSize',12,'LineWidth',2);
        hold off;
        drawnow;
    end
end

%% Final result
[maxPower,linearIndex] = max(powerMatrix(:));
[bestVelocityIndex,bestDirectionIndex] = ...
    ind2sub(size(powerMatrix),linearIndex);

bestVelocity = velocityList(bestVelocityIndex);
bestDirection = directionList(bestDirectionIndex);

fprintf('\n============================================================\n');
fprintf('Traditional F-K result\n');
fprintf('============================================================\n');
fprintf('Frequency             = %.6f Hz\n',frequency);
fprintf('Propagation direction = %.3f deg\n',bestDirection);
fprintf('Phase velocity        = %.3f m/s\n',bestVelocity);
fprintf('Maximum beam power    = %.6e\n',maxPower);

%% Final power map
figure('Color','w','Name','Traditional F-K result');
imagesc(directionList,velocityList,powerMatrix);
axis xy;
xlabel('Propagation direction (degree)');
ylabel('Phase velocity (m/s)');
title(sprintf('Traditional F-K beam power at %.3f Hz',frequency));
colorbar;
hold on;
plot(bestDirection,bestVelocity,'kx','MarkerSize',14,'LineWidth',2);
text(bestDirection,bestVelocity, ...
    sprintf('  %.1f deg, %.0f m/s',bestDirection,bestVelocity), ...
    'VerticalAlignment','bottom');
hold off;

%% Direction spectrum
figure('Color','w','Name','Direction spectrum');
plot(directionList,powerMatrix(bestVelocityIndex,:),'LineWidth',1.2);
grid on;
xlabel('Propagation direction (degree)');
ylabel('Beam power');
title(sprintf('Direction spectrum at %.1f m/s',bestVelocity));

%% Save
result = struct();
result.frequency = frequency;
result.bestDirection = bestDirection;
result.bestVelocity = bestVelocity;
result.maxPower = maxPower;
result.directionList = directionList;
result.velocityList = velocityList;
result.powerMatrix = powerMatrix;
result.xEast = xEast;
result.yNorth = yNorth;
result.station = station;
result.commonStart = commonStart;
result.commonEnd = commonEnd;
result.analysisStart = actualStart;
result.segmentSeconds = segmentSeconds;
result.fs = fs;

outputFile = fullfile(dataDir,'traditional_fk_result.mat');
save(outputFile,'result','-v7.3');

fprintf('\nResult saved to:\n%s\n',outputFile);
end

function tf = isMiniSEEDFile(files)
n = length(files);
tf = false(n,1);
for i = 1:n
    [~,~,ext] = fileparts(files(i).name);
    ext = lower(ext);
    tf(i) = isempty(ext) || any(strcmp(ext,{'.mseed','.miniseed','.msd','.seed'}));
end
end

function out = readChannelMiniSEED(filePath,targetChannel)
mseed = rdmseed(filePath);
out = struct('data',[],'fs',[],'startTime',[],'endTime',[]);
if isempty(mseed)
    return;
end

ch = {mseed.ChannelIdentifier};
idx = find(strcmp(ch,targetChannel));
if isempty(idx)
    return;
end

rec = mseed(idx);
recordTime = [rec.RecordStartTimeMATLAB];
[~,order] = sort(recordTime);
rec = rec(order);

fs = rec(1).SampleRate;
parts = {};

for r = 1:length(rec)
    if abs(rec(r).SampleRate-fs) > 1e-6
        error('Sample rate changes inside %s.',filePath);
    end

    x = double(rec(r).d(:));
    if isempty(x)
        continue;
    end

    if ~isempty(parts)
        prev = rec(r-1);
        prevStart = datetime(prev.RecordStartTimeMATLAB,'ConvertFrom','datenum');
        prevEnd = prevStart + seconds((length(prev.d)-1)/fs);
        recStart = datetime(rec(r).RecordStartTimeMATLAB,'ConvertFrom','datenum');
        gapSamples = seconds(recStart-(prevEnd+seconds(1/fs)))*fs;

        if gapSamples > 1
            error('Data discontinuity in %s: %.3f samples.',filePath,gapSamples);
        elseif gapSamples < -1
            overlapSamples = round(-gapSamples);
            if overlapSamples >= length(x)
                continue;
            end
            x = x(overlapSamples+1:end);
        end
    end

    parts{end+1} = x;
end

if isempty(parts)
    return;
end

out.data = vertcat(parts{:});
out.fs = fs;
out.startTime = datetime(rec(1).RecordStartTimeMATLAB,'ConvertFrom','datenum');

last = rec(end);
lastStart = datetime(last.RecordStartTimeMATLAB,'ConvertFrom','datenum');
out.endTime = lastStart + seconds((length(last.d)-1)/fs);
end

function [x,tActual] = extractTimeRange(s,startTime,endTime)
if endTime < s.startTime || startTime > s.endTime
    x = [];
    tActual = [];
    return;
end

startTime = max(startTime,s.startTime);
endTime = min(endTime,s.endTime);
fs = s.fs;

i1 = round(seconds(startTime-s.startTime)*fs)+1;
i2 = round(seconds(endTime-s.startTime)*fs)+1;
i1 = max(i1,1);
i2 = min(i2,length(s.data));

if i1 > i2
    x = [];
    tActual = [];
    return;
end

x = s.data(i1:i2);
tActual = s.startTime + seconds((i1-1)/fs);
end
