filename = '431 JU298 on food R_2011_03_03__11_17_44___2___1.hdf5';
info = h5info(filename);

% get the stage motion data
stageData = h5read(filename, '/stage_data');

% for relatively small videos, we can just load everything into memory
images = h5read(filename, '/mask');

% set parameters
timeDiff = 1; % how many frames between aligned images?
dS = 5; % pixel downsampling factor (2 means half size)

% estimate transformation from one image frame to another
xShift = NaN(size(images, 3)-timeDiff, 1);
yShift = NaN(size(images, 3)-timeDiff, 1);
for ii = 1+timeDiff:size(images, 3)
    disp(ii/size(images, 3))
    
    % estimate shift between images
    transMat = imregcorr(images(1:dS:end, 1:dS:end, ii), ...
        images(1:dS:end, 1:dS:end, ii - timeDiff), 'translation');
    xShift(ii - timeDiff) = transMat.T(3, 1);
    yShift(ii - timeDiff) = transMat.T(3, 2);
    
end

figure; plot(xShift)
figure; plot(yShift)
figure; plot(abs(xShift) + abs(yShift))

% from the plot, we estimate that a total shift of greater than 0.5 is a
% peak
motionFrames = abs(xShift) + abs(yShift) > 0.5;
motionStarts = find(diff(motionFrames) == 1) + 1; % shift to account for diff
motionEnds = find(diff(motionFrames) == -1);

% check if the beginning or end frames are motion frames, if so, add starts
% and ends as appropriate
if motionFrames(1) == 1
    motionStarts = [1, motionStarts];
end
if motionFrames(end) == 1
    motionEnds = [motionEnds, length(motionFrames)];
end

% get the size of the x and y shifts
xShiftsIntegrated = NaN(numel(motionStarts), 1);
yShiftsIntegrated = NaN(numel(motionStarts), 1);

% loop through motion bouts and add up distances
for ii = 1:numel(motionStarts)
    xShiftsIntegrated(ii) = sum(xShift(motionStarts(ii):motionEnds(ii)));
    yShiftsIntegrated(ii) = sum(yShift(motionStarts(ii):motionEnds(ii)));
end

% correct shifts for downsampling, time delay, and pixel to micron
% calibration
pixel2micron = 20.9974/94.16126278325967;
xShiftsIntegrated = xShiftsIntegrated / pixel2micron * dS / timeDiff;
yShiftsIntegrated = yShiftsIntegrated / pixel2micron * dS / timeDiff;

% compare x and y shifts between video and log file to check for shifts
% visually. Note switch of x and y between log file and video.
figure
plot(1:length(diff(stageData.stage_x)), diff(stageData.stage_x), ...
    1:length(yShiftsIntegrated), yShiftsIntegrated)

figure
plot(1:length(diff(stageData.stage_y)), diff(stageData.stage_y), ...
    1:length(xShiftsIntegrated), xShiftsIntegrated)