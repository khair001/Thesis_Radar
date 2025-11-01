%% MATLAB Script to Chunk Aggregated Data into Fixed-Size Frames
% WARNING: This method assumes a CONSTANT number of objects per frame.
% It will FAIL to logically separate frames if the object count was variable.

clear; close all; clc;

% --- CONFIGURATION ---
% 1. Specify the path to your existing .mat file (the one with 'allLoggedObjects')
inputFileName = 'radar_logs/baselinelutterremoval500frame_test_data.mat'; 

% 2. CRITICAL: Set the fixed number of detections per assumed frame.
% We will use 5, but set this to whatever you believe is the 
% fixed number from your experiment configuration.
ASSUMED_OBJECTS_PER_FRAME = 5; 
% ---------------------

outputVariableName = 'fixedChunkFrameData';
outputFileName = 'fixed_chunk_restructured_data.mat';

try
    % 1. Load the aggregated data
    load(inputFileName, 'allLoggedObjects');
    
    if ~exist('allLoggedObjects', 'var')
        error('Variable ''allLoggedObjects'' not found in the loaded .mat file.');
    end
    
    totalObjects = length(allLoggedObjects);
    
    % 2. Validation Check and Calculation
    remainder = mod(totalObjects, ASSUMED_OBJECTS_PER_FRAME);
    if remainder ~= 0
        warning(['Total objects (%d) is not cleanly divisible by the assumed objects per frame (%d). ',...
                 'The last %d objects will be discarded to ensure full chunks.'], ...
                 totalObjects, ASSUMED_OBJECTS_PER_FRAME, remainder);
        % Truncate data to process only full frames
        allLoggedObjects = allLoggedObjects(1 : end - remainder);
        totalObjects = length(allLoggedObjects);
    end
    
    numFrames = totalObjects / ASSUMED_OBJECTS_PER_FRAME;
    fprintf('Loaded %d total objects. Chunking into %d frames with %d objects per frame.\n', ...
            totalObjects, numFrames, ASSUMED_OBJECTS_PER_FRAME);

    % 3. Initialize the new frame-based structure
    fixedChunkFrameData = struct('frameNumber', cell(1, numFrames), ...
                                 'numDetections', cell(1, numFrames), ...
                                 'detections', cell(1, numFrames));
    
    % 4. Iterate and chunk the data into frames
    for frameIdx = 1:numFrames
        % Calculate start and end indices for the current frame
        startIndex = (frameIdx - 1) * ASSUMED_OBJECTS_PER_FRAME + 1;
        endIndex = frameIdx * ASSUMED_OBJECTS_PER_FRAME;
        
        % Extract the slice of objects for the current frame
        frameDetections = allLoggedObjects(startIndex:endIndex);
        
        % Populate the new frame structure
        fixedChunkFrameData(frameIdx).frameNumber = frameIdx;
        fixedChunkFrameData(frameIdx).numDetections = length(frameDetections);
        fixedChunkFrameData(frameIdx).detections = frameDetections;
    end

    % 5. Save the restructured data
    save(outputFileName, outputVariableName);
    fprintf('\nSuccessfully chunked data into %d frames.\n', numFrames);
    fprintf('Saved new frame-indexed data to: %s (variable: %s)\n', outputFileName, outputVariableName);

catch ME
    disp(['An error occurred: ', ME.message]);
end
