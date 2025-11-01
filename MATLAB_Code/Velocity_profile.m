%% MATLAB Script for Single-Frame 3D Position and Velocity Profile
% This code isolates a single frame and plots the detections in 3D (x, y, z) 
% with arrows (quivers) representing the radial velocity for clear visual
% inspection of directional consistency and potential interference.

clear all; close all; clc;

% --- CONFIGURATION ---
% 1. Input File: The file saved by the fixed-chunking script.
inputFileName = 'fixed_chunk_restructured_data.mat'; 

% 2. CRITICAL: The Frame ID we want to visualize.
FRAME_TO_ANALYZE = 1;  

% 3. Visualization Scale: Factor to visually scale the velocity arrows.
% Adjust this to make the arrows visible without overlapping.
VELOCITY_SCALE_FACTOR = 0.5; 

% 4. Max Range Display (m): Set the plot limits.
MAX_RANGE_DISPLAY = 5.0; 
% ---------------------

% NOTE: The variable name loaded will be 'fixedChunkFrameData'
try
    % 1. Load the frame-indexed data
    load(inputFileName, 'fixedChunkFrameData');
    
    if ~exist('fixedChunkFrameData', 'var')
        error('Variable ''fixedChunkFrameData'' not found. Check the input file.');
    end
    
    % 2. Check and Isolate the Specific Frame
    frameDataArray = fixedChunkFrameData;
    
    if FRAME_TO_ANALYZE > length(frameDataArray) || FRAME_TO_ANALYZE < 1
        error('Frame ID %d is invalid. Data contains only %d frames.', FRAME_TO_ANALYZE, length(frameDataArray));
    end
    
    frameData = frameDataArray(FRAME_TO_ANALYZE);
    detections = frameData.detections;
    
    if isempty(detections)
        disp(sprintf('Frame %d has no detected objects. Cannot generate 3D plot.', FRAME_TO_ANALYZE));
        return;
    end
    
    % 3. Extract and Prepare 3D Data
    % Q: Position vectors (x, y, z)
    Qx = double([detections.x]);
    Qy = double([detections.y]); % Forward direction
    Qz = double([detections.z]);
    
    % V: Velocity data (magnitudes)
    V = double([detections.velocity]);
    
    % R: Radial distances (for normalization)
    R = double([detections.range]); 
    
    % 4. Calculate 3D Velocity Components (Ux, Uy, Uz)
    % The radial velocity (V) is projected onto the x, y, z axes based on the 
    % normalized position vector (Qx/R, Qy/R, Qz/R).
    
    % Ensure no division by zero for points exactly at the origin (R=0)
    R(R == 0) = 1e-6; 
    
    % Calculate the vector components of velocity
    % Ux = V * (Qx / R) * SCALE
    % Uy = V * (Qy / R) * SCALE
    % Uz = V * (Qz / R) * SCALE
    Ux = V .* (Qx ./ R) * VELOCITY_SCALE_FACTOR;
    Uy = V .* (Qy ./ R) * VELOCITY_SCALE_FACTOR;
    Uz = V .* (Qz ./ R) * VELOCITY_SCALE_FACTOR;
    
    % --- 5. Generate 3D Plot with Velocity Quivers ---
    
    hFigure = figure('Name', sprintf('3D Position & Velocity Profile - Frame %d', FRAME_TO_ANALYZE), ...
                     'NumberTitle', 'off', 'Position', [100 100 900 900]);
    ax = axes('Parent', hFigure, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
    view(ax, 3); % Set 3D view
    hold(ax, 'on');
    
    % Plot the detection points
    scatter3(ax, Qx, Qy, Qz, 50, 'r', 'filled', 'DisplayName', 'Detection Point');
    
    % Plot the velocity vectors using quiver3
    hQuiver = quiver3(ax, Qx, Qy, Qz, Ux, Uy, Uz, 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.2, 'DisplayName', 'Radial Velocity Vector');
    
    % Add a sphere/marker for the radar origin
    scatter3(ax, 0, 0, 0, 100, 'k', 's', 'filled', 'DisplayName', 'Radar Origin');
    
    % Set plot limits and labels
    ax.XLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
    ax.YLim = [0, MAX_RANGE_DISPLAY]; % Y-axis is forward from radar
    ax.ZLim = [-MAX_RANGE_DISPLAY, MAX_RANGE_DISPLAY];
    
    xlabel(ax, 'X (meters) - Left/Right');
    ylabel(ax, 'Y (meters) - Forward');
    zlabel(ax, 'Z (meters) - Height');
    title(ax, sprintf('3D Velocity Profile: Frame ID %d', FRAME_TO_ANALYZE));
    legend(ax, 'Location', 'best');
    axis equal; % Ensure proper 3D scaling
    
    % Interference Check Guidance
    disp('*** VISUAL INTERFERENCE CHECK ***');
    disp('Normal movement (e.g., a walking target) should show velocity vectors that are:');
    disp('1. Clustered (arrows grouped together).');
    disp('2. Parallel (arrows pointing in the same general direction).');
    disp('Interference or noise often results in single points with:');
    disp('3. Random direction and position (arrows pointing haphazardly).');
    fprintf('Successfully generated 3D plot for visual interference check on Frame ID: %d.\n', FRAME_TO_ANALYZE);

catch ME
    disp(['An error occurred: ', ME.message]);
end
