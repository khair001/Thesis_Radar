clear all; close all; clc;

%% --- CONFIGURATION ---
BASELINE_FILES = {
    'radar_logs/baseline500frame_test_data.mat', 
    'radar_logs/baseline750frame_test_data.mat',
    'radar_logs/baseline1000frame_test_data.mat',
    'radar_logs/baseline2000frame_test_data.mat' 
};

BASELINE_VAR_NAME = 'allLoggedObjects'; 

% --- TEMPORARY DBSCAN PARAMETERS ---
TEMP_DBSCAN_EPSILON = 0.2; % Initial guess for cluster radius (m)
TEMP_DBSCAN_MINPTS = 5;    % Initial guess for cluster density

% --- STATISTICAL SETTINGS ---
JITTER_PERCENTILE = 99.7; % 99.7th percentile (approx. 3-sigma)

% --- OUTPUT VARIABLES ---
all_spatial_jitters = []; % Stores positional variation (in meters)
all_velocity_jitters = [];% Stores absolute velocity (in m/s)
all_static_centroids = [];% Stores centroids of static objects

%% --- STEP 1: LOAD REFERENCE BASELINE AND CREATE STATIC MAP ---
disp('--- Processing Reference Baseline Data ---');

if isempty(BASELINE_FILES)
    error('BASELINE_FILES list is empty. Please list clean data files.');
end

currentFile = BASELINE_FILES{1};

% --- DEBUGGING START (Confirmation that file loading is correct) ---
disp('--- DEBUGGING START ---');
disp(['Attempting to access file: ', currentFile]);
disp(['File existence (2=yes, 0=no): ', num2str(exist(currentFile, 'file'))]);
if exist(currentFile, 'file') ~= 2
    error(['File not found. Please check the path: ', currentFile]);
end
disp('Variables found in file:');
whos('-file', currentFile); 
disp('--- DEBUGGING END ---');

% Load the file (Reference)
refData = load(currentFile);

% Check and extract the object data
if isfield(refData, BASELINE_VAR_NAME)
    refObjects = refData.(BASELINE_VAR_NAME);
else
    error(['Variable ''', BASELINE_VAR_NAME, ''' not found. Fatal error.']);
end

% Cluster the reference data
[refClusters, ~] = clusterData(refObjects, TEMP_DBSCAN_EPSILON, TEMP_DBSCAN_MINPTS);

% Identify the STATIC objects in the reference map
TEMP_VELOCITY_THRESHOLD = 0.1; % Use a conservative threshold for mapping
refClusters = simpleClassification(refClusters, TEMP_VELOCITY_THRESHOLD);

% Extract centroids for all STATIC objects (our ground truth map)
refStaticClusters = refClusters([refClusters.isStatic]);
refStaticCentroids = vertcat(refStaticClusters.centroid);

if isempty(refStaticCentroids)
    warning('No static objects found in the reference baseline. Check TEMP_DBSCAN_EPSILON/MINPTS or the data.');
    return;
end
disp(['Successfully identified ', num2str(length(refStaticCentroids)), ' static objects in the reference map.']);

%% --- STEP 2: ITERATE AND COMPARE SUBSEQUENT BASELINE RUNS (JITTER CALCULATION) ---
disp(newline);
disp('--- Step 2: Calculating Spatial and Velocity Jitter ---');

for fileIndex = 2:length(BASELINE_FILES)
    currentFile = BASELINE_FILES{fileIndex};
    fprintf('  -> Analyzing file: %s\n', currentFile);
    
    try
        currentData = load(currentFile);
        currentObjects = currentData.(BASELINE_VAR_NAME);
    catch ME
        warning('Failed to load file %s. Skipping. Error: %s', currentFile, ME.message);
        continue;
    end
    
    % Cluster the current baseline data
    [currentClusters, ~] = clusterData(currentObjects, TEMP_DBSCAN_EPSILON, TEMP_DBSCAN_MINPTS);
    
    % --- Match Current Clusters to Reference Static Centroids ---
    
    for i = 1:length(currentClusters)
        currentCentroid = currentClusters(i).centroid;
        currentAvgVel = currentClusters(i).avgVelocity;
        
        % Calculate distance to all static reference centroids
        distances = sqrt(sum((refStaticCentroids - currentCentroid).^2, 2));
        
        % Find the closest static reference object
        [minDist, ~] = min(distances);
        
        % If the current cluster is spatially close to a known static object
        if minDist < 1.0 % TEMPORARY WIDE MATCHING THRESHOLD (1 meter)
            
            % 1. Record Spatial Jitter 
            all_spatial_jitters = [all_spatial_jitters; minDist];
            
            % 2. Record Velocity Jitter 
            all_velocity_jitters = [all_velocity_jitters; currentAvgVel];

            all_static_centroids = [all_static_centroids; currentCentroid];
        end
    end
end

if isempty(all_spatial_jitters)
    disp('No common static objects were consistently matched across the baseline runs.');
    return;
end

%% --- STEP 3: STATISTICAL ANALYSIS AND RECOMMENDATION ---
disp(newline);
disp('--- Step 3: Statistical Analysis and Final Recommendation ---');

% --- SPATIAL THRESHOLD CALCULATION ---
spatial_max_jitter = max(all_spatial_jitters);
spatial_percentile_value = prctile(all_spatial_jitters, JITTER_PERCENTILE);

% Round up to the next 0.1 meter
recommended_match_threshold = ceil(spatial_percentile_value * 10) / 10; 

fprintf('Total static object position measurements collected: %d\n', length(all_spatial_jitters));
fprintf('Max observed positional jitter: %.3f meters\n', spatial_max_jitter);
fprintf('%gth Percentile Positional Jitter: %.3f meters\n', JITTER_PERCENTILE, spatial_percentile_value);
fprintf('=================================================================\n');
fprintf('  --> Recommended BASELINE_MATCH_THRESHOLD: **%.1f** meters\n', recommended_match_threshold);
fprintf('=================================================================\n');

% --- VELOCITY THRESHOLD CALCULATION ---
velocity_max_jitter = max(all_velocity_jitters);
velocity_percentile_value = prctile(all_velocity_jitters, JITTER_PERCENTILE);

% Round up to the next 0.05 m/s
recommended_velocity_threshold = ceil(velocity_percentile_value * 20) / 20; 

fprintf('\nTotal static object velocity measurements collected: %d\n', length(all_velocity_jitters));
fprintf('Max observed absolute velocity (noise): %.3f m/s\n', velocity_max_jitter);
fprintf('%gth Percentile Velocity Jitter: %.3f m/s\n', JITTER_PERCENTILE, velocity_percentile_value);
fprintf('=================================================================\n');
fprintf('  --> Recommended STATIC_VELOCITY_THRESHOLD: **%.2f** m/s\n', recommended_velocity_threshold);
fprintf('=================================================================\n');

disp(newline);
disp('NOTE: DBSCAN_EPSILON should be slightly smaller than BASELINE_MATCH_THRESHOLD for effective clustering.');
disp(['RECOMMENDED DBSCAN_EPSILON: **', num2str(recommended_match_threshold - 0.1, '%.1f'), '** meters (or less).']);


%% --- STEP 4: VISUALIZATION (HISTOGRAMS) ---
figure('Name', 'Threshold Calibration Jitter Distribution', 'NumberTitle', 'off', 'Position', [100 100 1200 600]);

% 1. Spatial Jitter Histogram
subplot(1,2,1);
histogram(all_spatial_jitters, 'BinMethod', 'auto', 'FaceColor', [0 0.4470 0.7410]);
xline(spatial_percentile_value, 'r--', 'LineWidth', 2, 'DisplayName', '99.7th Percentile');
xline(recommended_match_threshold, 'g-', 'LineWidth', 2, 'DisplayName', 'Recommended Threshold');
title('Distribution of Static Object Positional Jitter');
xlabel('Positional Jitter (Distance from Reference Centroid) [m]');
ylabel('Frequency');
legend('show');
grid on;

% 2. Velocity Jitter Histogram
subplot(1,2,2);
histogram(all_velocity_jitters, 'BinMethod', 'auto', 'FaceColor', [0.8500 0.3250 0.0980]);
xline(velocity_percentile_value, 'r--', 'LineWidth', 2, 'DisplayName', '99.7th Percentile');
xline(recommended_velocity_threshold, 'g-', 'LineWidth', 2, 'DisplayName', 'Recommended Threshold');
title('Distribution of Static Object Velocity Jitter (Noise)');
xlabel('Absolute Average Velocity of Static Clusters [m/s]');
ylabel('Frequency');
legend('show');
grid on;

%% =========================================================================
%   LOCAL HELPER FUNCTIONS
% =========================================================================

function [clusteredObjects, outlierIdx] = clusterData(allDetections, epsilon, minPts)
    if isempty(allDetections)
        clusteredObjects = [];
        outlierIdx = [];
        return;
    end
    points = [[allDetections.x]', [allDetections.y]', [allDetections.z]'];
    clusterIdx = customDBSCAN(points, epsilon, minPts);
    outlierIdx = (clusterIdx == -1);
    numClusters = max(clusterIdx);
    
    clusteredObjects = repmat(struct('detections', [], 'avgVelocity', 0, 'centroid', [0,0,0], 'isStatic', false, 'isMoving', false), 1, numClusters);
    
    for i = 1:numClusters
        clusterDetections = allDetections(clusterIdx == i);
        avgVel = mean(abs([clusterDetections.velocity]));
        avg_x = mean([clusterDetections.x]); 
        avg_y = mean([clusterDetections.y]); 
        avg_z = mean([clusterDetections.z]);
        centroid = [avg_x, avg_y, avg_z];
        clusteredObjects(i).detections = clusterDetections;
        clusteredObjects(i).avgVelocity = avgVel;
        clusteredObjects(i).centroid = centroid;
    end
end

function clusteredObjects = simpleClassification(clusteredObjects, velocityThreshold)
    for i = 1:length(clusteredObjects)
        isMoving = clusteredObjects(i).avgVelocity > velocityThreshold;
        clusteredObjects(i).isMoving = isMoving;
        clusteredObjects(i).isStatic = ~isMoving;
    end
end

function clusterIdx = customDBSCAN(points, epsilon, minPts)
    numPoints = size(points, 1);
    clusterIdx = zeros(numPoints, 1);
    clusterCounter = 0;
    
    % --- Distance Calculation (Toolbox-independent) ---
    Z = zeros(numPoints, numPoints);
    for i = 1:numPoints
        for j = i+1:numPoints
            dist = sqrt(sum((points(i,:) - points(j,:)).^2));
            Z(i,j) = dist;
            Z(j,i) = dist;
        end
    end
    % --- End Distance Calculation ---

    for i = 1:numPoints
        if clusterIdx(i) == 0
            neighbors = find(Z(i, :) <= epsilon); 
            
            if length(neighbors) < minPts
                clusterIdx(i) = -1; % Noise point
            else
                clusterCounter = clusterCounter + 1;
                clusterIdx(i) = clusterCounter;
                k = 1;
                while k <= length(neighbors)
                    neighborIdx = neighbors(k);
                    
                    if clusterIdx(neighborIdx) == 0 % Unvisited
                        clusterIdx(neighborIdx) = clusterCounter;
                        
                        newNeighbors = find(Z(neighborIdx, :) <= epsilon);
                        
                        if length(newNeighbors) >= minPts
                            neighbors = [neighbors, setdiff(newNeighbors, neighbors)];
                        end
                    elseif clusterIdx(neighborIdx) == -1 % Previously marked as Noise
                        clusterIdx(neighborIdx) = clusterCounter;
                    end
                    k = k + 1;
                end
            end
        end
    end
end
