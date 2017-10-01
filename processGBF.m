% PROCESS_GBF
% 
% Name:       process_GBF.m
% Purpose:    Process the output files of the GBF model into probabilistic
%             hazard maps for ballistic impacts
% Author:     Sebastien Biass, Jean-Luc Falcone, Costanza Bonadonna
% Created:    April 2015
% Updated:    July 2015
% Copyright:  S Biass, JL Falcone, C Bonadonna - University of Geneva, 2015
% License:    GNU GPL3
% 
%         "You shake my nerves and you rattle my brain
%         Too much love drives a man insane
%         You broke my will, oh what a thrill
%         Goodness gracious great balls of fire"
%                                     -- J.L. Lewis
% 
% This is a free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     It is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with it. If not, see <http://www.gnu.org/licenses/>.

function project = processGBF(varargin)
addpath('dependencies/');

%{ 
Debug
    fl = 'vulcano_blocks.dat';
    inBal = struct('name', 'test',...
        'res', 100,...             % Resolution
        'eT', [60,100,4000,8000],...     % Vector of energy thresholds
        'pT', [10, 25, 50, 75, 90],...       % Vector of probability thresholds
        'dI', 250,...                      % Distance interval
        'rI', 20,...                    % Radial interval
        'vE', 496670,...                      % Vent easting
        'vN', 4250690,...                    % Vent easting
        'vZ', 33);                         % Vent UTM zone 
%}


%% Retrieve input parameters
% Case using the GUI
if nargin == 0
    
    [FileName, PathName] = uigetfile('*.*'); % Choose file
    if FileName == 0; return; end   % Check non-null value
    
    fl2load     = [PathName, filesep, FileName];    % Define file path
	input_val   = {'test_run','50', '100, 1000', '10, 25, 50, 75, 90', '1000', '22.5', '496670', '4250690', '33'}; % Default GUI values

    % GUI
    answer      = inputdlg({'Name:','Grid resolution (m):','Energy thresholds (J), comma delimited', 'Probability threshold (%), comma delimited', 'Distance interval (m)', 'Radial sector interval (degrees)', 'Vent easting (m)', 'Vent northing (m)', 'UTM Zone'},...
        'Input', [1 35; 1 35; 1 35; 1 35; 1 35; 1 35; 1 35; 1 35; 1 35],...
        input_val);
    
    % Check non-null value
    if isempty(answer); return; end
    
    % Fillup the input structure
    inBal = struct('name', answer{1},...
        'gridRes', str2double(answer{2}),...                % Resolution
        'eT', str2double(strsplit(answer{3}, ',')),...      % Vector of energy thresholds
        'pT', str2double(strsplit(answer{4}, ',')),...      % Vector of probability thresholds
        'dI', str2double(answer{5}),...                     % Distance interval
        'rI', str2double(answer{6}),...                     % Radial interval
        'vE', str2double(answer{7}),...                     % Vent easting
        'vN', str2double(answer{8}),...                     % Vent easting
        'vZ', str2double(answer{9}));                       % Vent UTM zone
        
% Case using a structure    
elseif nargin == 2
    % Case second argument is not a structure
    if ~isstruct(varargin{2})
        error('The second input argument should be a structure')
    end
    
    fl2load = varargin{1};
    inBal   = varargin{2};
    
% Case entering all input parameters separately    
elseif nargin == 10
    fl2load = varargin{1};
    inBal = struct('name', varargin{2},...
        'gridRes', varargin{3},...                          % Resolution
        'eT', varargin{4},...                               % Vector of energy thresholds
        'pT', varargin{5},...                               % Vector of probability thresholds
        'dI', varargin{6},...                               % Distance interval
        'rI', varargin{7},...                               % Radial interval
        'vE', varargin{8},...                               % Vent easting
        'vN', varargin{9},...                               % Vent easting
        'vZ', varargin{10});                                % Vent UTM zone
% Else error
else
    error('Wrong number of input arguments')
end

% New definition of the grid resolution by a 2x1 vector representing the
% resolution of the pixel-based approach and the resolution of the
% zone-based approach

%% Load file
% Create directory
if exist(inBal.name, 'dir') == 7
    choice = questdlg('The output name already exists, overwrite?', ...
	'Output name', ...
	'No','Yes','No');
    switch choice
        case 'No'
            return
        case 'Yes'
            rmdir(inBal.name,'s');
    end
end
mkdir(inBal.name);
pthout  = [inBal.name, filesep];         % Set output folder

display(sprintf('\tLoading file'));
data    = dlmread(fl2load, '', 1, 0);

if size(data,2) < 6
    error('The input file should have at least 6 columns organized as easting (m), northing (m), altitude (m asl), mass (kg), diameter (m) and kinetic energy (kJ)') 
end

% Get results from file
x       = data(:,1);                    % Easting
y       = data(:,2);                    % Northing
e       = data(:,6).*1000;              % Energy (J)
n       = size(data,1);                 % Number of bombs
d       = sqrt((x-inBal.vE).^2 + (y-inBal.vN).^2);  % Distance between vent and bombs
id      = zeros(size(x,1),2);           % Vectors of indices for cartesian grid

data    = data(:,3:end);                % Colums:
                                        % 1: Landing altitude (m a.s.l.)
                                        % 2: Mass (kg)
                                        % 3: Diameter (m)
                                        % 4: Kinetic energy (kJ)
                                        % 5: Landing angle (deg)
                                        % 6: Ejection andgle(deg)
                                        % 7: Flight time (sec)

%% Define the cartesian grid                                                                               
% Create cartesian grid
[gridEast, gridNorth]   = meshgrid(min(x)-inBal.gridRes(1):inBal.gridRes(1):max(x)+inBal.gridRes(1),...
    min(y)-inBal.gridRes(1):inBal.gridRes(1):max(y)+inBal.gridRes(1));
gridNorth               = flipud(gridNorth);

[gridLat,gridLon]       = utm2ll(gridEast,gridNorth,ones(size(gridEast)).*inBal.vZ);
[vLat,vLon]             = utm2ll(inBal.vE, inBal.vN, inBal.vZ);


%% 2 Conversion to cartesian grid and probability calculations
display(sprintf('\tConverting to cartesian grid'));
% Storage matrices
stor_part = zeros(size(gridNorth,1), size(gridEast,2));                     % Number of particles
stor_en   = zeros(size(gridNorth,1), size(gridEast,2), length(inBal.eT));   % Probability to exceed a given energy
stor_prb  = zeros(size(gridNorth,1), size(gridEast,2), length(inBal.pT));   % Energy for a probability of occurrence

% Retrieve the X position of each bomb in the grid
for iX = 1:size(gridEast,2)-1
    id(x>=gridEast(1,iX) & x<gridEast(1,iX+1),1)    = iX;
end
% Retrieve the Y position of each bomb in the grid
for iY = 1:size(gridEast,1)-1
    id(y<=gridNorth(iY,1) & y>gridNorth(iY+1,1),2)  = iY;
end

% Parse the data and calculate probabilities
count = 1;
h = waitbar(0,'Computing probabilities...');
for iX = 1:size(gridEast,2)
    for iY = 1:size(gridNorth,1)
        tmpI = id(:,1)==iX & id(:,2)==iY;
        
        % Number of particles
        stor_part(iY, iX) = length(x(tmpI));
        
        % Number of particles in pixel with energy > inBal.eT
        for j = 1:length(inBal.eT)
            stor_en(iY, iX, j) = sum(e(tmpI)>=inBal.eT(j));
        end
        
        % Energy for a given probability of occurrence
        for j = 1:length(inBal.pT)
            stor_prb(iY, iX, j) = prctile(e(tmpI), 100-inBal.pT(j)); % Here, the probability used to calculate the energy is considered as 100-percentile
        end 

        count = count+1;
        waitbar(count / (size(gridEast,2)*size(gridEast,1)));
    end
end
close(h)

%% 3 Probability calculations per distance enveloppes and radial sectors
% First, define matrices used for histograms
% a - distance
dMat        = sqrt((gridEast-inBal.vE).^2 + (gridNorth-inBal.vN).^2);               % Euclidian distance from the vent
dVec        = 0:inBal.dI:ceil(max(d)/1000)*1000;                            % Distance vector
dHist       = zeros(length(dVec), length(inBal.eT),2);                      % Matrix to plot histograms

% b - radial sector
rMat        = atan2d(gridEast-inBal.vE,gridNorth-inBal.vN);                         % Angle from the vent
idx         = rMat<0;                                                       % Angle correction to be in the interval [0 360]
rMat(idx)   = 360+rMat(idx);
rVec        = 0:inBal.rI:360;
rHist       = zeros(length(rVec), length(inBal.eT),2);                      % Matrix to plot histograms

% Second, define matrices to plot distance/radial sectors on a map
% Dim 1 = lon
% Dim 2 = lat
% Dim 3 = Energy threshold
% Dim 4 = Prob over total number of particles (1) or Prob over particles within enveloppe (2)
dP  = zeros(size(dMat,1), size(dMat,2), length(inBal.eT),2);            
rP  = zeros(size(dMat,1), size(dMat,2), length(inBal.eT),2); 

for iE = 1:length(inBal.eT)  % Loop over energy thresholds
    tmp_storD1       = zeros(size(dMat,1),size(dMat,2));   % Temp 2D storage matrix
    tmp_storD2       = zeros(size(dMat,1),size(dMat,2));   % Temp 2D storage matrix
    tmp_storR1       = zeros(size(rMat,1),size(rMat,2));   % Temp 2D storage matrix
    tmp_storR2       = zeros(size(rMat,1),size(rMat,2));   % Temp 2D storage matrix
    
    % a - distance
    for iD = 1:length(dVec)-1
        idx_dist        = dMat > dVec(iD) & dMat <= dVec(iD+1);             % Index of pixels within the distance increment
        tmp_prb         = stor_en(:,:,iE);                                  % 2D matrix containing the number of particles per pixel for a given energy threshold
        sum_impact      = sum(tmp_prb(idx_dist));                           % Sum of particles > inBal.eT and within distance increment
        nb_part_dist    = sum(stor_part(idx_dist));                         % Number of particles in the enveloppe
        
        tmp_storD1(idx_dist) = sum_impact/n.*100;                           % Probability over total number of bombs
        tmp_storD2(idx_dist) = sum_impact/nb_part_dist.*100;                % Probability over bombs in distance increment
        
        dP(:,:,iE,1)    = tmp_storD1;                                       % Fill final storage matrix
        dP(:,:,iE,2)    = tmp_storD2;
        
        dHist(iD,iE,1) = sum_impact/n;
        dHist(iD,iE,2) = sum_impact/nb_part_dist;  
    end
    
    % b - angle
    for iA = 1:length(rVec)-1
        idx_angle       = rMat > rVec(iA) & rMat <= rVec(iA+1);             % Index of pixels within the distance increment
        tmp_prb         = stor_en(:,:,iE);                                  % 2D matrix containing the number of particles per pixel for a given energy threshold
        sum_impact      = sum(tmp_prb(idx_angle));                          % Sum of particles > inBal.eT and within distance increment
        nb_part_angle   = sum(stor_part(idx_angle));                        % Number of particles in the enveloppe
        
        tmp_storR1(idx_angle) = sum_impact/n*100;                           % Probability over total number of bombs
        tmp_storR2(idx_angle) = sum_impact/nb_part_angle*100;               % Probability over bombs in distance increment
        
        rP(:,:,iE,1)    = tmp_storR1;                                       % Fill final storage matrix
        rP(:,:,iE,2)    = tmp_storR2;
        
        rHist(iA,iE,1)  = sum_impact/n;
        rHist(iA,iE,2)  = sum_impact/nb_part_angle; 
    end
end

%% 3 Save variables
display(sprintf('\tProcessing results'));

stor_en                 = stor_en./n.*100;
% Remove 0 for plotting
stor_part(stor_part==0) = nan;
stor_en(stor_en==0)     = nan;
stor_prb(stor_prb==0)   = nan;


% Write the final storage structure
project.eT              = inBal.eT;     % Energy thresholds
project.pT              = inBal.pT;     % Probability thresholds
% Cartesian data
project.gridP           = stor_prb;     % Probability (absolute)
project.gridE           = stor_en;      % Energy
project.gridN           = stor_part;    % Number of VBP per pixel
project.gridRes         = inBal.gridRes;    % Resolution of the cartesian grid
% Concentric approach
project.dI              = inBal.dI;     % Distance interval
project.dVec            = dVec;         % Vector of distances
project.dHist           = dHist;        % Probabilities as histogram
project.dP              = dP;           % Probabilities on a grid for plotting
% Radial approach
project.rI              = inBal.rI;     % Radial interval
project.rVec            = rVec;         % Vector of radial sectors
project.rHist           = rHist;        % Probabilities as histogram
project.rP              = rP;           % Probabilities on a grid for plotting
% Other data, used mostly as supporting variables for plotting
project.n               = n;            % Number of simulated VBPs
project.x               = x;            % Easting
project.y               = y;            % Northing
project.d               = d;            % Distance from the vent
project.e               = e;            % Kinetic energy
project.lat             = gridLat;          % Latitude
project.lon             = gridLon;          % Longitude
project.vLat            = vLat;         % Vent latitude
project.vLon            = vLon;         % Vent longitude
project.data            = data;         % All data loaded from input file
project.input           = inBal;        % Original input structure


save([pthout, inBal.name, '.mat'], 'project');

%% 4 Write results
display(sprintf('\tWriting results'));
% Number of particles
writeBAL([pthout, 'nb_part.txt'], gridEast, gridNorth-inBal.gridRes, stor_part);

% Probability of a given energy threshold
for i = 1:length(inBal.eT)
    writeBAL([pthout, 'prob_pixel_', num2str(inBal.eT(i)), 'J.txt'], gridEast, gridNorth-inBal.gridRes, stor_en(:,:,i));
    writeBAL([pthout, 'prob_distance_all', num2str(inBal.eT(i)), 'J.txt'], gridEast, gridNorth-inBal.gridRes, dP(:,:,i,1));
    writeBAL([pthout, 'prob_distance_zone', num2str(inBal.eT(i)), 'J.txt'], gridEast, gridNorth-inBal.gridRes, dP(:,:,i,2));
    writeBAL([pthout, 'prob_radial_all', num2str(inBal.eT(i)), 'J.txt'], gridEast, gridNorth-inBal.gridRes, rP(:,:,i,1));
    writeBAL([pthout, 'prob_radial_zone', num2str(inBal.eT(i)), 'J.txt'], gridEast, gridNorth-inBal.gridRes, rP(:,:,i,2));
end

% Energy for a given probability of occurrence
for i = 1:length(inBal.pT)
    writeBAL([pthout, 'en_', num2str(inBal.pT(i)), '%.txt'], gridEast, gridNorth-inBal.gridRes, stor_prb(:,:,i));
end



display(sprintf('Run finished at %s\n', datestr(now, 'HH:MM:SS')));


function writeBAL(out_name, X, Y, Z)
Z(isnan(Z)) = 0;
Z           = log10(Z);
Z(isinf(Z)) = 0;
Z(isnan(Z)) = -9999;

fid         = fopen(out_name,'w');
fprintf(fid,'%s\n',...
    ['ncols         ' num2str(size(X,2))],...
    ['nrows         ' num2str(size(X,1))],...
    ['xllcorner     ' num2str(min(X(1,:)))],...
    ['yllcorner     ' num2str(min(Y(:,1)))],...
    ['cellsize      ' num2str(X(1,2)-X(1,1))],...
    ['NODATA_value  ' num2str(-9999)]);
fclose(fid);

dlmwrite(out_name,Z,'-append','delimiter',' ', 'Precision', 10)