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
        'gridRes', [100,10],...             % Resolution
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
if length(inBal.gridRes) == 1
    % If no grid resolution specified for the grid used for the
    % concentric/radial approach is defined, use the same as the pixel
    % approach
    inBal.gridRes = [inBal.gridRes, inBal.gridRes];
end

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
r       = atan2d(x-inBal.vE,y-inBal.vN);
r(r<0)  = 360+r(r<0);
id      = zeros(size(x,1),2);           % Vectors of indices for cartesian grid

data    = data(:,3:end);                % Colums:
                                        % 1: Landing altitude (m a.s.l.)
                                        % 2: Mass (kg)
                                        % 3: Diameter (m)
                                        % 4: Kinetic energy (kJ)
                                        % 5: Landing angle (deg)
                                        % 6: Ejection andgle(deg)
                                        % 7: Flight time (sec)

%% Define storage
% Grid
grid        = struct;   % Storage for the grid data
% grid.pixelP
% grid.pixelE
coorRef     = struct;   % Storage for the reference coordinates of the grid
concentric  = struct;   % Storage for the concentric distance
radial      = struct;   % Storage for the radial sectors
                                        
%% Define the cartesian grid                                                                               
% Create cartesian grid for the pixel approach

[coorRef.pixelEast, coorRef.pixelNorth]     = meshgrid(min(x)-inBal.gridRes(1):inBal.gridRes(1):max(x)+inBal.gridRes(1),...
    min(y)-inBal.gridRes(1):inBal.gridRes(1):max(y)+inBal.gridRes(1));
coorRef.pixelNorth                          = flipud(coorRef.pixelNorth);


% New approach using histcounts2
%% --------------------------------------------
% Define grid in projected and geographic coordinates
[pixel.east, pixel.north] = meshgrid(min(x)-inBal.gridRes(1):inBal.gridRes(1):max(x)+inBal.gridRes(1),...
    min(y)-inBal.gridRes(1):inBal.gridRes(1):max(y)+inBal.gridRes(1));
pixel.northing            = flipud(coorRef.pixelNorth);
[pixel.lat, pixel.lon]    = utm2ll(pixel.east, pixel.north, ones(size(coorRef.pixelEast)).*inBal.vZ);

% Define storage
pixel.E     = zeros(size(pixel.east,1), size(pixel.east,2), length(inBal.pT));
pixel.N     = zeros(size(pixel.east,1), size(pixel.east,2), length(inBal.eT));
pixel.Pabs  = zeros(size(pixel.east,1), size(pixel.east,2), length(inBal.eT));
pixel.Prel  = zeros(size(pixel.east,1), size(pixel.east,2), length(inBal.eT));

% Bin the total number of particles
[pixel.Nt,~,~,xi,yi] = histcounts2(x,y,pixel.easting(1,:), pixel.northing(:,1));

% Calculate the number of particles per pixel per energy threshold and the
% absolute and relative probabilities (in %)
for iE = 1:length(inBal.eT)
    pixel.N(:,:,iE)     = histcounts2(x(e>inBal.eT(iE)), y(e>inBal.eT(iE)), pixel.east(1,:), pixel.north(:,1));
    pixel.Pabs(:,:,iE)  = pixel.N(:,:,iE) ./ n .* 100;
    pixel.Prel(:,:,iE)  = pixel.N(:,:,iE) ./ pixel.Nt .* 100;
end

% Calculate the energy for a given exceedance probability
for iE = 1:size(pixel.N, 2)
    for iN = 1:size(pixel.N, 1)
        for iP = 1:length(inBal.pT)
            pixel.E(iN, iE, iP) = prctile(e(xi==iE & yi==iN), 100-inBal.pT(iP));
        end
    end
end






%% --------------------------------------------



% Create cartesian grid for the concentric/radial approaches
[coorRef.crEast, coorRef.crNorth]           = meshgrid(min(x)-inBal.gridRes(2):inBal.gridRes(2):max(x)+inBal.gridRes(2),...
    min(y)-inBal.gridRes(2):inBal.gridRes(2):max(y)+inBal.gridRes(2));
coorRef.crNorth                             = flipud(coorRef.crNorth);

% Convert to geographical coordinates for plotting
[coorRef.pixelLat, coorRef.pixelLon]        = utm2ll(coorRef.pixelEast, coorRef.pixelNorth, ones(size(coorRef.pixelEast)).*inBal.vZ);
[coorRef.crLat, coorRef.crLon]              = utm2ll(coorRef.crEast, coorRef.crNorth,ones(size(coorRef.crEast)).*inBal.vZ);

% Same for the vent
[vLat,vLon]                                 = utm2ll(inBal.vE, inBal.vN, inBal.vZ);

% Adjust sizes of storage matrices to grids
grid.pixelP         = zeros(size(coorRef.pixelEast, 1), size(coorRef.pixelEast, 2), length(inBal.eT), 2);   % Pixel approach, probability
grid.pixelE         = zeros(size(coorRef.pixelEast, 1), size(coorRef.pixelEast, 2), length(inBal.pT));      % Pixel approach, energy
grid.pixelN         = zeros(size(coorRef.pixelEast, 1), size(coorRef.pixelEast, 2), length(inBal.pT));      % Pixel approach, number of particles per energy threshold
grid.pixelNt        = zeros(size(coorRef.pixelEast, 1), size(coorRef.pixelEast, 2));                        % Pixel approach, total number of particles
grid.concentricP    = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2), length(inBal.eT), 2);         % Concentric approach, probability
grid.concentricE    = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2), length(inBal.eT));            % Concentric approach, energy
grid.concentricN    = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2));                              % Concentric approach, number of particles
grid.concentricI    = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2));                              % Concentric approach, index
grid.radialP        = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2), length(inBal.eT), 2);         % Radial approach, probability
grid.radialE        = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2), length(inBal.eT));            % Radial approach, energy
grid.radialN        = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2));                              % Radial approach, number of particles
grid.radialI        = zeros(size(coorRef.crEast, 1), size(coorRef.crEast, 2));                              % Radial approach, index

%% 2 Conversion to cartesian grid and probability calculations
display(sprintf('\tConverting to cartesian grid'));
% Storage matrices
%stor_part = zeros(size(pixelNorth,1), size(pixelEast,2));                     % Number of particles
%stor_en   = zeros(size(pixelNorth,1), size(pixelEast,2), length(inBal.eT));   % Probability to exceed a given energy
%stor_prb  = zeros(size(pixelNorth,1), size(pixelEast,2), length(inBal.pT));   % Energy for a probability of occurrence

% Retrieve the X position of each bomb in the grid
for iX = 1:size(coorRef.pixelEast,2)-1
    id(x>=coorRef.pixelEast(1,iX) & x<coorRef.pixelEast(1,iX+1),1)    = iX;
end
% Retrieve the Y position of each bomb in the grid
for iY = 1:size(coorRef.pixelEast,1)-1
    id(y<=coorRef.pixelNorth(iY,1) & y>coorRef.pixelNorth(iY+1,1),2)  = iY;
end

% Parse the data and calculate probabilities
%count = 1;
%h = waitbar(0,'Computing probabilities...');
for iX = 1:size(coorRef.pixelEast,2)
    for iY = 1:size(coorRef.pixelNorth,1)
        tmpI = id(:,1)==iX & id(:,2)==iY;
        
        % Number of particles
        grid.pixelNt(iY, iX) = length(x(tmpI));
        
        % Number of particles in pixel with energy eT
        for j = 1:length(inBal.eT)
            grid.pixelN(iY, iX, j) = sum(e(tmpI)>=inBal.eT(j));
        end
        
        % Energy for a given probability of occurrence
        for j = 1:length(inBal.pT)
            grid.pixelE(iY, iX, j) = prctile(e(tmpI), 100-inBal.pT(j)); % Here, the probability used to calculate the energy is considered as 100-percentile
        end 

%        count = count+1;
%        waitbar(count / (size(coorRef.pixelEast,2)*size(coorRef.pixelEast,1)));
    end
end
%close(h)

% Calculate probabilities in percent
%grid.pixelP(:,:,:,2) = grid.pixelP(:,:,:,1) ./ repmat(grid.pixelN,1,1,length(inBal.eT)) .* 100;   % Relative probabilities
%grid.pixelP(:,:,:,1) = grid.pixelP(:,:,:,1) ./ n .* 100;                                            % Absolute probabilities

%% 3 Probability calculations per distance enveloppes and radial sectors
% First, define matrices used for histograms
% a - distance
%dMat 
%dVec 
%dHist 
%grid.concentricI    = sqrt((coorRef.crEast-inBal.vE).^2 + (coorRef.crNorth-inBal.vN).^2);               % Euclidian distance from the vent 
% concentric.vec      = 0:inBal.dI:ceil(max(d)/1000)*1000;                            % Distance vector 
% concentric.hist     = zeros(length(concentric.vec), length(inBal.eT),2);                      % Matrix to plot histograms
% concentric.vec      = 0:inBal.dI:ceil(max(d)/1000)*1000;                            % Distance vector 

concentric.bin      = 0:inBal.dI:ceil(max(d)/1000)*1000;
concentric.Nt       = zeros(length(concentric.bin),1);
concentric.N        = zeros(length(concentric.bin), length(inBal.eT));
concentric.E        = zeros(length(concentric.bin), length(inBal.pT));
concentric.P        = zeros(length(concentric.bin), length(inBal.eT), 2);

for iC = 1:length(concentric.bin)-1
    concentric.Nt(iC) = nnz(d > concentric.bin(iC) & d <= concentric.bin(iC+1));
    for iE = 1:length(inBal.eT)
        concentric.N(iC,iE) = nnz(d > concentric.bin(iC) & d <= concentric.bin(iC+1) & e > inBal.eT(iE));
    end
    for iP = 1:length(inBal.pT)
        concentric.E(iC,iP) = prctile(e(d > concentric.bin(iC) & d <= concentric.bin(iC+1)), 100-inBal.pT(iP));
    end
end

% b - radial sector
%rMat 
%idx
%rMat
%rVec
%rHist
grid.radialI        = atan2d(coorRef.crEast-inBal.vE,coorRef.crNorth-inBal.vN);                         % Angle from the vent                                                       % Angle correction to be in the interval [0 360]
grid.radialI(grid.radialI<0) = 360+grid.radialI(grid.radialI<0);
radial.vec          = 0:inBal.rI:360;
radial.hist         = zeros(length(radial.vec), length(inBal.eT),2);                      % Matrix to plot histograms

% Second, define matrices to plot distance/radial sectors on a map
% Dim 1 = lon
% Dim 2 = lat
% Dim 3 = Energy threshold
% Dim 4 = Prob over total number of particles (1) or Prob over particles within enveloppe (2)
%dP  = zeros(size(dMat,1), size(dMat,2), length(inBal.eT),2);            
%rP  = zeros(size(dMat,1), size(dMat,2), length(inBal.eT),2); 

for iE = 1:length(inBal.eT)  % Loop over energy thresholds
%     tmp_storD1       = zeros(size(dMat,1),size(dMat,2));   % Temp 2D storage matrix
%     tmp_storD2       = zeros(size(dMat,1),size(dMat,2));   % Temp 2D storage matrix
%     tmp_storR1       = zeros(size(rMat,1),size(rMat,2));   % Temp 2D storage matrix
%     tmp_storR2       = zeros(size(rMat,1),size(rMat,2));   % Temp 2D storage matrix
    
    % a - distance
    for iD = 1:length(concentric.bin)-1
        iConcentric     = d > concentric.vec(iD) & d <= concentric.vec(iD+1);
        concentric.N    = nnz(iConcentric);
        
        
        idx_dist        = grid.concentricI > concentric.vec(iD) & grid.concentricI <= concentric.vec(iD+1);             % Index of pixels within the distance increment
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
project.lat             = pixelLat;          % Latitude
project.lon             = pixelLon;          % Longitude
project.vLat            = vLat;         % Vent latitude
project.vLon            = vLon;         % Vent longitude
project.data            = data;         % All data loaded from input file
project.input           = inBal;        % Original input structure


save([pthout, inBal.name, '.mat'], 'project');

%% 4 Write results
display(sprintf('\tWriting results'));
% Number of particles
writeBAL([pthout, 'nb_part.txt'], pixelEast, pixelNorth-inBal.gridRes, stor_part);

% Probability of a given energy threshold
for i = 1:length(inBal.eT)
    writeBAL([pthout, 'prob_pixel_', num2str(inBal.eT(i)), 'J.txt'], pixelEast, pixelNorth-inBal.gridRes, stor_en(:,:,i));
    writeBAL([pthout, 'prob_distance_all', num2str(inBal.eT(i)), 'J.txt'], pixelEast, pixelNorth-inBal.gridRes, dP(:,:,i,1));
    writeBAL([pthout, 'prob_distance_zone', num2str(inBal.eT(i)), 'J.txt'], pixelEast, pixelNorth-inBal.gridRes, dP(:,:,i,2));
    writeBAL([pthout, 'prob_radial_all', num2str(inBal.eT(i)), 'J.txt'], pixelEast, pixelNorth-inBal.gridRes, rP(:,:,i,1));
    writeBAL([pthout, 'prob_radial_zone', num2str(inBal.eT(i)), 'J.txt'], pixelEast, pixelNorth-inBal.gridRes, rP(:,:,i,2));
end

% Energy for a given probability of occurrence
for i = 1:length(inBal.pT)
    writeBAL([pthout, 'en_', num2str(inBal.pT(i)), '%.txt'], pixelEast, pixelNorth-inBal.gridRes, stor_prb(:,:,i));
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