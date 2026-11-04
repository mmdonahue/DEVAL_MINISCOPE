function linCoords = linearize_fig8_coords(coords, varargin)
% function linCoords = linearize_fig8_coords(coords, binSz)
% 
% PURPOSE:
%   Linearize the coordinates from the figure 8 maze. See note below for
%   which bins correspond to which parts of the maze for recommended bin size.
% 
% INPUT:
%   coords = CONVERTED 2D coords array with
%       coords(:,1) = time in seconds
%       coords(:,2) = x coordinates in cm
%       coords(:,3) = y coordinates in cm
%   binSz = size in cm of which linear bin
%       Default (recommended) is 4 cm
% 
% OUTPUT:
%   linCoords = array with linearized coordinates of the mouse in bins
% 
% NOTES:
%   - This code is based on the 56 x 56 cm figure 8 maze. This can be
%   changed in the initial lines of this code.
%   - The bins that correspond to various positions are described below:
%       1-14: stem
%       15-21: choice point -> East stem entrance
%       22-35: East stem -> East reward
%       36-42: East reward -> delay zone
%       43-49: choice point -> West stem entrance
%       50-63: West stem -> West reward
%       34-70: West reward -> delay zone
% 
% MMD
% 05/2025
% Ewell Lab


%% MAZE SIZE

if nargin ==1
    binSz = 4;
else
    binSz = varargin{1};
end %set bin size

x_arenaSize = 56;
y_arenaSize = 56;

%% BIN CENTERS FOR LINEARIZING

% I KNOW this looks crazy but it should work
stemX = linspace(binSz/2, x_arenaSize - (binSz/2), x_arenaSize/binSz)';
stemY = repmat(y_arenaSize/2, y_arenaSize/binSz, 1);

%turn W
WturnX = repmat(x_arenaSize - (binSz/2), x_arenaSize/binSz/2, 1);
WturnY = transpose(y_arenaSize/2 + binSz/2:binSz:y_arenaSize - binSz/2);
WstemX = flip(stemX);
WstemY = repmat(y_arenaSize - binSz/2, y_arenaSize/binSz, 1);
WreturnX = repmat(binSz/2, x_arenaSize/binSz/2, 1);
WreturnY = flip(WturnY);

%turn E
EturnX = repmat(x_arenaSize - (binSz/2), x_arenaSize/binSz/2, 1);
EturnY = transpose(flip(binSz/2:binSz:y_arenaSize/2 - binSz/2));
EstemX = flip(stemX);
EstemY = repmat(binSz/2, y_arenaSize/binSz, 1);
EreturnX = repmat(binSz/2, x_arenaSize/binSz/2, 1);
EreturnY = flip(EturnY);

%put it all together
binCtrs(:,1) = [stemX; EturnX; EstemX; EreturnX; WturnX; WstemX; WreturnX];
binCtrs(:,2) = [stemY; EturnY; EstemY; EreturnY; WturnY; WstemY; WreturnY];

%% DO CALCULATIONS

linCoords(:,1) = coords(:,1); %time
closeBins = dsearchn(binCtrs, coords(:,2:3));
linCoords(:,2) = closeBins;

end %function