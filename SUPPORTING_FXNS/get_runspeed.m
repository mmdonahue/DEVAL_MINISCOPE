function runSpeed = get_runspeed(coords)
% function runSpeed = get_runspeed(coords)
% 
% PURPOSE:
%   Get the runspeed from 2-d coords.
% 
% INPUT:
%   coords = (output from convert_coords) array with:
%      coords(:,1) = time in seconds
%      coords(:,2) = x coordinates in cm
%      coords(:,3) = y coordinates in cm
% 
% OUTPUT:
%   runSpeed = array with:
%       runSpeed(:,1) = time in s
%       runSpeed(:,2) = speed in cm/s
% 
% MMD for Ewell lab
% 05/2025
% based on code from JBT for Colgin lab

fps = 1 ./ (mean(diff(coords(:,1)))); %should be 29.97

xDiff = diff(coords(:,2)); %distance mouse traveled in x dim in between each frame
yDiff = diff(coords(:,3)); %distance mouse traveled in y dim in between each frame

tmpSpd = hypot(xDiff, yDiff); %speed in cm/frame
tmpSpd = tmpSpd * fps; %convert to cm/s
tmpSpd = [tmpSpd(1); tmpSpd]; %assume same speed for first time bin so same length as coords

runSpeed(:,1) = coords(:,1);
runSpeed(:,2) = tmpSpd;

end %function

