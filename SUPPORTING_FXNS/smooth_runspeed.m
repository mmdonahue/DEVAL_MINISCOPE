function smRs = smooth_runspeed(instRs, varargin)
% function smRs = smooth_runspeed(instRs, gWinDur, gWinStd)
%
% PURPOSE:
%  To apply Gaussian smoothing to the runspeed time-series.
%
% INPUT:
%   instRs = runspeed matrix as outputted by get_runspeed
%  gWinDur = Gaussian window size in seconds
%            - If absent or empty, default is 250 ms
%  gWinStd = Gaussian window standard deviation in seconds
%            - If absent or empty, default is gWinDur / 6
%
% OUTPUT:
%    smRs = runspeed matrix, same size as input, but with the
%           speed column (:,2) now convolved with a Gaussian kernel
%
% JB Trimper modified by MMD for Ewell lab
% 09/2019
% Colgin Lab

if nargin == 1 
    gWinSz = 250 / 1000;
    gWinStd = gWinSz ./ 6;
end

if nargin < 3 || isempty(gWinStd)
    gWinSz = varargin{1};
    gWinStd = gWinSz ./ 6;
end

if nargin == 3
    gWinSz = varargin{1};
    gWinStd = varargin{2};
end

fps = 1/(mean(diff(instRs(:,1)))); %should be 29.97

%Make the output matrix and copy the time-stamps
smRs = zeros(size(instRs));
smRs(:,1) = instRs(:,1);

%convert window duration to # of samples based on 30 Hz frame rate
gWinSz = round(fps*gWinSz);

%convert window std to # of samples
gWinStd = round(fps*gWinStd);

%Make the kernel and convolve it with the instantaneous runspeed
gKrnl = gausskernel(gWinSz, gWinStd);
smRs(:,2) = conv(instRs(:,2),gKrnl,'same');


end %function