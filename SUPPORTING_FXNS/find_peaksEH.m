function xmax=find_peaksEH(data,threshold)
% Helper function to find peaks in a given continuous valued time series x
% Usage: xmax=findpeaks(data,threshold)
% Input:
%      data     (data in time x channels/trials form or a single vector)
%      threshold (if specified returns locations of peaks at which data exceeds threshold) - optional
% Output:
%      xmax     (locations of local maxima of data in a structure array of dimensions channels/trials)

% Currently Meg's stopgap until I can get the signal processing toolbox
% installed.

if nargin < 1; error('Need data'); end;
data=change_row_to_column(data);
C=size(data,2);
pp1=[data(1,:);data(1:end-1,:)];
pp2=[data(2:end,:);data(end,:)];
xmax(1:C)=struct('loc',[]);
for ch=1:C,
   if nargin ==1
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0)];
   else
       if length(threshold) == 1
            xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0 & data(:,ch)>threshold)];
       elseif length(threshold) == 2
           xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0 & data(:,ch)>=min(threshold) & data(:,ch)<max(threshold))];
       end
   end
end


function data=change_row_to_column(data)
% Helper routine to transform 1d arrays into column vectors that are needed
% by other routines in Chronux
%
% Usage: data=change_row_to_column(data)
% 
% Inputs:
% data -- required. If data is a matrix, it is assumed that it is of the
% form samples x channels/trials and it is returned without change. If it
% is a vector, it is transformed to a column vector. If it is a struct
% array of dimension 1, it is again returned as a column vector. If it is a
% struct array with multiple dimensions, it is returned without change
% Note that the routine only looks at the first field of a struct array.
% 
% Ouputs:
% data (in the form samples x channels/trials)
%

dtmp=[];
if isstruct(data);
   C=length(data);
   if C==1;
      fnames=fieldnames(data);
      eval(['dtmp=data.' fnames{1} ';'])
      data=dtmp(:);
   end
else
  [N,C]=size(data);
  if N==1 || C==1;
    data=data(:);
  end;
end;