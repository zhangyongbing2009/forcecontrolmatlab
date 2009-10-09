function figh = CreateWindow (position,Width,Height,figtitle)
%CREATEWINDOW to create new window with given dimensions     
% figh = Create_Window (dx,dy,position,figtitle)
%
% figh     : figure handle
% Width    : window width in pixels
% Height   : window height in pixels
% position : string indicating position of window on screen (optional)
%              'cen' center
%              'ulc' upper left corner
%              'llc' lower left corner
%              'urc' upper right corner
%              'lrc' lower right corner
% figtitle : title of window (otional)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          OpenSees Navigator                          %%
%% Matlab Engineering Toolbox for Analysis of Structures using OpenSees %%
%%                                                                      %%
%%                   Andreas Schellenberg & Tony Yang                   %%
%%        andreas.schellenberg@gmail.com, yangtony2004@gmail.com        %%
%%                                                                      %%
%%    Department of Civil and Environmental Engineering, UC Berkeley    %%
%%   (C) Copyright 2004, The Regents of the University of California    %%
%%                         All Rights Reserved                          %%
%%                                                                      %%
%%   Commercial use of this program without express permission of the   %%
%%     University of California, Berkeley, is strictly prohibited.      %%
%%     See Help -> OpenSees Navigator Disclaimer for information on     %%
%%  usage and redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% $Revision: 217 $
% $Date: 2007-01-02 11:09:29 -0800 (Tue, 02 Jan 2007) $
% $URL: $ 

% get the screen size
SS = get(0,'screensize');

if nargin < 1
    position = 'cen';
    Width = SS(3)/2;
    Height = SS(4)/2;
elseif nargin < 2
    Width = SS(3)/2;
    Height = SS(4)/2;
elseif nargin < 3
    Height = SS(4)/2;
end

% get the display position
switch position
   case 'cen'
      Position = [(SS(3)-Width)/2, (SS(4)-Height)/2, Width, Height];
   case 'ulc'
      Position = [1, SS(4)-Height, Width, Height];
   case 'llc'
      Position = [1, 1, Width, Height];
   case 'urc'
      Position = [SS(3)-Width, SS(4)-Height, Width, Height];
   case 'lrc'
      Position = [SS(3)-Width, 1, Width, Height];
end

% create figure handle
figh = figure('Units','pixel','Position',Position);      

if (nargin==4)
   set(figh,'NumberTitle','off', ...
      'Name',figtitle);
end