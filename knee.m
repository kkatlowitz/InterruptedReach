function [turn_time,turn_speed]=knee(coords,t,start_t,lpfilt,time_int,space_int)

%This function returns the index of the coordinates where the turn
%occurs. It works by calculating the x and y components of the
%instantaneous velocity and lowpass filtering them. It then finds the
%location of the maximal rate of change of the velocity angle (basically
%angular acceleration).
%
%outputs:
%
%turn_ind: index of the turn coordinates
%turn_speed: value proportional to the speed of the turn
%
%inputs:
%
%coords: coordinates of the trajectory
%t: time
%start_t: trial start time
%lpfilt: lowpass filter to use. If empty, function will create a default
%        filter, but this is slow. Better to create this filter outside the
%        function and then pass it as an input.
%time_int: times interval on which to look for a turn, Default [-inf,inf]
%space_int: space interval on which to look for a turn. 1x2 vector. 1st
%value is the x interval. The function will look at coordinates which are
%that many pixels in either direction from 500. 2nd value is the y.
%interval. The function will look at coordinates which are greater than
%this value. Default [0, 0]
%
%
%These arguments work pretty well for finding the turn:
%knee(coords,t,start_t,lpfilt,[b.event_time(i,3),inf],[200,300]);



if isempty(lpfilt)
    
    lpfilt=designfilt('lowpassfir','FilterOrder',50,...
        'PassbandFrequency',1,'StopbandFrequency',25, ...
        'SampleRate',100);
    
    lpfilt=lpfilt.Coefficients;
    
end

if isempty(time_int)
    time_int=[0,inf];
end


if isempty(space_int)
    space_int=[0,0];
end

%calculate velocities and filter them
xvel=diff(coords(:,1))./diff(t);
yvel=diff(coords(:,2))./diff(t);
xvelfilt=filtfilt(lpfilt,1,xvel);
yvelfilt=filtfilt(lpfilt,1,yvel);

%align the time
t=t-start_t;

%set the search interval
interval=t>time_int(1) & t<time_int(2) & ...
    (coords(:,1)>700+space_int(1) | coords(:,1)<700-space_int(1))...
    & coords(:,2)>space_int(2);

%find the point of maximal change in the velocity angle that is on the
%interval
[turn_speed,turn_ind]=max(abs(diff(unwrap(atan2(yvelfilt,xvelfilt))).*...
    interval(1:end-2)));

turn_time=t(turn_ind);

%plot(coords(:,1),coords(:,2))
%line(coords(turn_ind,1),coords(turn_ind,2),'Marker','o')


end
