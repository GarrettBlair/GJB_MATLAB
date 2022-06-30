function xdiff = jumpsmooth1D(pstamps, x, gapthresh, jumpthresh, nanradius, boxwidth)   

%REMOVES LARGE JUMPS FROM 1D POSITION TRACKING TIME SERIES DATA, THEN DOES
%BOXCAR SMOOTHING ON THE RESULT

%pstamps = time stamps of position data
%x = position data
%gapthresh = largest gap (in samples) across which to interpolate missing values after jumpremoval
%jumpthresh = movement speed that defines a large jumps (in pixels per sample)
%nanradius = radius (in samples) surrounding jump event to fill in with NaNs
%boxwidth = width of boxcar for smoothing

boxkernel = ones(boxwidth,1)/boxwidth;

%%%filter out big position jumps
if size(x,1)<size(x,2)
    x=x';
end
xdiff=[NaN; diff(x)];
jmps = find(abs(xdiff)>jumpthresh);
xdiff(jmps)=NaN;
jmps = jmps(find(jmps>nanradius));
jmps = jmps(find(jmps<(end-nanradius)));
for i=-nanradius:nanradius
    xdiff(jmps+i)=NaN;  
end

%%%fill in missing x values
x_nandex=find(isnan(xdiff)); %%indices of NaNs
nd_x=[100; diff(x_nandex)]; %%how far back was previous NaN?
x_firstnansdex=find(nd_x>1); %%indices of indices of NaNs that are the first of a sequence
nanlengths=diff(x_firstnansdex); %%lengths of the sequences
x_firstnansdex=x_firstnansdex(find(nanlengths<gapthresh)); %%indicies of NaN sequences less than half a second long
nanlengths=nanlengths(find(nanlengths<gapthresh)); %%lengths of NaN sequences less than one second long
%%construct a list of indices for NaN values that can be filled in by interpolation
filldex=[];
for dd=1:length(nanlengths)-1  
    filldex=[filldex x_nandex(x_firstnansdex(dd)):x_nandex((x_firstnansdex(dd))+nanlengths(dd)-1)];
end
if ~isnan(filldex)
    fillx=interp1(pstamps(find(~isnan(xdiff))),xdiff(find(~isnan(xdiff))),pstamps(filldex),'linear'); %%interpolate missing values
    xdiff(filldex)=fillx;
end
xdiff = conv(xdiff,boxkernel);
xdiff = xdiff((round(1+boxwidth/2)-1):(round(end-boxwidth/2)+1));

