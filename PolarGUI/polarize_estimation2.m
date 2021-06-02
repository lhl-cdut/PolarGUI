%"polarize.m" is a complicated times series program that reads in SAC (Seismic Analysis Code) binary files,
%does some filtering of the data, and then performs an eigenvalue/eigenvector
%decomposition of the particle motion to investigate estimates of azimuth and incidence angle
%from arrivals in the three component waveforms.  Plots include the input time series,
%filtered time series, eigenvalues as a function of time,
%azimuth and incidence angle calculations from the eigenvectors, and rose plots of the same.
%Links are given to functions called in "polarize.m".



function [azi,inci,maxeig,Dlp,Dpp] = polarize_estimation2(X,Y,Z,delt,ttot,twin)
%
%  function polarize(station,delt,ttot,twin,hilb,flp,fhi)
%
%  Program to read in 3 component waveform data
%  Create the covariance matrix for a moving time window
%  Find the principal components and infer polarization
%
% input:
%   station = station name for sacfile prefix
%   delt = sampling interval
%   ttot = total number of seconds to analyze in traces
%   twin = time window length, each time shift will be 1/2 of the
%      window length
%   hilb = 0, no hilbert transform of vertical component
%      = 1, hilbert transform
%   flp = low frequency corner frequency of a 2nd order butterworth
%     filter used to filter the data, if 0, then no filtering
%   fhi = hi frequency corner frequency of the filter

%Only the first part of the  filename is given.  This part reconstructs the rest.

clearvars -except X Y Z delt ttot twin;
e=dmean(X');           % remove the mean from each
n=dmean(Y');
z=dmean(Z');

% e = X';
% n = Y';
% z = Z';


nsall = fix(ttot/delt) ;   %  total number of samples to analyze
for m = 1:1:nsall %Change the X-axis t
    t(m) = m*delt;
end


% Apply a moving window analysis to the 3 component data.
%   Moving window loop

nswin = fix(twin/delt);  %  number of samples in a time window
% 
npshift = 1;  % number of samples to shift over (recommend)
kfin = fix((nsall-nswin)/npshift); % number of time windows considered

mxde1=0.;
mxde2=0.;
mxde3=0.;
for k=1:1:kfin
    nwinst = (k-1)*npshift+1;  % start of time window
    nwinend = nwinst+nswin-1;    % end of time window
    ta = csigm(e,n,z,nwinst,nwinend);  % signal matrix %Load a time portion of the signals into the signal matrix.   
    
    % remove local mean over the window length
%     ta(1,:) = a(1,:)-convsm1d(a(1,:),nwin);
%     ta(2,:) = a(2,:)-convsm1d(a(2,:),nwin);
%     ta(3,:) = a(3,:)-convsm1d(a(3,:),nwin);   
    
    aa = ta'; %Column->Row
    
    %     c = (a'*a)./(size(a,1)-1);
    aa = aa-repmat(mean(aa),nswin,1);
    c = (aa'*aa)./(nswin-1);   % Construct the covariance matrix.


    [v1,d1]=eig(c);      % eigenvalue/eigenvectors,Find the eigenvalues and eigenvectors of the signal matrix.
    %Order the resulting eigenvectors/eigenvalues in decending order.

    
    
    [v,d]=order(v1,d1);     % order the eigenvalues and eigenvectors,
    % in descending order
    
    vp(:,k) = v(:,1);  % used to plot polarization trajectory/largest eigenvalues and eigenvectors
    
    ang1(k)=atan2d(v(1,1)*sign(v(3,1)),v(2,1)*sign(v(3,1))); % azimuth for first  eigenvalue/  azi=atan(u11/21)
    
    ang2(k)=atan2(v(2,2),v(3,2)) * 180/pi; % azimuth for second eigenvalue
    ang3(k)=atan2(v(2,3),v(3,3)) * 180/pi; % azimuth for third  eigenvalue
    de1(k)=d(1);
    maxeig(k) = de1(k);
    de2(k)=d(2);
    de3(k)=d(3);
    mxde1=max(mxde1,de1(k));  % find the maximum values
    mxde2=max(mxde2,de2(k));
    mxde3=max(mxde3,de3(k));
    %Calculate the apparent incidence angles from the eigenvectors.
    
    vang1(k)=atan2d(sqrt(v(1,1)^2+ v(2,1)^2), v(3,1));
    
    vang2(k)=acos(abs(v(1,2)))* 180/pi;
    vang3(k)=acos(abs(v(1,3)))* 180/pi;
    t2(k)=delt*(nwinst-1);    % assign time for this window to the window start
    
    Dlp(k) = 1- (d(2)+d(3))/d(1);  %rectilinearity, r, or degree of linear polarization,as for pure body waves.
                                  %In the case of an elliptical polarized Rayleigh wave, 
                                  %the  lamuda1 and lamuda2 are much greater than lamuda3, 
                                  %so their degree of rectilinearity is small, (L > 0.5).
    
    Dpp(k) = 1 - 2*d(3)/(d(1)+d(2));% The degree of planarity,
    
%    dirVX(k) = rad2deg(atan2(vp(2,k),vp(1,k))); % Vector direction angle X/Y plane
%    dirVY(k) = rad2deg(atan2(vp(3,k),vp(1,k))); % Vector direction angle X/Z plane
%    dirVZ(k) = rad2deg(atan2(vp(3,k),vp(2,k))); % Vector direction angle Y/Z plane
    
   % calculate the ellipticity
    ellip(k) = real(sqrt(d(2)/d(1)));


    dirVX(k) = rad2deg(acos(vp(1,k)/norm(vp(:,k)))); % Vector direction cosine X-axis
    dirVY(k) = rad2deg(acos(vp(2,k)/norm(vp(:,k)))); % Vector direction cosine Y-axis
    dirVZ(k) = rad2deg(acos(vp(3,k)/norm(vp(:,k)))); % Vector direction cosine Z-axis
    
    if ang1(k) < 0
        ang1(k) = ang1(k)+360;
    elseif  ang1(k) > 360
        ang1(k) = ang1(k)-360;
    end
    
    if vang1(k) < 0
        vang1(k) = vang1(k)+180;
    elseif  vang1(k) > 180
        vang1(k) = vang1(k)-180;
    end

    
end;

azi = ang1;
inci = vang1;
%=======================================================================


