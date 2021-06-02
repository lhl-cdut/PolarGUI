%"polarize.m" is a complicated times series program that reads in SAC (Seismic Analysis Code) binary files,
%does some filtering of the data, and then performs an eigenvalue/eigenvector
%decomposition of the particle motion to investigate estimates of azimuth and incidence angle
%from arrivals in the three component waveforms.  Plots include the input time series,
%filtered time series, eigenvalues as a function of time,
%azimuth and incidence angle calculations from the eigenvectors, and rose plots of the same.
%Links are given to functions called in "polarize.m".

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polarize_estimation.m
%
% PROGRAMMER:
% Huailiang Li
%
% Last revision date:
% 09 Dec 2019
% Last modified by ()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azi,inci,maxeig,Dlp,Dpp] = polarize_estimation(X,Y,Z,delt,ttot,twin)
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


nsall = fix(ttot/delt);   %  total number of samples to analyze
for m = 1:1:nsall %Change the X-axis t
    t(m) = m*delt;
end


% Apply a moving window analysis to the 3 component data.
%   Moving window loop

nswin = fix(twin/delt);  %  number of samples in a time window

npshift = 1;  % number of samples to shift over (recommend)
kfin = fix((nsall-nswin)/npshift); % number of time windows considered

mxde1=0.;

for k=1:1:kfin
    nwinst = (k-1)*npshift+1;  % start of time window
    nwinend = nwinst+nswin-1;    % end of time window
    a = csigm(e,n,z,nwinst,nwinend);  % signal matrix %Load a time portion of the signals into the signal matrix.
    a = a'; % Row->Col
%     %       COV_XX  COV_XY  COV_XZ
%     % M =   COV_YX  COV_YY  COV_YZ
%     %       COV_ZX  COV_ZY  COV_ZZ
%     CovXX = cov(a(:,1),a(:,1),1);
%     CovYY = cov(a(:,2),a(:,2),1);
%     CovZZ = cov(a(:,3),a(:,3),1);
% %     CovXX = var(a(:,1),1);
% %     CovYY = var(a(:,2),1);
% %     CovZZ = var(a(:,3),1);
%     
%     CovXY = cov(a(:,1),a(:,2),1);
%     CovXZ = cov(a(:,1),a(:,3),1);
%     
%     CovYX = cov(a(:,2),a(:,1),1);
%     CovYZ = cov(a(:,2),a(:,3),1);
%     
%     CovZX = cov(a(:,3),a(:,1),1);
%     CovZY = cov(a(:,3),a(:,2),1);
%     
%     M = [CovXX,CovXY,CovXZ;CovYX,CovYY,CovYZ;CovZX,CovZY,CovZZ];
    
      M = cov(a);

    %Order the resulting eigenvectors/eigenvalues in decending order.
    [v1,d1]=eig(M);
    [v,d]=order(v1,d1);     % order the eigenvalues and eigenvectors,
    % in descending order
    
    vp(:,k) = v(:,1);  % used to plot polarization trajectory/largest eigenvalues and eigenvectors
    
    ang1(k)=atan2d(v(1,1)*sign(v(3,1)),v(2,1)*sign(v(3,1))); % azimuth for first  eigenvalue/  azi=atan(u11/21)
    
%     de1(k)=d(1);
    maxeig(k) = d(1);

    %mxde1=max(mxde1,de1(k));  % find the maximum values
    
    Dlp(k) = 1- (d(2)+d(3))/d(1);  %rectilinearity, r, or degree of linear polarization,as for pure body waves.
                                  %In the case of an elliptical polarized Rayleigh wave, 
                                  %the  lamuda1 and lamuda2 are much greater than lamuda3, 
                                  %so their degree of rectilinearity is small, (L > 0.5).
    
    Dpp(k) = 1 - 2*d(3)/(d(1)+d(2));% The degree of planarity,

    %Calculate the apparent incidence angles from the eigenvectors.
%     vang1(k)=acos(abs(v(1,1)))* 180/pi; %angle from the vertical /incidence = acos(|u11|)
    vang1(k)=atan2d(sqrt(v(1,1)^2+ v(2,1)^2), v(3,1));
    
    % calculate the ellipticity
    ellip(k) = real(sqrt(d(2)/d(1)));
    
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


