%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polar_analyticSig.m
%
% Last modified by Huailiang Li (09 Dec 2019)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Program polar_analyticSig is a Matlab function to perform complex-valued 
% polarization analysis using a analytic signal with its complex conjugate
% Modified by Li
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
% dtac          a matrix of size (3 x Nt), where Nt is the number of time 
%               samples; the three rows of the matrix are the vertical, 
%               east, and north components in that order
% wndo          the length of the window to calculate polarization 
%               parameters in samples
%
% Output:
% azim          time series of azimuth angles in clockwise degrees from 
%               north
% incd          time series of incidence angles in degrees from vertical
% ellip         time series of ellipticity - intermediate axis of ellipsoid 
%               divided by major axis
%
% Dependency:
% Uses program convsm1d.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [azim incd maxeig,Dlp,Dpp] = polar_analyticSig(dataX, dataY, dataZ,wndo)

clearvars -except dataX dataY dataZ wndo;

dtac = [dataZ'; dataX'; dataY'];

% unpack the matrix of three-component motion
pv6z = dtac(1,:);
pv6e = dtac(2,:);
pv6n = dtac(3,:);


% length of the time series
tln = length(pv6z);

nsall = tln ;   %  total number of samples to analyze

% Apply a moving window analysis to the 3 component data.
%   Moving window loop

nswin = wndo;  %  number of samples in a time window
npshift = 1;  % number of samples to shift over (recommend)
kfin = fix((nsall-nswin)/npshift); % number of time windows considered


% Apply a moving window analysis to the 3 component data.
%   Moving window loop
npshift = 1;  % number of samples to shift over (recommend)
kfin = fix((nsall-nswin)/npshift); % number of time windows considered

mxde1=0.;

for k=1:1:kfin
    
    nwinst = (k-1)*npshift+1;  % start of time window
    nwinend = nwinst+nswin-1;    % end of time window
    a = csigm(pv6e,pv6n,pv6z,nwinst,nwinend);  % signal matrix %Load a time portion of the signals into the signal matrix.
    a = a'; % Row--->Column
    
    
    zh = hilbert(a(:,3)); %construct analytic signal
    nh = hilbert(a(:,2));
    eh = hilbert(a(:,1));
    
    aa = [eh,nh,zh];

    %       COV_XX  COV_XY  COV_XZ
    % M =   COV_YX  COV_YY  COV_YZ
    %       COV_ZX  COV_ZY  COV_ZZ
    
%     CovXX = eh*(conj(eh))';
%     CovXY = eh*(conj(nh))';
%     CovXZ = eh*(conj(zh))';
%     
%     CovYX = nh*(conj(eh))';
%     CovYY = nh*(conj(nh))';
%     CovYZ = nh*conj(zh)';
%     
%     CovZX = zh*(conj(eh))';
%     CovZY = zh*(conj(nh))';
%     CovZZ = zh*(conj(zh))';
%     
%     M = [CovXX,CovXY,CovXZ;CovYX,CovYY,CovYZ;CovZX,CovZY,CovZZ];
%     
     M = cov(aa);
    [v1,d1]=eig(M);
    
    [v,d]=order(v1,d1);     % order the eigenvalues and eigenvectors,
    % in descending order
    
    vp(:,k) = v(:,1);  % used to plot polarization trajectory/largest eigenvalues and eigenvectors

    %Calculate azimuths from the eigenvector components

    %%
    v_norm = v(:,1)./abs(v(:,1));
%     v_norm = v(:,1);
    
    alph = 0:0.01:179.99;
    cis = cosd(alph)+ 1i*sind(alph);
    
    R_Temp = sqrt((real(v_norm(1)*cis)).^2+ (real(v_norm(2)*cis)).^2+ (real(v_norm(3)*cis)).^2);
    
    [~, Angle_index]= max(R_Temp);
    
    ans(:,1) = real(v(:,1));
    
    v(1,1) = v(1,1)* (cosd(alph(Angle_index))+ 1i*sind(alph(Angle_index)));
    v(2,1) = v(2,1)* (cosd(alph(Angle_index))+ 1i*sind(alph(Angle_index)));
    v(3,1) = v(3,1)* (cosd(alph(Angle_index))+ 1i*sind(alph(Angle_index)));
     
    ans(:,2) = real(v(:,1));

    ang1(k)=atan2d(real(v(1,1))*sign(real(v(3,1))),real(v(2,1))*sign(real(v(3,1))));
    
    de1(k)=d(1);
    maxeig(k) = de1(k);
    mxde1=max(mxde1,de1(k));  % find the maximum values

    Dlp(k) = 1- (d(2)+d(3))/d(1);  %rectilinearity, r, or degree of linear polarization,as for pure body waves.
                                  %In the case of an elliptical polarized Rayleigh wave, 
                                  %the  lamuda1 and lamuda2 are much greater than lamuda3, 
                                  %so their degree of rectilinearity is small, (L > 0.5).
    
    Dpp(k) = 1 - 2*d(3)/(d(1)+d(2));% The degree of planarity,
    
   % calculate the ellipticity
    ellip(k) = real(sqrt(d(2)/d(1)));
    
    %Calculate the apparent incidence angles from the eigenvectors.
    
    vang1(k)=atan2d(sqrt(real(v(1,1))^2+ real(v(2,1))^2), real(v(3,1)));
    
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

azim = ang1;
incd = vang1;




