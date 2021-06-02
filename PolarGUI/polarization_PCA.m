
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PROGRAM:
% polarization_PCA.m
%
% PROGRAMMER:
% Huailiang Li
%
% Last revision date:
% 09 Dec 2019
% Last modified by ()
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [azi,inci,maxeig,Dlp,Dpp] = polarization_PCA(X,Y,Z,delt,twin) 
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

clearvars -except X Y Z delt twin;

% e=dmean(X');           % remove the mean from each 
% n=dmean(Y'); 
% z=dmean(Z'); 


e = X';
n = Y';
z = Z';

L = size(z',1);

%========================================================================================
% Polarization/ Azimuth /Incidence
%========================================================================================

nsall = size(z',1);   %  total number of samples to analyze 
nswin = fix(twin/delt);  %  number of samples in a time window 
% npshift = fix(twin/(2*delt))+1;  % number of samples to shift over 
npshift = 1;     % number of samples to shift over (recommend)
kfin = fix((nsall-nswin)/npshift); % number of time windows considered 

for k=1:1:kfin
   nwinst = (k-1)*npshift+1;  % start of time window
   nwinend = nwinst+nswin-1;    % end of time window
   ptrace_cut =  csigm(e,n,z,nwinst,nwinend);  % signal matrix %Load a time portion of the signals into the signal matrix.
   
   [ptrace_Cen, mu] = centerRows(ptrace_cut);
%    %traces(floor(tpb(nt)/dt):floor(tpe(nt)/dt),(nt-1)*nc+1:(nt-1)*nc+nc);

   [coeff,score,latent] = pca(ptrace_Cen,'Algorithm','SVD');
   dirVect(:,k) = coeff(:,1);
   
   % calculate the azimuth and incidence angle
   az(k) = atan2d(dirVect(1,k)*sign(dirVect(3,k)),dirVect(2,k)*sign(dirVect(3,k)));
   
   inc(k)= atan2d(sqrt(dirVect(1,k)^2+dirVect(2,k)^2), dirVect(3,k));

   maxeig(k) = latent(1);
    
   d(1) = latent(1);
   d(2) = latent(2);
   d(3) = latent(3);
    Dlp(k) = 1- (d(2)+d(3))/d(1);  %rectilinearity, r, or degree of linear polarization,as for pure body waves.
                                  %In the case of an elliptical polarized Rayleigh wave, 
                                  %the  lamuda1 and lamuda2 are much greater than lamuda3, 
                                  %so their degree of rectilinearity is small, (L > 0.5).
    
    Dpp(k) = 1 - 2*d(3)/(d(1)+d(2));% The degree of planarity,
    
   % calculate the ellipticity
    ellip(k) = real(sqrt(d(2)/d(1)));

   % Check to ensure az is between 0 - 180 degrees

   if az(k) < 0
        az(k) = az(k)+360;
    elseif  az(k) > 360
        az(k) = az(k)-360;
   end
   
   
   if inc(k) < 0
        inc(k) = inc(k)+180;
    elseif  inc(k) > 180
        inc(k) = inc(k)-180;
   end

end

azi = az;
inci = inc;





