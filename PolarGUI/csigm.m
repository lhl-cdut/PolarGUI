function [a]=csigm(e,n,z,k1,k2) 
% 
% [a]=csigm(e,n,z,k1,k2) 
% 
% Construct the complex signal matrix consisting of the spectra 
% of windows of three component time series data. 
e1=e(k1:k2); 
n1=n(k1:k2); 
z1=z(k1:k2); 
a=[e1;n1;z1]; 
return;