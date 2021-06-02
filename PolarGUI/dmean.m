function [a]=dmean(b) 
% 
% [a]=dmean(b) 
% 
%     Remove the mean from a row vector 
m=mean(b); 
a=b-m; 
return;