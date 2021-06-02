function altwritesegy(sgyfile, datain, sampint, numsamps, numtraces, segfmt, byteorder)
% altwritesegy - a platform independent SEG-Y file writer
%
% This is a special version of writesegy meant for use with synth and logsec.
%
% segfmt should generally be 5, for IEEE floating point data generation
%
% function altwritesegy(sgyfile, datain, sampint, numsamps, numtraces, segfmt, byteorder)
%
%   sgyfile   - filename (should end in .sgy)
%   datain    - 2-D matrix of samples indexed by (sample, trace)
%   sampint   - sample interval in seconds
%   numsamps  - number of samples per trace (optional)
%   numtraces - number of traces (optional)
%   segfmt    - output segfmt (optional, default=5)
%               1=IBM floating point, 2=4-byte integer, 3=2-byte integer, 5=IEEE floating point
%   byteorder - output byte order (optional, default=be)
%               'be' = big-endian (standard), 'le' = little-endian (not conformant to standard,
%                                                    used by bad PC implementations)
%
% Example:
%  datain = ones(1024,10)   % 10 traces of 1024 samples
%  altwritesegy('file.sgy', datain, 0.002);
%
% BEGIN TERMS OF USE LICENSE
%
% This SOFTWARE is maintained by the CREWES Project at the Department
% of Geology and Geophysics of the University of Calgary, Calgary,
% Alberta, Canada.  The copyright and ownership is jointly held by
% its author (identified above) and the CREWES Project.  The CREWES 
% project may be contacted via email at:  crewesinfo@crewes.org
%
% The term 'SOFTWARE' refers to the Matlab source code, translations to
% any other computer language, or object code

% Terms of use of this SOFTWARE
%
% 1) Use of this SOFTWARE by any for-profit commercial organization is
%    expressly forbidden unless said organization is a CREWES Project
%    Sponsor.
%
% 2) A CREWES Project sponsor may use this SOFTWARE under the terms of the
%    CREWES Project Sponsorship agreement.
%
% 3) A student or employee of a non-profit educational institution may
%    use this SOFTWARE subject to the following terms and conditions:
%    - this SOFTWARE is for teaching or research purposes only.
%    - this SOFTWARE may be distributed to other students or researchers
%      provided that these license terms are included.
%    - reselling the SOFTWARE, or including it or any portion of it, in any
%      software that will be resold is expressly forbidden.
%    - transfering the SOFTWARE in any form to a commercial firm or any
%      other for-profit organization is expressly forbidden.
%
% END TERMS OF USE LICENSE

if nargin < 3 
   error('Invalid number of arguments: at least 3 required (segyfile, datain, sampint)');
end
if nargin < 4 
   numsamps = size(datain,1);
end
if nargin < 5 
   numtraces = size(datain,2);
end
if nargin < 6 
   segfmt = 5;
end
if nargin < 7 
   byteorder = 'be';
end

% Open or create the file for writing.

switch lower(byteorder)
case {'be', 'bigendian', 'ieee-be'}
   platform='ieee-be';
case {'le', 'littleendian', 'ieee-le'}
   platform = 'ieee-le';
otherwise
   error(['Invalid byte order parameter: ' byteorder]);
end


[fid, errmsg] = fopen(sgyfile, 'w', platform);

if (fid == -1) 
    error(['Unable to create ' sgyfile ': ' errmsg] );
end
% Build a text header
textLine{1} = char(['File generated ' datestr(now)]);

switch segfmt
case 1
  textLine{1} = 'IBM Floating point fmt (not implemented)';
  dformat = 'float32'; % 4 bytes, should be IBM floating point
case 2
  textLine{2} = '4-byte integer format';
  dformat = 'int32'; % 4 bytes, signed.
case 3
  textLine{3} = '2-byte integer format';
  dformat = 'int16'; % 2 bytes, signed.
case 4
  textLine{4} = 'WARNING BAD DATA: Fixed point w/gain code format (not implemented)'
  error('Can not write this format. (Fixed point with gain code.)');
case 5
  textLine{5} = 'IEEE Floating point format';
  dformat = 'float32'; % 4 bytes presumably IEEE floating point
case 8
  textLine{8} = '1-byte integer format';
  dformat = 'uchar8'; % 4 bytes presumably IEEE floating point
otherwise
  error(['invalid format specified: ', num2str(segfmt)]);
end

textLine{39} = 'SEG Y rev1';    % ooh, bleeding edge
textLine{40} = 'END EBCDIC';    %
% Next line should contain 80 spaces
textLine{41} = '                                                                                ';
for i = 1:40
   if isempty(textLine{i})
      textLine{i} = ' ';
   end
end
textHdrBlock = char(textLine);
textHdrBlock = [reshape(sprintf('C%2d ',1:40),4,40)' textHdrBlock(1:40,1:76)];
textHdrBlock = ascii2ebcdic(textHdrBlock);
% Write the text header.
count = fwrite(fid, textHdrBlock', 'uchar');

% Check if writing went successfully.
if count ~= 3200
    error(['EBCDIC header is too short. Size is only ', num2str(count), '.']);
end

% Now fill in the binary header.
temp = zeros(27,1);
temp(6) = sampint*1000000;
temp(8) = numsamps;
temp(10) = segfmt;

if (fwrite(fid, temp(1:3), 'ulong') ~= 3) | ...
    (fwrite(fid, temp(4:27), 'ushort') ~= 24) | ...
    (fwrite(fid, zeros(170,1), 'ushort') ~= 170) ;

    error('The binary header has been truncated.')
end

clear temp

% Write out the data traces and trace headers.
% Other initializations needed for the while loop.
tracewrit = 1;
trcount = 0;
traceheader = zeros(120,1);
traceheader(floor((29+1)/2)) = 1;   % seismic data 
traceheader(floor((35+1)/2)) = 1;   % production data
traceheader(floor((115+1)/2)) = numsamps;  % number of samples for this trace


while tracewrit <= numtraces
    trcount = fwrite(fid, traceheader, 'ushort');

    count = fwrite(fid, datain(:, tracewrit), dformat);
    if count ~= numsamps
       error(['Problem writing trace ', num2str(tracewrit), '.'])
    end
    tracewrit = tracewrit + 1;
end

fclose(fid);

