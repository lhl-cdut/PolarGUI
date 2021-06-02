function Z = loadAudio(paths)
% Syntax:   Z = loadAudio(paths);

% Parse inputs
p = numel(paths);

% Load audio
audio = cell(1,p);
for i = 1:p
    path = paths{i};
    [~, ~, ext] = fileparts(path);
    if strcmpi(ext,'.mat')
        % Load from .mat file
        data = load(path);
        fields = fieldnames(data);
        audio{i} = data.(fields{1});
    else
        % Load from .wav file
        try
            % Old MATLAB versions
            audio{i} = wavread(path,'double');
        catch %#ok
            % New MATLAB versions
            audio{i} = audioread(path,'double');
        end
    end   
end

% Determine sample length
nSamples = cellfun(@numel,audio);
n = min(nSamples);

% Construct audio matrix
Z = zeros(p,n);
for i = 1:p
    ni     = nSamples(i);
    gapL   = floor(0.5 * (ni - n));
    gapR   = ceil( 0.5 * (ni - n));
    Z(i,:) = audio{i}((gapL + 1):(ni - gapR));
end

% Normalize audio
Z = normalizeAudio(Z);
