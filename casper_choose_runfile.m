function runfile = casper_choose_runfile(folder, pattern)
% CASPER_CHOOSE_RUNFILE  List run .mat files (newest first) and prompt user.
% Usage:
%   f = casper_choose_runfile();                    % defaults to runs/, casper_run_*.mat
%   f = casper_choose_runfile('runs');
%   f = casper_choose_runfile('runs','casper_run_*.mat')

if nargin < 1 || isempty(folder),  folder = 'runs'; end
if nargin < 2 || isempty(pattern), pattern = 'casper_run_*.mat'; end

if ~exist(folder,'dir'), mkdir(folder); end

D = dir(fullfile(folder, pattern));    D = D(~[D.isdir]);
D2 = dir(fullfile('.',  'casper_run_*.mat')); D2 = D2(~[D2.isdir]);

allFiles = [D; D2];
if isempty(allFiles), error('No run files found in "%s" or current folder.', folder); end

% sort newest first
[~, idx] = sort([allFiles.datenum], 'descend');
allFiles = allFiles(idx);

fprintf('\nAvailable run files (newest first):\n');
for i = 1:numel(allFiles)
    fprintf('  %2d) %s  [%s]\n', i, allFiles(i).name, datestr(allFiles(i).datenum,'yyyy-mm-dd HH:MM:SS'));
end

sel = input(sprintf('Select file [1-%d] (default=1): ', numel(allFiles)), 's');
if isempty(sel), sel = '1'; end
k = str2double(sel);
if isnan(k) || k < 1 || k > numel(allFiles), error('Invalid selection.'); end

runfile = fullfile(allFiles(k).folder, allFiles(k).name);
fprintf('➡️  Using: %s\n\n', runfile);
end
