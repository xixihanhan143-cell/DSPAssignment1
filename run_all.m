%% Run all exercises of Assignment 1 (DFT spectral analysis)
% Exercise 2: frequency resolution, zero padding, window type
% Exercise 3: scale yueyin1.wav -> notes and key
% Exercise 4: chord hexian1.wav -> notes; 16th-note resolution problem
% Exercise 5: scale with harmonics yueyin2.wav -> harmonic composition
%
% Figures are saved to ./figures, the command-window output to results_log.txt
% (each exercise script starts with "clear", so no variables are kept here)
cd(fileparts(mfilename('fullpath')));
if exist('results_log.txt', 'file'), delete('results_log.txt'); end
diary('results_log.txt');
fprintf('\n########## Exercise 2 ##########\n');  run('ex2_resolution.m');
fprintf('\n########## Exercise 3 ##########\n');  run('ex3_scale.m');
fprintf('\n########## Exercise 4 ##########\n');  run('ex4_chord.m');
fprintf('\n########## Exercise 5 ##########\n');  run('ex5_harmonic.m');
diary off;
fprintf('\nAll done. Figures are in the "figures" folder, output in results_log.txt\n');
