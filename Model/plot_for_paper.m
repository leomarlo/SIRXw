%% the script for everything

% plot crit transition
disp('plot at rate')
clear('all')
plot_SIRXw_atrate

% plot rinf at rate
disp('plot crit transitions')
clear('all')
plot_SIRXw_2d

% plot imax in 2d
disp('plot imax')
clear('all')
plot_imax_2d_mf_sim

% plot the implicit function
disp('plot implicit function')
clear('all')
plot_implicit_rho_vs_simulation

% plot positivity plot
disp('plot positivity')
clear('all')
plot_positivity

% plot imax examples
disp('plot imax examples')
clear('all')
% plot_sample_paths