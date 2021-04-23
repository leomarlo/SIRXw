
clear('pl');
figure; 
pl(3) = errorbar(betarange,simulation_mean, simulation_sdev/sqrt(sisi(2)));
% pl(3) = plot(betarange,simulations);
pl(3).Color = 'blue';
pl(3).LineWidth = 1.75;
pl(3).LineStyle = 'none';
pl(3).Marker = 'o';
pl(3).MarkerSize = 5;
pl(3).DisplayName = 'Simulation';
hold on;
pl(1) = plot(betarange,mean_field);
pl(1).LineWidth = 2.75;
pl(1).Color = 'black';
pl(1).DisplayName = 'Pair approximation';
% hold on;
pl(2) = plot(betarange,estimated);
pl(2).Color = 'red';
pl(2).LineWidth = 2.75;
pl(2).LineStyle = '--';
pl(2).DisplayName = 'Implicit equation from initial cond.';
% 
% pl(2).DisplayName = 'Implicit equation from \rho_S(0) and \rho_{SI}(0)';

ax = gca;
ax.YLabel.String = 'r_\infty';
ax.XLabel.String = '\beta';
ax.YLim = [0,0.12];
ax.XLim = [0,0.00252];
ax.FontSize = fontsize;
ax.LineWidth = FrameThickness;
ax.TickLength=TickLengths;
ax.FontWeight = FontWeights;

% daspect([1 1 1])

lgd = legend(pl);
lgd.Location = 'northwest';
% lgd.BoxFace.ColorType='truecoloralpha';
lgd.FontSize=LegendFontsize;

% grid on
% daspect([1 1 1])


resolution=600;
folder='figures/';
filename=strcat('Rinf_comparison_', date);
direction=strcat(folder,filename,'.png');
saveas(gcf,direction)