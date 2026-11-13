figure
hold on
title("Net Velocity Changes")
group=repmat(["GSH","MET","Urea","FABs"],length(total_velocities_case1),1);
group=categorical(reshape(group,length(total_velocities_case1)*4,1));
data=reshape(total_velocities_case1,length(total_velocities_case1)*4,1);
c = zeros(length(data), 3); 
c(data > 0, :) = repmat([.25 0 .5], sum(data > 0), 1);
s = scatter(group, data, 15,c,'filled', 'jitter', 'on', 'JitterAmount', 0.15);
s.MarkerFaceAlpha = 0.5; % Set transparency
yline(0, 'r'); % Horizontal line at y=0
ylabel("net velocity change (mM/min)")
xlabel("module")
fontsize(gcf, "scale", 3.52);
hFigure = gcf;
hFigure.Color = [1,1,1];
hAxes = findobj(gcf,"Type","axes");
hAxes.FontSize = 18;
hConstantline = findobj(gcf,"Type","constantline");
hConstantline.LineWidth = 2;
hold off;





