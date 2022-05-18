function error_plot(err,text,alpha)
% function error_plot(err,text,alpha)
% Input - err: vector containing errors
%		  text: caption for legend
%		  alpha: current damping factor

semilogy(err,'DisplayName',text);

legend

xlabel ('Iterazioni','interpreter','latex');
ylabel ('Norma del residuo','interpreter','latex');

title (['$$\alpha$$ = ', num2str(alpha)],'interpreter','latex');

ax=gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
ax.XLabel.FontSize = 21;
ax.YLabel.FontSize = 21;
ax.Title.FontSize = 22;
ax.Legend.FontSize = 18;