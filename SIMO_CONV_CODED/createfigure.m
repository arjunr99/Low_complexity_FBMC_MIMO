function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMATRIX1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 21-Sep-2017 13:52:28

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to semilogy
semilogy1 = semilogy(X1,YMatrix1,'LineWidth',2);
set(semilogy1(1),'DisplayName','OFDM-100-ETU','LineStyle','--');
set(semilogy1(2),'DisplayName','FBMC-100-ETU','Marker','*');

% Create xlabel
xlabel('SNR');

% Create ylabel
ylabel('64-QAM symbol error rate');

%% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.04 1]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14,'YMinorTick','on','YScale','log');
% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.158488444988161 0.268960797918625 0.340101515530733 0.121212118016534],...
    'FontSize',14);

