#takes event files created by guvTrack.py
function bob = gasTransportFitFunc(eventNum1, eventNum2, eventNum3, thresh)

%% Fit: 'untitled fit 1'.
%timeStep = .0075
timeStep = .01;


%[TMFile,TMPath]=uigetfile('*.csv', 'Select Intesity File...');
%filename1 = [TMPath,TMFile];

%[TMFile,TMPath]=uigetfile('*.csv', 'Select Intesity File2...');
%filename2 = [TMPath,TMFile];

%[TMFile,TMPath]=uigetfile('*.csv', 'Select Intesity File2...');
%filename3 = [TMPath,TMFile];


filename1 = ['channel1Events/channel1_event'  num2str(eventNum1) '.csv'];
filename2 =0;
filename3 = 0;
if nargin > 1
    filename2 = ['channel2Events/channel2_event'  num2str(eventNum2) '.csv']
end
if nargin > 2
    filename3 = ['channel3Events/channel3_event'  num2str(eventNum3) '.csv']
end

data = readtable(filename1);
x = data{:,1};
y = data{:,2};
if size(data, 2) > 2
    r = data{:,3};
end
event1 = regexp(filename1, '\d+(?=\.csv)', 'match', 'once')

if filename2 > 0
    data = readtable(filename2);
    x = [x; data{:,1}];
    y = [y; data{:,2}];
    if size(data, 2) > 2
        r = [r; data{:,3}];
    end
    event2 = regexp(filename2, '\d+(?=\.csv)', 'match', 'once')
end

if filename3 > 0
    data = readtable(filename3);
    x = [x; data{:,1}];
    y = [y; data{:,2}];
    if size(data, 2) > 2
        r = [r; data{:,3}];
    end
    event3 = regexp(filename3, '\d+(?=\.csv)', 'match', 'once')
end

%chop off first x points if measuring too close to the curve
%x = x(10:end);
%y = y(10:end);

% normalize the time to t0 = 0, and convert to seconds
initFrame = x(1);
x = x-x(1);
x = x*timeStep

%normalize the intensity to I0 = 1
y = y/y(1)

[xData, yData] = prepareCurveData( x, y );

if thresh>0
    mask = yData>thresh;
    yData(mask==0) = [];
    xData(mask==0) = [];
end

%Run from here if re-editing data by hand

% Set up fittype and options.
ft = fittype( 'a*exp(-x/b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];

opts.StartPoint = [0.063547147818684 1 0.533028237317691];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

coeffs = coeffvalues(fitresult);
uncert = confint(fitresult);
offset = coeffs(3)
offsetUncert = 0.5*uncert(2,3) - 0.5*uncert(1,3)

prefactor = coeffs(1)
prefacUncert = 0.5*uncert(2,1) - 0.5*uncert(1,1)


tau = coeffs(2)
tauError = 0.5*uncert(2,2) - 0.5*uncert(1,2);

if size(data,2)>2
    radius = median(r);
    radiusError = std(r);
    %P = coeffs*radius;
    %Perror = ((radius*uncert)^2 + (coeffs*radiusError)^2)^(1/2)
end


%extract even numbers from filenames








% Plot fit with data.
figure1 = figure( 'Name', 'untitled fit 1' );
hold on;
h = plot( fitresult, xData, yData );
legend( h, 'I vs. t', 'fit', 'Location', 'NorthEast' );




if size(data, 2) > 2
    annotation(figure1,'textbox',...
        [0.7 0.75 0.3 0.06],...
        'String',{['t0 = ' num2str(tau) ' � ' num2str(tauError)],['r = ' num2str(radius) ' � ' num2str(radiusError)]},...
        'FontSize',12,...
        'FitBoxToText','off',...
        'EdgeColor','none');
else
    annotation(figure1,'textbox',...
        [0.7 0.75 0.3 0.06],...
        'String',['t0 = ' num2str(tau) ' � ' num2str(tauError)],...
        'FontSize',12,...
        'FitBoxToText','off',...
        'EdgeColor','none');
end

if filename3>0
    title(['Events ' num2str(event1) '-' num2str(event2) '-' num2str(event3)]);
elseif filename2>0
    title(['Events ' num2str(event1) '-' num2str(event2)]);
else
    title(['Event ' num2str(event1)])
end
% Label axes
% Create ylabel
ylabel('I (arb.)','FontSize',14);

% Create xlabel
xlabel('t (s)','FontSize',14);
grid on

bob  = [initFrame, radius, radiusError, tau, tauError, prefactor, prefacUncert, offset, offsetUncert, eventNum1, eventNum2, eventNum3]
%bob = [initFrame, coeffs, uncert,prefactor, prefacUncert, offset, offsetUncert, eventNum1, eventNum2, eventNum3]
prompt = 'Save figure (1 for yes)? ';
resp = input(prompt)

if resp == 1
    saveas(gcf,['A02_' num2str(eventNum1) '-' num2str(eventNum2) '-' num2str(eventNum3) '.fig'])
end
