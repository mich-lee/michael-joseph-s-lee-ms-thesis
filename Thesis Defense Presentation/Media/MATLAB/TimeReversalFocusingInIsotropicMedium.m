%% Initialization stuff
clear;
addpath('Colormaps/');

%% Constants and settings
epsilon = 8.854e-12;
mu = 4*pi*(1e-7);

v_prop = 1 / sqrt(epsilon*mu);

Fs = 50e9;
N = 2000;
BW = [3e9 6e9];

freqsAll = (0:(N-1))*(Fs/N);
% freqs = freqsAll((0:301) + 61); % freqsAll((0:101) + 61);
freqs = freqsAll((freqsAll >= BW(1)) & (freqsAll <= BW(2)));
timeDomainPeriod = N / Fs;

transmitFromTargets = true;

trackMaximumFieldMagnitude = false;
maximumTrackerMinimumY = 0.3;

w = 1;
xmin = -w/2;
xmax = w/2;
ymin = 0;
ymax = w;

d = 0.001; % 0.005;
d = min(d, v_prop / max(freqs) / 2); % Forces the resolution to be at least as high as half the smallest wavelength
xdelta = d;
ydelta = d;

numelems = 4*2;
elemspacing = 0.25/2;
elemYLoc = 0.05;

targLocs = [	0.15 0.5;
				-0.25 0.7;
							]; % [-0.25 0.7];
targNums = [1 2];

% tStart_TargetTx = -0.5e-9;
% tStop_TargetTx = tStart_TargetTx + 0.2*timeDomainPeriod;%3.155906406071504e-09;
% tStart_ElementTx = -4e-9;%-3.155906406071504e-09;
% tStop_ElementTx = tStart_ElementTx + 8e-9;

tStart = -0.5e-9;
tStop = tStart + 18e-9;

% Time when the last peak reaches an element after the targets transmit
%	The value -tDurPeak should approximately be when the elements start
%	transmitting when doing non-causal time-reversal.
%	Delaying the non-causal time-reversed received signals by tDurPeak should
%	result in the retransmitted fields starting at around t=0.
tDurPeak = 3.155906406071504e-09;

% Should approximately be the time between end of receiving and beginning
%	of transmitting.
retransmitDelay = 6.8441e-09;

retransmitOffset = 2*tDurPeak + retransmitDelay;
elemTxStart = retransmitOffset - tDurPeak; % Impulse should occur around this time when retransmitting
tStep = (1/mean(freqs)) / 1;

noncausalityReducingEnvelopeOffset = -tStart/2;
noncausalityReducingEnvelopeSlope = -1 / (5 * tStart/2);
noncausalityReducingEnvelopeOneCutoff = 0.99;

useGPU = true;

display_title_in_plot = false;
outputFigurePosition = [100 50];
outputImageSize = [900 900];
zlvl = 2500;

targetVideoLenSeconds = 10;


%% Main code
xvals = xmin:xdelta:xmax;
yvals = ymin:ydelta:ymax;

[xgrid, ygrid] = meshgrid(xvals, yvals);
ygrid = flipud(ygrid);

targLocs = targLocs(targNums, :);
elemlocs = [transpose(elemspacing * ((0:(numelems-1)) - ((numelems - 1)/2))) (zeros(numelems, 1) + elemYLoc)];

% if (transmitFromTargets)
% 	tStart = tStart_TargetTx;
% 	tStop = tStop_TargetTx;
% else
% 	tStart = tStart_ElementTx;
% 	tStop = tStop_ElementTx;
% end

maxAmplitudeLoc = [0 0];
maxAmplitude = 0;

tdFromTargetToElements = zeros(length(targNums), numelems);
for k = 1:numelems
	for l = 1:length(targNums)
		tdFromTargetToElements(l, k) = sqrt((elemlocs(k,1) - targLocs(l,1))^2 + (elemlocs(k,2) - targLocs(l,2))^2) / v_prop;
	end
end

[g_freq_data_targ_tx, g_elem_to_targ_data] = getFreqDomainGreensFunctions(freqs, xgrid, ygrid, targLocs, targNums, numelems, elemlocs, mu, epsilon, useGPU, true);
[g_freq_data_elem_tx, ~] = getFreqDomainGreensFunctions(freqs, xgrid, ygrid, targLocs, targNums, numelems, elemlocs, mu, epsilon, useGPU, false);

retransmitDelayResponse = zeros(1, 1, length(freqs));
retransmitDelayResponse(1, 1, :) = exp(-j*2*pi*freqs*retransmitOffset);
g_freq_data_elem_tx = g_freq_data_elem_tx .* retransmitDelayResponse;


%% Computing time-domain stuff and displaying
fprintf("Computing and displaying time-domain fields...\n");

fig = figure(1);

tic;
numIter = 0;

g_freq_data = zeros(size(g_freq_data_elem_tx));
tVec = tStart:tStep:tStop;
targetTxEnv = getSigmoidEnvelope(tVec, 0, noncausalityReducingEnvelopeSlope, noncausalityReducingEnvelopeOffset, noncausalityReducingEnvelopeOneCutoff);
elemTxEnv = getSigmoidEnvelope(tVec, retransmitOffset, noncausalityReducingEnvelopeSlope, noncausalityReducingEnvelopeOffset, noncausalityReducingEnvelopeOneCutoff);

if useGPU
	g_freq_data_targ_tx = gpuArray(g_freq_data_targ_tx);
	g_freq_data_elem_tx = gpuArray(g_freq_data_elem_tx);
	g_freq_data = gpuArray(g_freq_data);
end

% cmap = transpose(127:-1:0); cmap = [cmap zeros(128, 1) flipud(cmap); zeros(128, 1) flipud(cmap) cmap]; cmap = cmap / 127; cmap = log(cmap + 1); cmap = cmap / max(max(cmap)); cmap = cmap .* [1 1 0.25];
% colormap(cmap);

curTimeOfDaySecondsTotal = mod(now,1)*(24*60*60);
curTimeOfDaySeconds = mod(floor(curTimeOfDaySecondsTotal), 60);
curTimeOfDayMinutes = mod(floor(curTimeOfDaySecondsTotal / 60), 60);
curTimeOfDayHours = mod(floor(curTimeOfDaySecondsTotal / 3600), 24);
fileOutputName = ['fields' num2str(curTimeOfDayHours) num2str(curTimeOfDayMinutes) num2str(curTimeOfDaySeconds) '.mp4'];
fileOutputFPS = length(tVec) / targetVideoLenSeconds;
if exist('v')
	close(v);
end
v = VideoWriter(fileOutputName,'MPEG-4');
v.Quality = 75;
v.FrameRate = fileOutputFPS;
open(v);

for tInd = 1:length(tVec)
	t = tVec(tInd);

	clf;
	hold on;

	fig.Position = [outputFigurePosition outputImageSize];
	axes('Units', 'normalized', 'Position', [0 0 1 1]);
	hold on

	fprintf("t = %.2f ns\t|", t*(1e9));

	g_freq_data(:,:,:) = (targetTxEnv(tInd) * g_freq_data_targ_tx) + (elemTxEnv(tInd) .* g_freq_data_elem_tx);

	f = zeros(1, 1, length(freqs));
	f(1,1,:) = freqs;
	omega = 2*pi*f;
	gField = 2*real(sum(g_freq_data .* exp(j*omega*t), 3));

	imagesc([xvals(1) xvals(end)], [yvals(end) yvals(1)], gField);
	set(gca, 'YDir', 'normal');
% 	colorbar;
	caxis([-zlvl zlvl]);

	scatter(elemlocs(:,1), elemlocs(:,2), 500, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 4, 'DisplayName', 'Tx/Rx Elements');
	for l = 1:length(targNums)
		scatter(targLocs(l,1), targLocs(l,2), 200, 'd', 'LineWidth', 3, 'DisplayName', ['Point Source #' num2str(l)]);
	end
	
	if (trackMaximumFieldMagnitude)
		gFieldAbs = abs(gField);
		gFieldAbs(ygrid < maximumTrackerMinimumY) = -Inf;
		gFieldMaxTemp = max(max(gFieldAbs));
		if (gFieldMaxTemp >= maxAmplitude)
			maxAmplitude = gFieldMaxTemp;
			gFieldMaxTempX = xgrid(gFieldAbs == gFieldMaxTemp);
			gFieldMaxTempY = ygrid(gFieldAbs == gFieldMaxTemp);
			maxAmplitudeLoc = [gFieldMaxTempX gFieldMaxTempY];
		end
		scatter(maxAmplitudeLoc(1), maxAmplitudeLoc(2), 'o', 'LineWidth', 2);
	end
	
	axis tight;
	
	xlim([xvals(1) xvals(end)]);
	ylim([yvals(1) yvals(end)]);
	
	xticks([]);
	yticks([]);
	
	if (display_title_in_plot)
		title(['t = ' num2str(t * (1e9)) 'ns']);
	end

	legend('Location', 'northeast', 'FontSize', 16);

	if (t < ((tDurPeak + elemTxStart)/2))
		text(-0.45, 0.9, {'Scatterers/point sources emit impulses.', 'Tx/Rx elements record received signals.'}, 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 24);
	else
		text(-0.45, 0.9, {'Tx/Rx elements transmit time-reversed','versions of received signals.'}, 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 24);
	end
	if (t >= 12e-9)
		text(-0.45, 0.2, {'Focusing on scatterers/point sources occurs.'}, 'HorizontalAlignment', 'left', 'FontWeight', 'bold', 'FontSize', 24);
	end
	
% 	drawnow limitrate;
	drawnow;
	hold off;
	
	frame = getframe(gcf);
	writeVideo(v,frame);

% 	commandwindow;
% 	pause;

	tempTime = toc;
	numIter = numIter + 1;
	avgTimePerIteration = tempTime / numIter;
	fprintf("\tApproximate rate: %f seconds/iteration\n", avgTimePerIteration);
end

close(v);




%% Functions
function [g_freq_data, elementsToTargetGreensFuncData] = getFreqDomainGreensFunctions(freqs, xgrid, ygrid, targLocs, targNums, numelems, elemlocs, mu, epsilon, useGPU, transmitFromTargets)
	% Precomputing some data
	fprintf("Computing R matrices...");
	if (transmitFromTargets)
		R_data = zeros([size(xgrid) length(targNums)]);
		for l = 1:length(targNums)
			r0 = targLocs(l,:);
			dx = xgrid - r0(1);
			dy = ygrid - r0(2);
			R_data(:,:,l) = sqrt((dx.^2) + (dy.^2));
		end
	else
		R_data = zeros([size(xgrid) numelems]);
		for k = 1:numelems
			r0 = elemlocs(k,:);
			dx = xgrid - r0(1);
			dy = ygrid - r0(2);
			R_temp = sqrt((dx.^2) + (dy.^2));
			R_data(:,:,k) = R_temp;
		end
	end
	if useGPU
		R_data = gpuArray(R_data);
	end
	fprintf("Done!\n");
	
	fprintf("Computing frequency domain fields...\n");
	fprintf("\t0%%...");
	elementsToTargetGreensFuncData = zeros(length(freqs), length(targNums), numelems);
	g_freq_data = zeros([size(xgrid) length(freqs)]);
	for fInd = 1:length(freqs)
		if (mod(fInd, 20) == 0)
			fprintf("%.2f%%...", (fInd - 1) / length(freqs) * 100);
		end
		
		f = freqs(fInd);
		
		beta = 2*pi*f*sqrt(mu*epsilon);
		omega = 2*pi*f;
	
		elemphases = zeros(1, numelems);
		for k = 1:numelems
			r0 = elemlocs(k,:);
			for l = 1:length(targNums)
				targLocTemp = targLocs(l,:);
				dx = targLocTemp(1) - r0(1);
				dy = targLocTemp(2) - r0(2);
				R = sqrt((dx^2) + (dy^2));
				elementsToTargetGreensFuncData(fInd,l,k) = exp(-j*beta*R) / R;
			end
		end
	
		g_temp = zeros(size(xgrid));
		if (transmitFromTargets)
			for l = 1:length(targNums)
				R = R_data(:,:,l);
				g_temp = g_temp + (exp(-j*beta*R) ./ R);
			end
		else
			for k = 1:numelems
				R = R_data(:,:,k);
				tx_sig_temp = 0;
				for l = 1:length(targNums)
					tx_sig_temp = tx_sig_temp + conj(elementsToTargetGreensFuncData(fInd,l,k));
				end
				g_temp = g_temp + (tx_sig_temp * (exp(-j*beta*R) ./ R));
			end
		end
		
		g_freq_data(:,:,fInd) = g_temp;
	end
	fprintf("100%%\nDone!\n");
	
	fprintf("\n");
end


function env = getSigmoidEnvelope(t, t0, alpha, t0_shift, oneCutoff)
	env = 1 ./ (1 + exp(-alpha * (t - t0 - t0_shift)));
	env(env >= oneCutoff) = 1;
end