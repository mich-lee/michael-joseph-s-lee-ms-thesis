%% Initialization stuff
clear;
addpath('Colormaps/');
addpath('CW Field Library/');

%% Settings
lambda = 0.075;
n_medium = 1;

numElements = 4;
elementSpacing = 0.2;
elemYLoc = 0.1;
txLocs = [elementSpacing*transpose((-(numElements-1)/2):1:(numElements-1)/2) (zeros(numElements,1) + elemYLoc) zeros(numElements,1)];

rxLocs = zeros(size(txLocs));
rxLocs(:,2,:) = txLocs(:,1,:) + 0.5 + elementSpacing/2;
rxLocs(:,1,:) = 0.5 - elemYLoc;
rxLocs = flipud(rxLocs);

targLocs = [	0.15 0.5;
				-0.25 0.7;
		];
targRadii = [10e-3 10e-3];
targScatteringAmplitudesUnadjusted = [-1 -1] * 100;	% Like the reflection coefficients for the electric field (this variable does not describe proportion of power scattered)

zCut = 0;
deltaXY = 0.001;
xCoords = -0.5:deltaXY:0.5;
yCoords = 0:deltaXY:1;


%% Code
[xGrid, yGrid] = meshgrid(xCoords, yCoords);

targRCS = zeros(1, size(targLocs, 1));
targRCS(1,:) = calculateScatteringCrossSections(lambda, targRadii, j*12312313);
targScatteringAmplitudes = transpose(reshape(targScatteringAmplitudesUnadjusted, 1, []) .* sqrt(targRCS) / (4*pi));

if exist('v')
	close(v);
end
v = VideoWriter('newfile.avi','Motion JPEG AVI');
v.Quality = 75;
v.FrameRate = 15;
open(v);

fig = figure(1);
clf;
% outputFigurePosition = [100 50];
% outputImageSize = [1600 900];
% fig.Position = [outputFigurePosition outputImageSize];

for inputInd = 1:size(txLocs, 1)
	inputVector = zeros(size(txLocs, 1), 1);
	inputVector(inputInd) = 1;
	
	u_tx = get_cw_field_patterns(lambda, xGrid, yGrid, txLocs, inputVector, 'UseObliquity', true, 'apertureNormal', [0 1 0]);
	u_scatterers = get_cw_field_patterns(lambda, targLocs(:,1), targLocs(:,2), txLocs, inputVector, 'UseObliquity', true, 'apertureNormal', [0 1 0]);
	u_scattered = get_cw_field_patterns(lambda, xGrid, yGrid, [targLocs zeros(size(targLocs,1),1)], targScatteringAmplitudes .* u_scatterers, 'UseObliquity', false);

	u_rx = get_cw_field_patterns(lambda, rxLocs(:,1), rxLocs(:,2), [targLocs zeros(size(targLocs,1),1)], targScatteringAmplitudes .* u_scatterers, 'UseObliquity', false);
	u_rx = u_rx + get_cw_field_patterns(lambda, rxLocs(:,1), rxLocs(:,2), txLocs, inputVector, 'UseObliquity', true, 'apertureNormal', [0 1 0]);
	
	u_total = u_tx + u_scattered;
	
	% clf;
	% % colormap('turbo');
	% plot_field_pattern(u_total, xGrid, yGrid, 'plot_scale', 'linear', 'plot_type', 'magnitude');
	% caxis([0 10]);

	plotTitle = ['Measuring response from Tx' num2str(inputInd) '\rightarrow \{'];
	for ind0=1:size(rxLocs,1)
		plotTitle = [plotTitle 'Rx' num2str(ind0) ', '];
	end
	plotTitle = [plotTitle(1:end-2) '\}'];
	
	framesPerPeriod = 15;
	txEnvelopeSchedule = [1 5 45 5 1];
	scatteredEnvelopeSchedule = txEnvelopeSchedule;
	txEnvelope = getFadeInFadeOutEnvelope(txEnvelopeSchedule);
	scatteredEnvelope = getFadeInFadeOutEnvelope(scatteredEnvelopeSchedule);
	
	if length(txEnvelope) ~= length(scatteredEnvelope)
		error("asdf");
	end
	
	u_total_td = get_time_domain_fields(u_total, (0:length(txEnvelope)-1)*((1/3.99723277e9)/framesPerPeriod), lambda);
	u_tx_td = get_time_domain_fields(u_tx, (0:length(txEnvelope)-1)*((1/3.99723277e9)/framesPerPeriod), lambda);
	u_scattered_td = get_time_domain_fields(u_scattered, (0:length(txEnvelope)-1)*((1/3.99723277e9)/framesPerPeriod), lambda);
	u_rx_td_phasor = get_time_domain_fields(u_rx, (0:length(txEnvelope)-1)*((1/3.99723277e9)/framesPerPeriod), lambda, 'return_complex', true);
	u_rx_td_phasor_normed = u_rx_td_phasor ./ sqrt(sum(abs(u_rx_td_phasor).^2, 1));
	delta_phi = 2*pi/framesPerPeriod;

	for ind=1:length(txEnvelope)
		u_temp = txEnvelope(ind) * u_tx_td(:,:,ind) + scatteredEnvelope(ind) * u_scattered_td(:,:,ind);
		clf;
		subplot(1, 3, [1 2]);
		handles1 = plot_field_pattern(u_temp, xGrid, yGrid, 'plot_scale', 'linear', 'plot_type', 'magnitude');
		handles2 = plot_tx_rx_scatterers(txLocs, rxLocs, targLocs);
		caxis([0 20]);
		delete(handles1.Colorbar);
		handles2.ScatterPlots{1}{2}.MarkerEdgeColor = [0 1 0];
		title(plotTitle, 'FontSize', 22);
		set(gca, 'xtick', []);
		set(gca, 'ytick', []);

		numSideplots = size(rxLocs,1) + 1;
		phasorPlotHandles = [];
		for ind2=1:size(rxLocs,1)
			subplot(numSideplots, 3, ind2*3);
			tempAngle = angle(u_rx_td_phasor_normed(ind2,1,ind));
			tempMag = abs(u_rx_td_phasor_normed(ind2,1,ind)) * scatteredEnvelope(ind);
			polarplot([tempAngle tempAngle], [0 tempMag], 'LineWidth', 2);
			rlim([0 1]);
			title({'Field Phasor', ['Receiver #' num2str(ind2)]}, 'FontSize', 12);
			phasorPlotHandles = [phasorPlotHandles gca];
		end

		subplot(numSideplots, 3, numSideplots*3);
		tempAngle = (ind - 1) * delta_phi; % Assumes that the transmitter starts with phase 0 at t = 0
		tempMag = 1 * txEnvelope(ind);
		polarplot([tempAngle tempAngle], [0 tempMag], 'LineWidth', 2, 'Color', 'red');
		rlim([0 1]);
		title({'Field Phasor', ['Transmitter #' num2str(inputInd)]}, 'FontSize', 12);
		phasorPlotHandles = [phasorPlotHandles gca];

		for ind2=1:length(phasorPlotHandles)
			tempHandle = phasorPlotHandles(ind2);
			tempHandle.RTickLabel(:) = '';
			tempHandle.ThetaTickLabel(:) = '';
			set(tempHandle, 'OuterPosition', [0.51 tempHandle.Position(2:end)]);
		end

		drawnow;
		frame = getframe(gcf);
    	writeVideo(v,frame);
	end
end

close(v);