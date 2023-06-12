%% NOTE
% This script simulates DORT experiments in a homogeneous, isotropic, nondispersive medium.  It treats the sources as isotropic point sources
% with radiation of the form (e^(-j*beta*R))/R.  Note that the (e^(-j*beta*R))/R expression IS TECHNICALLY ONLY VALID FOR THE FARFIELD, at least for electromagnetics.
% I am unsure about how badly that approximation breaks down in the near-field and Fresnel region, but it is certainly less valid in
% those regions.




%% Initialization stuff
clear;

%% Settings
% l_synth = 12.8e-6;
l_synth = 10e-3;
l1 = (1 / sqrt((8.854e-12)*(4*pi*(1e-7)))) / (300e9);
l2 = l1 * l_synth / (l_synth - l1);
lambda0 = [l1 l2]; %[854.22e-9 854.38e-9];
% lambda0 = l_synth;
n_medium = 1;

numElements = 8;
elementSpacing = 0.5e-3;%0.0508/2;
txLocs = [elementSpacing*transpose((-(numElements-1)/2):1:(numElements-1)/2) zeros(numElements,1) 2.5e-3 + zeros(numElements,1)];

rxLocs = txLocs;
% rxLocs = txLocs + [0 0.05 0];
% % rxLocs = [txLocs(:,2) txLocs(:,1) txLocs(:,3)] + [-0.025 0.025 0];
% % rxLocs = txLocs(480+1:end-480, :);
% % rxLocs = [rxLocs(:,2) rxLocs(:,1) rxLocs(:,3)] + [-0.025 0.025 0];


% targLocs = [	-0.01 0.03 0;
% 				0.01 0.035 0;
% 				0 0.025 0;
% 		];
% targRadii = [2e-6 2e-6 2e-6];
% targScatteringAmplitudesUnadjusted = [-1 -1 -1];	% Like the reflection coefficients for the electric field (this variable does not describe proportion of power scattered)

targLocs = [0.0081*2	0.08	0;
			-0.0106*2	0.076	0;
		];
targRadii = [0.02e-3 0.02e-3];
targScatteringAmplitudesUnadjusted = [-1 -1 -1 -1 -1];	% Like the reflection coefficients for the electric field (this variable does not describe proportion of power scattered)

	% Adjusting size
		targRadii = targRadii(1:size(targLocs, 1));
		targScatteringAmplitudesUnadjusted = targScatteringAmplitudesUnadjusted(1:size(targLocs, 1));

addNoiseToTransferMatrices = false;
transferMatrixSNR_dB = -10;

doSyntheicWavelengthStuff = false;

%% Settings for plotting interference patterns
zCut = 0;
xCoords = -0.05:0.000125/4:0.05; %-0.005:0.000025:0.005; %-0.5:0.01:0.5;
yCoords = 0:0.000125/4:0.1; %0:0.000025:0.01;

txRxFocusingSelector = 1;
wavelengthNum = 1;

yCutoff = 0.0; % Prevents the plot scale from being skewed by the factor of 1/R that gets large when close to the antenna elements

normalizeFocusingVectorsToNumberOfElements = false;

autoscale_plot_caxis = true;
plot_in_dB = true;
if (normalizeFocusingVectorsToNumberOfElements)
	max_plot_dB = 20;
	plotDynamicRange_dB = 40;
	linearScaleLims = [0 1];
else
	max_plot_dB = 50;
	plotDynamicRange_dB = 30;
	linearScaleLims = [0 30];
end


%% Code
N_wavelengths = length(lambda0);
lambda = lambda0;

v0 = 1 / sqrt((8.854e-12)*(4*pi*(1e-7)));
f0 = v0 ./ lambda0;

% Calculating scattering coefficients based on scattering cross sections
% obtained ROUGHLY from Mie theory (or Rayleigh scattering if you set stuff
% small enough, but it is probably best to keep stuff in the Mie regime for
% now...?)
targRCS = zeros(N_wavelengths, size(targLocs, 1));
for k = 1:N_wavelengths
	if (any((2*pi*targRadii / lambda0(k)) < 0.1))
		error("Target radius too small---in Rayleigh scattering regime.")
	end
	targRCS(k,:) = calculateScatteringCrossSections(lambda0(k), targRadii, j*12312313); % Last argument (index of refraction) does not matter because the code throws an error for Rayleigh scattering sizes and Rayleigh is the only one that makes use of the index of refraction.
end
	% Off by a factor of sqrt(4*pi)???  Maybe???
	%	targScatteringAmplitudes = reshape(targScatteringAmplitudesUnadjusted, 1, []) .* sqrt(targRCS);
targScatteringAmplitudes = getReflectionFactor(targRCS, reshape(targScatteringAmplitudesUnadjusted, 1, []));

% targScatteringAmplitudes = targScatteringAmplitudesUnadjusted;

H = get_CW_transfer_matrices(lambda0, n_medium, txLocs, rxLocs, targLocs, targScatteringAmplitudes);

if (addNoiseToTransferMatrices)
	H_no_noise = H;
	H = awgn(H, transferMatrixSNR_dB, 'measured');
end

if 0
	% Doing something experimental
	% Seems to have a flaw in that multiplication results in cross terms
	% (think multiplying Green's functions of the form exp{-jkr}/R) which get
	% detected when doing the SVD.
	H = H(:,:,1) .* conj(H(:,:,2));
	lambda = lambda0(1) * lambda0(2) / abs(lambda0(1) - lambda0(2));
end

if doSyntheicWavelengthStuff
	% Doing SWH stuff
	lambda = lambda0(1) * lambda0(2) / abs(lambda0(1) - lambda0(2));
	lambdaSynthConstituentWavelengths = [lambda0(1) lambda0(2)];

	[U1,S1,V1] = get_SVD_matrices_from_transfer_matrices(H(:,:,1));
	[U2,S2,V2] = get_SVD_matrices_from_transfer_matrices(H(:,:,2));
	U = U1 .* conj(U2);
	S = sqrt(S1 .* S2);
	V = V1 .* conj(V2);
else
	% Not doing SWH stuff
	[U,S,V] = get_SVD_matrices_from_transfer_matrices(H);
end

%% Plotting interference patterns
numPlots = size(targLocs,1);

figure(1);
clf;
tiledlayout(ceil(numPlots / ceil(sqrt(numPlots))), ceil(sqrt(numPlots)), "TileSpacing", "compact");
	% tiledlayout(2, numPlots, "TileSpacing", "compact");
	% tiledlayout(size(targLocs,1),size(targLocs,1));

% g_asdf = zeros(length(yCoords), length(xCoords));

fieldPatterns = zeros(length(yCoords), length(xCoords), length(lambda), size(targLocs, 1));

if (doSyntheicWavelengthStuff)
	wavelengthNum = 1;
end

% for k = 1:size(targLocs,1)^2
for k = 1:numPlots
	if (txRxFocusingSelector == 1)
		focusingVector = V(:,k,:);
	elseif (txRxFocusingSelector == 2)
		focusingVector = conj(U(:,k,:));
	else
		error("'txRxFocusingSelector' variable should be 1 for 'Tx' and 2 for 'Rx'.")
	end
% 	focusingVector = exp(j * angle(focusingVector)) / sqrt(size(focusingVector,1)); % Does phase-only focusing

	if (normalizeFocusingVectorsToNumberOfElements)
		focusingVector = focusingVector / size(focusingVector, 1);
	end

% 	focusingVector = focusingVector .* exp(j*2*pi*0.5*(rand(size(focusingVector,1),1) - 0.5));
	
	[g, xGrid, yGrid] = get_field_patterns(lambda, n_medium, txLocs, rxLocs, txRxFocusingSelector, focusingVector, xCoords, yCoords, zCut);
% 	g_asdf = g_asdf + (1./abs(g));
	fieldPatterns(:,:,:,k) = g;

	nexttile;
	
	plot_handles_temp = plot_field_pattern(g(:,:,wavelengthNum), txLocs*1000, rxLocs*1000, targLocs*1000, xCoords*1000, yCoords*1000, yGrid, yCutoff, autoscale_plot_caxis, plot_in_dB, max_plot_dB, plotDynamicRange_dB, linearScaleLims);
	plot_handles_temp.Colorbar.Label.FontSize = 20;
	plot_handles_temp.Legend.FontSize = 24; % 12;
	plot_handles_temp.Colorbar.FontSize = 20;

	ax = gca;
	ax.FontSize = 14;
	xlabel('X-Axis Location (mm)', 'FontSize', 24);
	ylabel('Y-Axis Location (mm)', 'FontSize', 24);

	if (txRxFocusingSelector == 1)
		title(['Field Patterns - Right Singular Vector #' num2str(k)], 'FontSize', 24);
		plot_handles_temp.Legend.Location = 'northeast';
	else
		title(['Field Patterns - Left Singular Vector #' num2str(k)], 'FontSize', 24);
		plot_handles_temp.Legend.Location = 'southwest';
	end
end
% if (doSyntheicWavelengthStuff)
% 	sgtitle({['\lambda_1 = ' sprintf('%.2f', lambdaSynthConstituentWavelengths(1)*1e9) 'nm, \lambda_2 = ' sprintf('%.2f', lambdaSynthConstituentWavelengths(2)*1e9) 'nm'], ['\Lambda_{synth} = ' sprintf('%.2f', lambda*1e3) 'mm'], ''}, 'FontSize', 24, 'FontWeight', 'bold');
% else
% % 	sgtitle({['\lambda = ' sprintf('%.2f', lambda(wavelengthNum)*1e9) 'nm'], ''}, 'FontSize', 24, 'FontWeight', 'bold');
% 	sgtitle({['Background Medium: Free Space (n = 1) — Frequency: ' sprintf('%.2f', f0(wavelengthNum)/1e9) ' GHz']}, 'FontSize', 36, 'FontWeight', 'bold');
% end
colormap('turbo');


% if 0
% 	figure(2);
% 	clf;
% 	tiledlayout(ceil(numPlots / ceil(sqrt(numPlots))), ceil(sqrt(numPlots)));
% 	for k = 1:size(targLocs,1)
% 		nexttile;
% 		tempSurfData = abs(fieldPatterns(:,:,wavelengthNum,k)) .* (yGrid >= 0.0075); % .* (yGrid >= 0.02) .* (yGrid <= 0.04);
% 	% 	if (plot_in_dB)
% 	% 		tempSurfData = 20*log10(tempSurfData);
% 	% 	end
% 		surf(tempSurfData, 'EdgeAlpha', 0);
% 		axis tight;
% 		view(2);
% 	end
% end


% clf;
% plot_field_pattern(g_asdf, k, txLocs, rxLocs, targLocs, xCoords, yCoords, yGrid, yCutoff, autoscale_plot_caxis, plot_in_dB, max_plot_dB, plotDynamicRange_dB, linearScaleLims);

% clf;
% asdf = get_field_patterns(repelem(lambda, 30), n_medium, txLocs, rxLocs, txRxFocusingSelector, reshape(V(:,20+(1:30),:), [1920 1 30]), xCoords, yCoords, zCut);
% plot_field_pattern(1./sqrt(sum(abs(asdf) .^ 2, 3)), k, txLocs, rxLocs, targLocs, xCoords, yCoords, yGrid, yCutoff, autoscale_plot_caxis, plot_in_dB, max_plot_dB, plotDynamicRange_dB, linearScaleLims);

% clf;
% asdf = get_field_patterns(lambda, n_medium, txLocs, rxLocs, txRxFocusingSelector, V(:,4,:), xCoords, yCoords, zCut);
% plot_field_pattern(1./asdf(:,:,1), k, txLocs, rxLocs, targLocs, xCoords, yCoords, yGrid, yCutoff, autoscale_plot_caxis, true, max_plot_dB, plotDynamicRange_dB, linearScaleLims);

%% Switching back to figure 1
figure(1)

%% Empty space

















%% Functions
function H = get_CW_transfer_matrices(lambda0, n_medium, txLocs, rxLocs, targLocs, targScatteringAmplitudes)
	k_lambda = n_medium .* (2*pi ./ lambda0);

% 	numTargs = size(targLocs, 1);
% 	N_wavelengths = length(k_lambda);
% 	numTxElements = size(txLocs, 1);
% 	numRxElements = size(rxLocs, 1);	
% 	H = gpuArray(zeros(numRxElements, numTxElements, N_wavelengths));

	k_lambda_temp = gpuArray(reshape(k_lambda, [1 1 1 length(k_lambda)]));

	if (isrow(targScatteringAmplitudes) || iscolumn(targScatteringAmplitudes))
		targScatteringAmplitudesTemp = gpuArray(reshape(targScatteringAmplitudes, [1 1 length(targScatteringAmplitudes)]));
	else
		targScatteringAmplitudesTemp = transpose(targScatteringAmplitudes);
		targScatteringAmplitudesTemp = gpuArray(reshape(targScatteringAmplitudes, [1 1 size(targScatteringAmplitudesTemp, 1) size(targScatteringAmplitudesTemp, 2)]));
	end
	targLocsTemp = gpuArray(reshape(permute(targLocs, [2 1]), [1 size(targLocs, 2) size(targLocs, 1)]));
	d1 = gpuArray(permute(sqrt(sum((txLocs - targLocsTemp) .^ 2, 2)), [2 1 3]));
	d2 = gpuArray(sqrt(sum((rxLocs - targLocsTemp) .^ 2, 2)));
	
% 	% Adds random path length variations
% 	fprintf("asdf\n");
% 	pause;
% 	d1 = d1 + (20e-6)*(rand(size(d1)) - 0.5);

	temp1 = targScatteringAmplitudesTemp .* (exp(j*k_lambda_temp.*d1) ./ d1);
	temp2 = (exp(j*k_lambda_temp.*d2) ./ d2);
	H = gather(squeeze(sum(temp2 .* temp1, 3)));


	if 0
		% For debugging.  Comparing H with H2 in this, it looks like the
		% vectorized GPU code gives the same results.
		H2 = zeros(size(rxLocs, 1), size(txLocs, 1), length(k_lambda));
		for m = 1:size(txLocs, 1)
			for k = 1:size(rxLocs, 1)
				H2_km = zeros(1,length(k_lambda));
				for l = 1:size(targLocs, 1)
					targLoc = targLocs(l,:);
					
					d1 = sqrt(sum((txLocs(m,:) - targLoc) .^ 2));
					d2 = sqrt(sum((rxLocs(k,:) - targLoc) .^ 2));
					d = d1 + d2;
		
					H2_km_temp = (exp(j*k_lambda*d2) / d2) .* targScatteringAmplitudes(l) .* (exp(j*k_lambda*d1) / d1);
					H2_km = H2_km + H2_km_temp;
				end
				H2(k,m,:) = H2_km;
			end
		end
	end


	% Old code, does not seem to separate transmitters and receivers:
% 	for k = 1:numElements
% 		for m = 1:numElements
% 			H_km = zeros(1,length(k_lambda));
% 			for l = 1:numTargs
% 				targLoc = targLocs(l,:);
% 				
% 				d1 = sqrt(sum((elementLocs(k,:) - targLoc) .^ 2));
% 				d2 = sqrt(sum((elementLocs(m,:) - targLoc) .^ 2));
% 				d = d1 + d2;
% 	
% 				H_km_temp = (exp(j*k_lambda*d2) / d2) .* targScatteringAmplitudes(l) .* (exp(j*k_lambda*d1) / d1);
% 				H_km = H_km + H_km_temp;
% 			end
% 			H(k,m,:) = H_km;
% 		end
% 	end
end


function [U,S,V] = get_SVD_matrices_from_transfer_matrices(H)
	N_wavelengths = size(H, 3);
	M = size(H, 1);
	N = size(H, 2);

	U = zeros(M, M, N_wavelengths);
	S = zeros(M, N, N_wavelengths);
	V = zeros(N, N, N_wavelengths);

	for k = 1:N_wavelengths
		[u,s,v] = svd(H(:,:,k));
		U(:,:,k) = u;
		S(:,:,k) = s;
		V(:,:,k) = v;
	end
end


function scaleFactor = getReflectionFactor(RCS, targScatteringAmplitudesUnadjusted)
	% IMPORTANT NOTE: I make no guarantees that this is correct.
	%
	% 'unadjustedScatteringAmplitude' is going to be treated like a reflection
	% coefficient for the field, i.e. outgoing field is proportional to
	% 'unadjustedScatteringAmplitude'.  It is not being treated as a constant
	% of proportionality for reflected energy---it is being treated as a constant of
	% proportionality for reflected for field magnitude.
	%
	% Incident energy = RCS * (E^2)/eta		(where RCS is radar cross section)
	% Reflected energy = RCS * ((C*E)^2)/eta	(where C = some number in targScatteringAmplitudesUnadjusted)
	% Assume exp(...)/R Green's functions
	% Scattered field: A * exp(...)/R	(assume isotropic)
	% Energy density at R is |A * exp(...)/R|^2 / eta = A^2 / R^2 / eta
	% At R, total energy is 4*pi*(R^2)*energy density = 4*pi*(A^2)/(R^2)/eta
	% Letting R = 1, one gets 4*pi*(A^2)/eta
	% Set equal to reflected energy
	% 4*pi*(A^2)/eta = RCS * ((C*E)^2)/eta
	% Then...
	% A/E = C * sqrt(RCS / (4*pi))

	scaleFactor = targScatteringAmplitudesUnadjusted .* sqrt(RCS / (4*pi));
end


function sigmas = calculateScatteringCrossSections(lambda, r, n)
	% References:
	%	https://en.wikipedia.org/wiki/Mie_scattering#/media/File:Radar_cross_section_of_metal_sphere_from_Mie_theory.svg
	%	https://en.wikipedia.org/wiki/Rayleigh_scattering

	% Test code:
		% aaar = (0.0001:0.0001:3)*lambda0(1)/(2*pi);
		% aaa = calculateScatteringCrossSections(lambda0(1), aaar, 1.0003); % 1.003 assumes air.  If 1 is used, smaller radii will result in zero.
		% loglog(2*pi*aaar/lambda0(1), aaa ./ (pi*(aaar.^2)));

	areas = pi*(r.^2);
	relFreq = 2*pi*r / lambda;
	sigmas = ones(1, length(r));
	
	sigmas(relFreq > 0.6) = areas(relFreq > 0.6);
	sigmas((relFreq <= 0.6) & (relFreq >= 0.1)) = (7.166244588 * (relFreq((relFreq <= 0.6) & (relFreq >= 0.1)) .^ 3.855291627)) .* areas((relFreq <= 0.6) & (relFreq >= 0.1));

	lowerInterpLim = 5e-3;
	upperInterpLim = 1e-1;	
	interpFactor = zeros(1, length(relFreq(relFreq < 0.1)));
	interpFactor((relFreq >= lowerInterpLim) & (relFreq < upperInterpLim)) = (relFreq((relFreq >= lowerInterpLim) & (relFreq < upperInterpLim)) - lowerInterpLim) / (upperInterpLim - lowerInterpLim);
	
	interpFactor = interpFactor .^ 4;
	interpFactor = interpFactor - min(interpFactor);
	interpFactor = interpFactor / max(interpFactor);
	
	temp1 = (7.166244588 * (relFreq(relFreq < 0.1) .^ 3.855291627)) .* areas(relFreq < 0.1);
	temp2 = (2*(pi^5)/3) * (((2*r(relFreq < 0.1)).^6) / (lambda^4)) * ((((n^2) - 1) / ((n^2) + 2)) ^ 2);
	temp3 = (temp1 .* interpFactor) + (temp2 .* (1 - interpFactor));
	sigmas(relFreq < 0.1) = temp3;
end


function [g, xGrid, yGrid] = get_field_patterns(lambda, n_medium, txLocs, rxLocs, txRxFocusingSelector, focusingVector, xCoords, yCoords, zCut)
	errFlagFocusingVector = false;
	if (length(size(focusingVector)) == 2)
		if (~iscolumn(focusingVector))
			errFlagFocusingVector = true;
		end
	elseif (length(size(focusingVector)) == 3)
		if (size(focusingVector, 2) ~= 1)
			errFlagFocusingVector = true;
		end
	else
		errFlagFocusingVector = true;
	end
	if (errFlagFocusingVector)
		error("'focusingVector' should be a column vector or an Nx1xM vector with the last dimension corresponding to wavelength.");
	end

	if (~isvector(lambda))
		error("'lambda0' should be a vector or a scalar.")
	end
	if (~isvector(n_medium))
		error("'n_medium' should be a vector or a scalar.")
	end

	if (txRxFocusingSelector == 1)
		elementLocs = txLocs;
	elseif (txRxFocusingSelector == 2)
		elementLocs = rxLocs;
	else
		error("'txRxFocusingSelector' argument should be 1 for 'Tx' and 2 for 'Rx'.")
	end

	[xGrid, yGrid] = meshgrid(xCoords, yCoords);
	g = zeros([size(xGrid) length(lambda)]);

	k_lambda = n_medium .* (2*pi ./ lambda);
	k_lambda = gpuArray(reshape(k_lambda, [1 1 1 length(k_lambda)]));
	
	dx = gpuArray(xGrid - reshape(elementLocs(:,1), [1 1 size(elementLocs, 1)]));
	dy = gpuArray(yGrid - reshape(elementLocs(:,2), [1 1 size(elementLocs, 1)]));
	dz = gpuArray(zCut - reshape(elementLocs(:,3), [1 1 size(elementLocs, 1)]));
	R = sqrt((dx.^2) + (dy.^2) + (dz.^2));

	clear("dx", "dy", "dz");
	wait(gpuDevice);

	for m = 1:length(k_lambda)
		focusingVecTemp = gpuArray(reshape(focusingVector(:,1,m), [1 1 size(focusingVector, 1) 1]));
		g(:,:,m) = sum(focusingVecTemp .* (exp(j*k_lambda(1,1,1,m).*R) ./ R), 3);
	end


% 	numElements = size(elementLocs, 1);

% 	for k = 1:numElements
% 		elemLoc = elementLocs(k,:);
% 		dx = xGrid - elemLoc(1);
% 		dy = yGrid - elemLoc(2);
% 		dz = zCut - elemLoc(3);
% 		
% 		R = sqrt((dx.^2) + (dy.^2) + (dz^2));
% 		g = g + (focusingVector(k,1,:) .* (exp(j*k_lambda.*R) ./ R));
% 	end
	
% 	scalarFieldsAtTargets = zeros(1, numTargs);
% 	for l = 1:numTargs
% 		for k = 1:numElements
% 			elemLoc = elementLocs(k,:);
% 			dx = targLocs(l,1) - elemLoc(1);
% 			dy = targLocs(l,2) - elemLoc(2);
% 			dz = targLocs(l,3) - elemLoc(3);
% 	
% 			R = sqrt((dx^2) + (dy^2) + (dz^2));
% 			scalarFieldsAtTargets(l) = scalarFieldsAtTargets(l) + focusingVector(k) * (exp(j*beta_fc*R) / R);
% 		end
% 	end
end


function handles = plot_field_pattern(g, txLocs, rxLocs, targLocs, xCoords, yCoords, yGrid, yCutoff, autoscale_plot_caxis, plot_in_dB, max_plot_dB, plotDynamicRange_dB, linearScaleLims)
	g_plot = g;
	g_plot = abs(g_plot);
	if (plot_in_dB)
		g_plot = 20*log10(g_plot);
% 		if (autoscale_plot_caxis)
% 			g_plot = min(g_plot, max_plot_dB);
% 			g_plot = max(g_plot, max(max(g_plot)) - plotDynamicRange_dB);
% 		end
	end
	g_plot(yGrid <= yCutoff) = min(min(g_plot));
	
	hold on;
	im = imagesc([min(xCoords) max(xCoords)], [min(yCoords) max(yCoords)], g_plot);
	set(gca, 'YDir', 'normal');
	c = colorbar;
	
% 	colormap('turbo');
	
	if (plot_in_dB)
		if (autoscale_plot_caxis)
			caxis([(max_plot_dB - plotDynamicRange_dB) max_plot_dB]);
		end
		c.Label.String = 'Field Magnitude (dB)';
	else
		if (autoscale_plot_caxis)
			caxis(linearScaleLims);
		end
		c.Label.String = 'Field Magnitude (Linear)';
	end
	
	scatterPlots = cell(1, 4);
	scatterPlots{1} = scatter(txLocs(:,1), txLocs(:,2), 100, 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'blue', 'MarkerEdgeAlpha', 1, 'DisplayName', 'Transducer Element');
	scatterPlots{3} = scatter(targLocs(:,1), targLocs(:,2), 350, 'd', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1, 'LineWidth', 2, 'DisplayName', 'Scatterer');
	scatterPlots = {{scatterPlots}};
	
	axis equal;
	xlim([min(xCoords) max(xCoords)]);
	ylim([min(yCoords) max(yCoords)]);
	
	lgd = legend('FontSize', 8, 'Location', 'best');
	
	hold off;

	handles = struct('Image', im, 'Legend', lgd, 'Colorbar', c, 'ScatterPlots', scatterPlots);
end