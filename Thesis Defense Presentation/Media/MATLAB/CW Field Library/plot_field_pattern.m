function handles = plot_field_pattern(field, xGrid, yGrid, varargin)
	p = inputParser;
	addRequired(p, 'field', @(x) isnumeric(x) && ismatrix(x));
	addRequired(p, 'xGrid', @(x) isnumeric(x) && ismatrix(x));
	addRequired(p, 'yGrid', @(x) isnumeric(x) && ismatrix(x));
	addParameter(p, 'plot_scale', 'linear', @(x) any(strcmp(x, {'linear','dB'})));
	addParameter(p, 'plot_type', 'magnitude', @(x) any(strcmp(x, {'magnitude','phase'})));
	addParameter(p, 'autoscale_caxis', false, @(x) isscalar(x) && islogical(x))
	addParameter(p, 'magnitude_caxis_mask', true(size(field)), @(x) islogical(x) && isequal(size(x), size(field)));
	parse(p, field, xGrid, yGrid, varargin{:});

	plot_in_dB = strcmp(p.Results.plot_scale, 'dB');
	plot_phase = strcmp(p.Results.plot_type, 'phase');
	autoscale_caxis = p.Results.autoscale_caxis;
	magnitude_caxis_mask = p.Results.magnitude_caxis_mask;

	g_plot = field;
	if ~plot_phase
		g_plot = abs(g_plot);
		if (plot_in_dB)
			g_plot = 20*log10(g_plot);
		end
	else
		g_plot = angle(g_plot);
	end
	
	xlims = [min(min(xGrid)) max(max(xGrid))];
	ylims = [min(min(yGrid)) max(max(yGrid))];

	hold on;
	im = imagesc(xlims, ylims, g_plot);
	set(gca, 'YDir', 'normal');
	c = colorbar;
	
	if (plot_in_dB)
		c.Label.String = 'Field Magnitude (dB)';
	else
		c.Label.String = 'Field Magnitude (Linear)';
	end

	c.Label.FontSize = 16;

	axis equal;
	xlim(xlims);
	ylim(ylims);

	if autoscale_caxis && ~plot_phase
		maxVal = max(max(g_plot(magnitude_caxis_mask)));
		minVal = floor(min(min(g_plot(magnitude_caxis_mask))));
		caxis([minVal maxVal]);
	end

	hold off;

	handles = struct('Image', im, 'Colorbar', c);
end