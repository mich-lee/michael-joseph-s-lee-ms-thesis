function handles = plot_tx_rx_scatterers(txLocs, rxLocs, scattererLocs, varargin)
	hold on;
	
	scatterPlots = cell(1, 4);
	if ~isempty(txLocs)
		scatterPlots{1} = scatter(txLocs(:,1), txLocs(:,2), 650, 'x', 'MarkerEdgeColor', 'red', 'LineWidth', 4, 'DisplayName', 'Tx Element');
	end
	if ~isempty(rxLocs)
		scatterPlots{2} = scatter(rxLocs(:,1), rxLocs(:,2), 650, 'x', 'MarkerEdgeColor', 'blue', 'LineWidth', 4, 'DisplayName', 'Rx Element');
	end
	if ~isempty(scattererLocs)
		scatterPlots{3} = scatter(scattererLocs(:,1), scattererLocs(:,2), 500, 'd', 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'white', 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 1, 'LineWidth', 4, 'DisplayName', 'Scatterer', 'HandleVisibility', 'off');
		scatterPlots{4} = scatter(scattererLocs(:,1), scattererLocs(:,2), 500, 'd', 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'white', 'MarkerFaceAlpha', 0, 'MarkerEdgeAlpha', 1, 'LineWidth', 1, 'DisplayName', 'Scatterer');
	end
	scatterPlots = {{scatterPlots}};
	
	lgd = legend('FontSize', 16, 'FontWeight', 'bold', 'Location', 'northwest');
	
	hold off;

	handles = struct('Legend', lgd, 'ScatterPlots', scatterPlots);
end