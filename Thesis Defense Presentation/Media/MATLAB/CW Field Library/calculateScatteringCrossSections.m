function sigmas = calculateScatteringCrossSections(lambda, r, n)
	% References:
	%	https://en.wikipedia.org/wiki/Mie_scattering#/media/File:Radar_cross_section_of_metal_sphere_from_Mie_theory.svg
	%	https://en.wikipedia.org/wiki/Rayleigh_scattering

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