function fieldsOut = get_time_domain_fields(field, times, lambda, varargin) % Lambda is in free space
	p = inputParser;
	addRequired(p, 'field', @(x) isnumeric(x) && ismatrix(x));
	addRequired(p, 'times', @(x) isnumeric(x) && isvector(x));
	addRequired(p, 'lambda', @isvector);
	addParameter(p, 'return_complex', false, @(x) isscalar(x) && islogical(x))
	parse(p, field, times, lambda, varargin{:});

	return_complex = p.Results.return_complex;
	
	c = 1/sqrt((8.854e-12)*(4*pi*1e-7));
	freq = c / lambda;

	t = zeros(1, 1, length(times));
	t(1,1,:) = times;

	fieldsOut = field .* exp(j*2*pi*freq*t);
	if ~return_complex
		fieldsOut = 2 * real(fieldsOut);
	end
end