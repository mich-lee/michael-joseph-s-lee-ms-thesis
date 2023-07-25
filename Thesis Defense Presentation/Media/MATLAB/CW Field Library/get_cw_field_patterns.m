function g = get_cw_field_patterns(lambda, xGrid, yGrid, elementLocs, inputVector, varargin)
	function output = gpuArrayHelper(x, gpuUseBool)
		if gpuUseBool
			output = gpuArray(x);
		else
			output = x;
		end
	end

	function ret = validateInputVector(inputVec)
		% 'inputVector' should be a column vector or an Nx1xM vector with the last dimension corresponding to wavelength.
		if (length(size(inputVec)) == 2)
			if (~iscolumn(inputVec))
				ret = false;
				return;
			end
		elseif (length(size(inputVec)) == 3)
			if (size(inputVec, 2) ~= 1)
				ret = false;
				return;
			end
		else
			ret = false;
			return;
		end
		ret = true;
	end

	p = inputParser;
	addRequired(p, 'lambda', @isvector);
	addRequired(p, 'xGrid', @(x) isnumeric(x) && ismatrix(x));
	addRequired(p, 'yGrid', @(x) isnumeric(x) && ismatrix(x));
	addRequired(p, 'elementLocs', @(x) ismatrix(x) && (size(x, 2) == 3));
	addRequired(p, 'inputVector', @validateInputVector)
	addParameter(p, 'n_medium', 1, @(x) isnumeric(x) && isscalar(x));
	addParameter(p, 'zCut', 0, @(x) isnumeric(x) && isscalar(x));
	addParameter(p, 'UseObliquity', false, @(x) islogical(x) && isscalar(x));
	addParameter(p, 'apertureNormal', false, @(x) isnumeric(x) && (length(x) == 3));
	addParameter(p, 'UseGPU', true, @(x) islogical(x) && isscalar(x));
	parse(p, lambda, xGrid, yGrid, elementLocs, inputVector, varargin{:});

	n_medium = p.Results.n_medium;
	zCut = p.Results.zCut;
	useObliquity = p.Results.UseObliquity;
	apertureNormal = p.Results.apertureNormal;
	useGPU = p.Results.UseGPU;

	if useObliquity && islogical(apertureNormal)
		error("Must specify 'apertureNormal' if 'UseObliquity' is set to true.");
	end



	g = zeros([size(xGrid) length(lambda)]);

	k_lambda = n_medium .* (2*pi ./ lambda);
	k_lambda = gpuArrayHelper(reshape(k_lambda, [1 1 1 length(k_lambda)]), useGPU);
	
	dx = gpuArrayHelper(xGrid - reshape(elementLocs(:,1), [1 1 size(elementLocs, 1)]), useGPU);
	dy = gpuArrayHelper(yGrid - reshape(elementLocs(:,2), [1 1 size(elementLocs, 1)]), useGPU);
	dz = gpuArrayHelper(zCut - reshape(elementLocs(:,3), [1 1 size(elementLocs, 1)]), useGPU);
	R = sqrt((dx.^2) + (dy.^2) + (dz.^2));

	if (useObliquity)
		apertureNorm = apertureNormal / sqrt(sum(apertureNormal .^ 2));
		obliquity = ((dx * apertureNorm(1)) + (dy * apertureNorm(2)) + (dz * apertureNorm(3))) ./ R;
	else
		obliquity = 1;
	end

	clear("dx", "dy", "dz");
	wait(gpuDevice);

	for m = 1:length(k_lambda)
		focusingVecTemp = gpuArrayHelper(reshape(inputVector(:,1,m), [1 1 size(inputVector, 1) 1]), useGPU);
		g(:,:,m) = sum(focusingVecTemp .* (exp(-j*k_lambda(1,1,1,m).*R) ./ R) .* obliquity, 3);
	end
end