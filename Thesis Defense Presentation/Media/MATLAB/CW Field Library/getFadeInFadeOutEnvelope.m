% For fading fields in and out in time-domain plots.
function env = getFadeInFadeOutEnvelope(lengths)
	% Lengths specify the lengths of these stages (in order):
	%	Off --> Fade in --> Normal --> Fade out --> Off

	% Off
	env = zeros(1, lengths(1));

	% Fade in
	temp = linspace(0, 1, lengths(2) + 2);
	env = [env temp(2:end-1)];

	% Normal
	env = [env ones(1, lengths(3))];

	% Fade out
	temp = linspace(1, 0, lengths(4) + 2);
	env = [env temp(2:end-1)];

	% Off
	env = [env zeros(1, lengths(5))];
end