function x_delayed = delay(x, delay_samples)
% delay a signal by the specified number of samples.
% pads the start with zeroes

delay_samples = floor(delay_samples); % needs to be integer

assert(delay_samples < length(x), "cannot delay by longer than input length")
assert(delay_samples > 0, "cannot delay by a negative number")

x_delayed = zeros(size(x));
x_delayed(delay_samples+1:end) = x(1:end-delay_samples);

end