function [b, a] = my_iirpeak(f_0, Q, f_sampling)
    % lifted straight from scipy source code: https://github.com/scipy/scipy/blob/v1.8.0/scipy/signal/_filter_design.py#L4946-L5023
    w0 = f_0 / (f_sampling / 2); % normalise to nyquist frequency

    assert(w0 > 0)
    assert(w0 < 1)

    w0 = w0*pi;

    % Compute beta: formula 11.3.19 (p.579) from reference [1]
    beta = tan(w0/(2*Q));

    % Compute gain: formula 11.3.6 (p.575) from reference [1]
    gain = 1 / (1 + beta);

    % Compute numerator b and denominator a
    % formulas 11.3.7 (p.575) and 11.3.21 (p.579)
    % from reference [1]
    b = (1 - gain) .* [1, 0, -1];
    a = [1.0, -2*gain*cos(w0), 2*gain - 1];
end