function f_iwave = do_iwave_and_resonator(t, noise, signal, f_0_start, band_width, f_sampling, limiting_band, tau)

f_guess = f_0_start;

STATE_PS = iwave_setup(f_sampling, tau, f_guess, 0);
state_res = resonator_setup();
f_iwave = zeros(size(t));

f_0_current = f_0_start;

for i = 1 : length(t)
    resonator_input = noise(i) + signal(i);
    [state_res, resonator_output] = resonator_step2(state_res, resonator_input, f_0_current, band_width, f_sampling);
    iwave_input = resonator_output;
    [STATE_PS, pdout, pqout, paout, peout, pfout, pfstateout] = iwave_step(STATE_PS, iwave_input);
    f_0_current = pfout;
    assert(~isnan(f_0_current), 'iwave output was nan')
    if f_0_current > limiting_band(2)
        f_0_current = limiting_band(2);
    elseif f_0_current < limiting_band(1)
        f_0_current = limiting_band(1);
    end
    f_iwave(i) = f_0_current;
end

end

