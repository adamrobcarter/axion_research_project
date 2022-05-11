f_sampling = 2e9; % Hz
T = 10e-3; % K

h = 6.62607015e-34; % m2 kg / s
k_B = 1.38065e-23; % m2 kg s-2 K-1

dPdf = @(f) 0.5*h*f + h*f/(exp(h*f/k_B/T)-1);
noise_power = integral(dPdf, 0, f_sampling/2) * 10;

axion_power = 1.52e-21; % W, from eq 2.28 of Daw thesis

mult1 = [2.2336    1.1489    1.5633    5.0113    3.4929]; % tau 5e-4
mult10 = [0.9777    1.8556    1.0693    0.5650    0.2224]; % tau 20e-4

x = [axion_power/noise_power, axion_power/noise_power/10];
y = [mean(mult1), mean(mult10)];
err = [std(mult1), std(mult10)];

%set(0,'defaulttextinterpreter','latex')

scatter(x, y)
hold on
errorbar(x, y, err/2, 'LineStyle','none'); % err/2 cause matlab is silly, check the docs
xlabel('$P_{sig} / P_{noise}$')
ylabel('$SNR_{dynamic} / SNR_{static}$')
adstyle(8, 8)
%saveas(gcf,'figures/main.eps')