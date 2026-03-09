close all; clear all; clc;

period = 0:12;
frequency = [0.038003841,0.22533590278469698,0.10825336526608549, ...
             0.086756237,0.09673705, 0.086756237, ...
             0.067562372,0.042226478282407345, 0.034165057, ...
             0.020729367816909986, 0.009596906, 0.018809958, ...
                0.0065259];

% Ensure sum of frequences is 1 (it's not - likely error in paper graph)
frequency_sum = sum(frequency);


% Plot distribution
edges = 0:13;
figure;
histogram('BinEdges', edges, 'BinCounts', frequency,'FaceColor', 'r','LineWidth', 2)
xlabel('period (min)')
ylabel('frequency')

% sale such that it sums to 1
multiplier = 1/frequency_sum;
frequency_scaled = frequency*multiplier;
figure;
histogram('BinEdges', edges, 'BinCounts', frequency_scaled,'FaceColor', 'r','LineWidth', 2)
xlabel('period (min)')
ylabel('frequency')

%% Fit distribution to data
N = 91;
counts = round(frequency_scaled*N);
x = 0.5:1:12.5;
data = repelem(x,counts);
distribution = 'gamma';
x_pdf = linspace(min(edges), max(edges), 100);


switch distribution
    case "normal"
        pd = fitdist(data', 'normal');
        PDF = normpdf(x_pdf,pd.mu,pd.sigma);
        rv = normrnd(pd.mu, pd.sigma, [1 100]);
    case "lognormal"
        pd = fitdist(data','lognormal');
        PDF = lognpdf(x_pdf,pd.mu,pd.sigma);
        rv = lognrnd(pd.mu, pd.sigma, [1 100]);
        mu = exp(pd.mu+0.5*pd.sigma^2);
        sigma = sqrt(mu^2*(exp(pd.sigma^2)-1));
    case "burr"
        pd = fitdist(data','burr');
        PDF = pdf('burr',x_pdf,pd.alpha,pd.c,pd.data);
        rv = random('burr',pd.alpha,pd.c,pd.data,[1 100]);
        mu = pd.data * beta(pd.data - 1/pd.c,1+1/pd.c);
        %mu = gamma(pd.alpha-1/pd.c)*gamma(1+1/pd.c)*pd.data^(1/pd.c)/gamma(pd.alpha)
        %sigma = (gamma(pd.alpha - 2/pd.c)*gamma(1+2/pd.c)-(gamma(pd.alpha-1/pd.c)*gamma(1+1/pd.c))^2/gamma(pd.alpha))*pd.data^(2/pd.c)/gamma(pd.alpha)

    case "gamma"
        pd = fitdist(data','gamma');
        PDF = pdf('gamma',x_pdf,pd.a,pd.b);
        rv = random('gamma',pd.a,pd.b,[1 N]);
    case "gev"
        pd = fitdist(data', 'gev');
        PDF = pdf('gev', edges,pd.data,pd.sigma,pd.mu);
        rv = random('gev', pd.data,pd.sigma,pd.mu,[1 N]);
    case "tLocationScale"
        pd = fitdist(data', 'tLocationScale');
        PDF = pdf('tLocationScale',edges,pd.mu,pd.sigma, pd.nu);
        rv = random('tLocationScale',pd.mu,pd.sigma,pd.nu, [1 N]);
    case "stable"
        pd = fitdist(data','stable');
        PDF = pdf('Stable',edges,pd.alpha,pd.beta, pd.gam,pd.delta);
        rv = random('stable',pd.alpha,pd.beta,pd.gam,pd.delta, [1 N]);
    case "loglogistic"
        pd = fitdist(data','loglogistic');
        PDF = pdf('Loglogistic',x_pdf,pd.mu,pd.sigma);
        rv = random('Loglogistic',pd.mu,pd.sigma, [1 N]);
        alpha = exp(pd.mu);
        beta = 1/pd.sigma;
        mu = alpha*pi/(sin(pi/beta)*beta);
        sigma = alpha*sqrt(2*pi/sin(2*pi/beta)/beta-(pi/beta)^2/sin(pi/beta)^2);
    case "kernel"
        pd = fitdist(data','kernel');
        PDF = pdf(pd,x_pdf);
        rv = random(pd, [1 N]);

    otherwise
        error("wrong input")
end

figure(3);
hold on;
histogram('BinEdges', edges, 'BinCounts', counts, 'Normalization', 'pdf','facecolor', 'g','facealpha',.4,'edgealpha',0,'edgecolor','r');
plot(x_pdf,  PDF, 'r-', 'LineWidth', 2);

max_rv = max(rv);
edges_rv = 0:1:(max_rv+1);
histogram(rv,edges_rv,'Normalization', 'pdf','facecolor', 'r','facealpha',.4,'edgealpha',0,'edgecolor','r');