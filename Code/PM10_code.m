%%DATI PM10 MILANO VERZIERE
data = importdata('PM10.txt');
Y=data.data
t1 = datetime(2018,1,1,9,0,0);
t2 = datetime(2021,12,30,9,0,0);
X = t1:t2;

n = length(Y)
i = n

while(i > 0)
    if Y(i) < 0
        Y(i) = mean(Y)
    end
    i=i-1
end

i = n
while(i > 0)
    if Y(i) > 3*mean(Y)
        Y(i) = mean(Y)
    end
    i=i-1
end

i = n
while(i > 0)
    if Y(i) < 0
        Y(i) = 0
    end
    i=i-1
end

figure
plot(Y, 'Color' , '#7E2F8E', 'LineWidth', 1.0)
title('PM10')

%Ym = Y - mean(Y);% observation mean removal
N = length(Y)
N2 = N/2;                   % max length for ecf evaluation
ecf = xcorr(pm10,N2, 'unbiased');			% ecf 

figure
plot(ecf,'.r')
title('ecf unbiased')

ecf = ecf(N2+1:2*N2+1);

%covariance model parameter estimation
% setting up the starting value of A
Ast01 = 2*ecf(2) - ecf(3); % line passing on the first two covariances
Ast02 = 3*ecf(2) - 3*ecf(3) + ecf(4); % parabola passing on the first 3 values

if ecf(3) < (ecf(2)+ecf(4))/2 % concavity up
Ast0 = Ast01;
else
Ast0 = Ast02;
end

sigma2st0_v = ecf(1) - Ast0;

ecf_05 = Ast0/2;
i_lcorr = find((ecf-ecf_05)<=0, 1 ) ; % index of the first f(t) >= A/2;
t_lcorr = t(i_lcorr-1); % correletion lenght

ast01 = log(2)/t_lcorr;
ast02 = log(Ast0/ecf(i_lcorr-1))/t_lcorr;
ecm01 = Ast0 * exp(-ast01*t); % covariance function model
ecm02 = Ast0 * exp(-ast02*t); % covariance function model

figure
plot (0:199,ecf(1:200),'.r'); hold on
plot (0:199,ecm01(1:200),'Color' , '#D95319', 'LineWidth', 1.5); hold on
plot (0:199,ecm02(1:200),'Color' , '#7E2F8E', 'LineWidth', 1.5); hold on
plot (1, 0.36, 'ob'); grid
legend('Empirical Values', 'Straight Line Interpolation', 'Polynomial of Degree 2', 'correlation length')
title('Comparison of Interpolation Methods')

fprintf('\n');
fprintf(['Approximated values of the parameters (I selected the ''*'' one)\n ' ...
'where Ast01 is for line interpolation and Ast02 for parabola: \n']);
fprintf(' Ast01 = %f',Ast01);

if (Ast01 == Ast0)
fprintf(' (*) \n');
else
fprintf('\n');
end
fprintf(' Ast02 = %f',Ast02);

if (Ast02 == Ast0)
fprintf(' (*) \n');
else
fprintf('\n');
end

fprintf(' ast01 = %f (*) \n',ast01);
fprintf(' ast02 = %f \n',ast02);
fprintf(' sigmast0_v = %f \n',sqrt(sigma2st0_v));
fprintf('\n');
fprintf('Correlation length: \n');
fprintf(' t = %f \n',t_lcorr);

%% number of points to use
i_up = 200;
t_up = t(i_up);
tQ = t(2:i_up); % interpolation points
ecfQ = ecf(2:i_up); % values to interpolate
Ast = Ast0; % iteration zero
ast = ast01;
ecm = Ast * exp(-ast*t);% covariance model - it is the same of before


figure
plot (t(1:200),ecf(1:200),'.r');
hold on
plot (t(1:200),ecm(1:200),'Color', '#EDB120', 'LineWidth', 2.0);
hold on;
plot (t(1:200),ecmx(1:200),'Color', '#7E2F8E', 'LineWidth', 2.0);
hold on;
plot (t(1:200),ecmy(1:200),'Color', '#D95319', 'LineWidth', 2.0);
hold on;
plot (t(1:200),ecmy(1:200),'Color', '#0072BD', 'LineWidth', 2.0);
grid on;
title("LS interpolation")
legend('ecf','LS interpolated function', '1st interpolation', '2nd interpolation', '3rd interpolation')

%We only iterate once, because otherwise we would end up with a null matrix
%after the first iteration
A = [exp(-ast*tQ) -Ast*tQ.*exp(-ast*tQ)]; % design matrix
a = Ast*exp(-ast*tQ); % array of the known terms
xst = (A'*A)\(A'*(ecfQ-a)); % estimation parameter
Astx = Ast+xst(1); % iteration i
astx = ast+xst(2);
sigma2st_vx = ecf(1) - Astx;

ecmx = Astx * exp(-astx*t);

%SECOND ITERATION
A1 = [exp(-astx*tQ) -Astx*tQ.*exp(-astx*tQ)]; % design matrix
a1 = Astx*exp(-astx*tQ); % array of the known terms
xst = (A1'*A1)\(A1'*(ecfQ-a1)); % estimation parameter
Asty = Astx+xst(1); % iteration i
asty = astx+xst(2);
sigma2st_vy = ecf(1) - Asty;

ecmy = Asty * exp(-asty*t);

%%THIRD ITERATION
A2 = [exp(-asty*tQ) -Asty*tQ.*exp(-asty*tQ)]; % design matrix
a2 = Asty*exp(-asty*tQ); % array of the known terms
yst = (A2'*A2)\(A2'*(ecfQ-a2)); % estimation parameter
Astz = Asty+yst(1); % iteration i
astz = asty+yst(2);
sigma2st_vz = ecf(1) - Asty;

ecmz = Astz * exp(-astz*t);

fprintf('\n');
fprintf('Estimated values of the parameters: \n');
fprintf(' Ast = %f \n',Asty);
fprintf(' ast = %f \n',asty);
fprintf(' Estimation of the variance of the noise = %f \n',sqrt(sigma2st_vy));
fprintf('\n');

%% estimation with collocation

n = 200
Cyy = toeplitz(ecmx(1:n)); % covariance matrix of the signal
%Cyy(1:5,1:5)
Cy0y0 = Cyy + eye(n)*sigma2st_vx; % covariance matrix of the observations
%Cy0y0(1:5,1:5)
yst = Cyy * (Cy0y0\pm10(1:n)); % filtered signal
est = pm10(1:n) - yst;

figure
plot(pm10(1:n),'.-b'); hold on
plot(yst,'.-k');
title ('Signal vs. Filtered Signal (1st iteration)')
legend('signal', 'filtered signal');

fprintf('\n');
fprintf('estimation error: \n');
fprintf('   mean  = %f \n',mean(est));
fprintf('   std = %f \n',std(est));
fprintf('\n');