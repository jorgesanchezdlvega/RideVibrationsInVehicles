clear
clc

%% Parametric Values

m1 = 600;
m2 = 45;
m3 = m2;
J = 1020;
L = 2.5;
p = 0.52;
q = 0.1;
k1 = 500e3;
k2 = 20e3;
k3 = k2;
c1 = 0;
c2 = 1e3;
c3 = c2;

m = 300;
mt = 45;
kt = 500e3;
k = 20e3;
ct = 0;
c = 1e3;

lambda = 16;
zeta = 4e-2;

%% Sub-Task 1

% For the four degree-of-freedom model

% 1)
data = readmatrix("bil0404.dat");
V = data(:,1); % in km/h
Vms = V*(5/18); % in m/s
z1 = data(:,2);
theta = data(:,3);
zf = data(:,4);

[maxZ1, maxZ1Index] = max(z1);
vMaxZ1 = V(maxZ1Index);
[maxTheta, maxThetaIndex] = max(rad2deg(theta));
vMaxTheta = V(maxThetaIndex);
[maxZf, maxZfIndex] = max(zf);
vMaxZf = V(maxZfIndex);

% 2)
figure1 = figure;
plot(V,z1,vMaxZ1,maxZ1,"ok")
xlabel('Speed [m/s]');
ylabel('Displacement [m]');
title("Car chassis displacement (four degree-of-freedom model)")
yline(maxZ1, "-.k", ['Max displacement: ', num2str(round(maxZ1,3)), ...
    ' m'], "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","bottom");
xline(vMaxZ1, "-.k", ['V for max displacement: ', ...
    num2str(round(vMaxZ1,2)), ' km/h'], ...
    "LabelHorizontalAlignment", "center", ...
    "LabelVerticalAlignment","middle");
figure2 = figure;
plot(V,rad2deg(theta),vMaxTheta,maxTheta,"ok", "LineWidth",1.2)
xlabel('Speed [m/s]');
ylabel('Rotation [deg]');
title("Car chassis rotation (four degree-of-freedom model)")
yline(maxTheta, "-.k", ['Max displacement: ', num2str(round(maxTheta,3)), ...
    ' m'], "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","bottom");
xline(vMaxTheta, "-.k", ['V for max displacement: ', ...
    num2str(round(vMaxTheta,2)), ' km/h'], ...
    "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","middle");
figure3 = figure;
plot(V,zf,vMaxZf,maxZf,"ok")
xlabel('Speed [m/s]');
ylabel('Displacement [m]');
title("Driver displacement (four degree-of-freedom model)")
yline(maxZf, "-.k", ['Max displacement: ', num2str(round(maxZf,3)), ...
    ' m'], "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","bottom");
xline(vMaxZf, "-.k", ['V for max displacement: ', ...
    num2str(round(vMaxZf,2)), ' km/h'], ...
    "LabelHorizontalAlignment", "center", ...
    "LabelVerticalAlignment","middle");

% For quarter car model

% 3) an 4)
f = Vms/lambda;
om = 2*pi*f;

z = ((k+om*c*1i).*(kt+om*ct*1i)*zeta)./((k+om*c*1i-om.*om*m).* ...
    ((k+kt)+om*1i*(c+ct)-om.*om*mt)-(k+om*c*1i).^2);

[maxZ, maxZIndex] = max(abs(z));
vMaxZ = V(maxZIndex);

figure4 = figure;
plot(V,abs(z),vMaxZ,maxZ,"ok")
xlabel('Speed [m/s]');
ylabel('Displacement [m]');
title("Driver displacement (quarter car model)")
yline(maxZ, "-.k", ['Max displacement: ', num2str(round(maxZ,3)), ...
    ' m'], "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","bottom");
xline(vMaxZ, "-.k", ['V for max displacement: ', ...
    num2str(round(vMaxZ,2)), ' km/h'], ...
    "LabelHorizontalAlignment", "center", ...
    "LabelVerticalAlignment","middle");

%% Sub-Task 2

% 1)
K = [k2+k3 -k2 -k3 (p*k2-(1-p)*k3)*L;
    -k2 k1+k2 0 -p*k2*L; -k3 0 k1+k3 (1-p)*k3*L;
    (p*k2-(1-p)*k3)*L -p*k2*L (1-p)*k3*L p*p*L*L*k2+L*L*k3*(1-p)*(1-p)];
M = [m1 0 0 0; 0 m2 0 0; 0 0 m3 0; 0 0 0 J];

% 2)
[Veig, D] = eig(K, M);
frecs = sqrt(diag(D))/(2*pi);
vcrits = frecs*lambda;

% 3)
vCrit1 = vcrits(1)/(5/18); % v to km/h
vCrit2 = vcrits(2)/(5/18); % v to km/h

% Only the two first values are relevant as the last two are much higher
% than 120 km/h .

labelV1 = "V_c_r_i_t_1: " + round(vCrit1,2) + " km/h";
labelV2 = "V_c_r_i_t_2: " + round(vCrit2,2) + " km/h";

figure(figure1)
xline(vCrit1, "--r", labelV1, ...
    "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","middle");
xline(vCrit2, "--r", labelV2, ...
    "LabelHorizontalAlignment", "right", ...
    "LabelVerticalAlignment","middle");
figure(figure2)
xline(vCrit1, "--r", labelV1, ...
    "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","middle");
xline(vCrit2, "--r", labelV2, ...
    "LabelHorizontalAlignment", "right", ...
    "LabelVerticalAlignment","middle");
figure(figure3)
xline(vCrit1, "--r", labelV1, ...
    "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","middle");
xline(vCrit2, "--r", labelV2, ...
    "LabelHorizontalAlignment", "right", ...
    "LabelVerticalAlignment","middle");

% 4)
eigenvector1 = Veig(:,1)
eigenvector2 = Veig(:,2)
eigenvector3 = Veig(:,3)
eigenvector4 = Veig(:,4)

%% Sub-Task 3

% 1)
c2 = linspace(1e3,10e3,length(V));
c3 = c2;
v = 40; % in km/h
fVSpecific = (v*(5/18))/lambda; % V in m/s
omVSpecific = 2*pi*fVSpecific;
fv = [0; 1; exp(1i*2*pi*L/lambda); 0]*(zeta)*(k1+1i*omVSpecific*c1);
maxDisp = 0.04; % in m

% 2)
z = zeros(4,length(c2));
for i = 1:length(c2)
    C = [c2(i)+c3(i) -c2(i) -c3(i) (p*c2(i)-(1-p)*c3(i))*L;
        -c2(i) c1+c2(i) 0 -p*c2(i)*L;
        -c3(i) 0 c1+c3(i) (1-p)*c3(i)*L;
        (p*c2(i)-(1-p)*c3(i))*L -p*c2(i)*L (1-p)*c3(i)*L ...
        p*p*L*L*c2(i)+L*L*c3(i)*(1-p)*(1-p)];
    z(:,i) = (K+1i*omVSpecific*C -omVSpecific*omVSpecific*M)\fv;
end

% 3)
zf = z(1,:) - q*L*z(4,:);
abszf = abs(zf);

% 4)
figure5 = figure;
plot(c2*1e-3,abszf)
xlabel('Damping coeficient value [kNs/m]');
ylabel('Displacement [m]');
title("Driver displacement")
yline(maxDisp, "--r", [num2str(maxDisp*1e2), ' cm'], ...
    "LabelHorizontalAlignment", "left", ...
    "LabelVerticalAlignment","bottom");

% 5)
dampCoefC2 = min(c2(abszf <= maxDisp))
figure(figure5)
hold on
plot(dampCoefC2*1e-3,maxDisp,"o")
xline(dampCoefC2*1e-3, "-.k", ['Damping Coeficient: ' ...
    , num2str(round(dampCoefC2*1e-3,2)), ' kNs/m'], ...
    "LabelHorizontalAlignment", "right", ...
    "LabelVerticalAlignment","top");
hold off