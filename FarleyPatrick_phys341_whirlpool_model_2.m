% FarleyPatrick_phys341_whirlpool_model_2.m
% Oct 25, 2014
%
% Plot the shape of a whirlpool in a body of water, whose paramters are 
% based on the description of Charybdis in Homer's The Odyssey:


clear all;

%%
L =  9.5*10^8;
Q = 10; %m^3 /s
rho = 1020;
d = 50;
r = 0:.1:25;
R = 75
P = 101325;
g = 9.8;
V = pi*R^2 * d;
Gamma = 2*pi*L/V;
%% calculate circular water velocity as a function of radius:
% see function Vtheta

%% find the height of the water as it relates to theta-velocity:
% The next step is to model the height of the water as a function of its
% velocity. As whirlpools form, the higher-velocity water has a
% lower pressure because of the Bernoulli principle. As a result, this
% water yeilds to atmospheric pressure, and the whirlpool forms its iconic
% sunken hole shape. By solving the Bernoulli equation:
% (P + .5pv^2 + pgh = constant), we get h = -(v^2)/(2g) + const, and the
% constant turns out to be the starting height of the water. so:

v_t = Vtheta(L,V,rho,d,r);
v_r = Q./(2*d*pi.*r);
v_tot = sqrt(v_r.^2 + v_t.^2);
h = -P/(rho*g)  -(v_tot.^2) ./ (2*g) + (d+10.136); % height of water surface from floor (m)

% NOTICE: there is a simplification made here. The velocity used in the
% Bernoulli equation was only the theta-direciton component. There is also
% z and r direction velocity, but these are overlooked for the time being,
% since the theta-direction velocity is dominant.

% now we have related the height of the water to the radius, through the
% theta-velocity variable. 

%% plot water height vs radius:
% We can plot the height of the whirlpool with respect to its radius, giving
% an accurate cross-sectional view of the whirlpool:

figure;
plot(r,h); hold on;

% add the mirrored function to get a better picture:
plot(-r,h); hold off;

% format plot:
axis([-10 10 48 51]);
title('Cross-Section of Charybdis whirlpool');
xlabel('distance from whirlpool axis (m)');
ylabel('height of water (m)');


set(gca,'FontSize',16,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

%% Summary:
% we now have an accurate function to model the shape of this particular
% whirlpool, and we have the theta-directional velocity as a function of location.
% Future steps would be to delve into the z- and r- direction velocities of 
% the water, and one possible appliction would be finding the whirlpool's 
% "event horizon" for a ship with a defined top-speed. This would reveal 
% just how great a hazard charybdis' whirlpool would have been.