% FarleyPatrick_phys341_whirlpool_model.m
% Oct 25, 2014
%
% Plot the shape of a whirlpool in a body of water, whose paramters are 
% based on the description of Charybdis in Homer's The Odyssey:


clear all;

%% Calculate the velocity profile of the water before whirlpool begins:
% Whirlpools do not create circulation in a body of water. Rather, they 
% concentrate the existing circulation to a smaller and smaller radius. 
% Since angular momentum is conserved, the fluid simply flows faster as 
% it gets nearer to the axis of the whirlpool.

% Based on some information in the story (which I will explain in the
% paper), We will say the whirlpool occurs 75 meters away from shore,
% in a strait that is 200 meters across and 75 meters deep, with steep walls.

% For the current of the water, I will start with a value taken from a report on a
% similar strait of water: about .1 m/s. But, the water flow will not be
% uniform. Rather, it will be slower on the sides and faster in the middle 
% due to viscous effects. The current flow distribution is estimated as follows:

x = 0:1:200; % distance from the left shore (m)
waterVelocity = RiverVelocity(x); % veloctiy of water flow (m/s)

% the water velocity profile can then be visualized:
plot(x,waterVelocity); hold on;

% add in a marker for the whirlpool location:
plot(75,0:.01:.1, '*r'); hold off;

% format plot:
axis([0 200 0 .3]);
title('Default Current Flow in River')
xlabel('location across river (m)')
ylabel('veloctity of water (m/s)')
text(80,.1,'whirlpool axis')
text(140,.12,'velocity curve')

set(gca,'FontSize',16,'fontWeight','bold')

set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')

%% Find angular momentum in the water before whirlpool begins
% The next step is to calculate the angular momentum of the water with
% respect to the whirlpool's axis. In this case, what we actually want is
% angular momentum per unit of volume, since we are dealing with
% indeterminate volumes. L = m(r x v) becomes L/V = p(r x v) where p is the
% density of the water. 

v_cross_r_left = integral(@RiverV_cross_R,0,75);
v_cross_r_right = integral(@RiverV_cross_R,75,200);

%noting that:
density_water = 1020; % kg / m^3


% After integrating the velocities with respect to
% their horizontal dispacement from the axis, and multiplying by density, we get:
L_per_V_left = density_water * v_cross_r_left; % kg / ms
L_per_V_right = density_water * v_cross_r_right; % kg / ms

% the net angular momentum per volume about the axis will be the difference in
% angular momentum per volume on either side
L_per_V_net = abs(L_per_V_right - L_per_V_left); % kg / ms

% the direction of this vectory quantity is not important.

%%
L = 9.5*10^8
%% Find theta direction velocity:
r = 0:75;
f = @(r) L_per_V_net./(1000.*r.^2) .* 50*2*pi.*r
omg = integral(f,75,74);
%% calculate circular water velocity as a function of radius:
% Since we now know the angular momentum per volume, and we know that is is conserved
% in this system, we can get a value for the theta-velocity of the water
% in terms of different radii from the axis of the whirlpool, once the whirlpool
% starts. To clarify, this whirlpool is caused by a large potential drop
% over a small area at the bottom of the water (a creature sucking water
% in). This has the effect of pulling large volumes of water horizontally
% toward the axis of the whirlpool (r decreases). As r decreases, the
% velocity of this water must increase to conserve angular momentum.

% since L/V = p(r x v), and the added velocity will be tangent to the
% radius, we can approximate L/V = prv and thus v = L/pVr
r = 1:.5:150; % we will use a wide range for sample r values (m)
visc = .5;
v_theta = visc .* L_per_V_net ./ (density_water .* r); % veloctiy in the theta direction (m/s)

%% find the height of the water as it relates to theta-velocity:
% The next step is to model the height of the water as a function of its
% velocity. As whirlpools form, the higher-velocity water has a
% lower pressure because of the Bernoulli principle. As a result, this
% water yeilds to atmospheric pressure, and the whirlpool forms its iconic
% sunken hole shape. By solving the Bernoulli equation:
% (P + .5pv^2 + pgh = constant), we get h = -(v^2)/(2g) + const, and the
% constant turns out to be the starting height of the water. so:
g = 9.8; % gravitational constant (m/s^2)
h = -101325./(density_water*g)  -(Vtheta(L,rho,h,r)^2) ./ (2*g) + 50; % height of water surface from floor (m)

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
axis([-75 125 0 50]);
title('Cross-Section of Charybdis whirlpool');
xlabel('distance from whirlpool axis (m)');
ylabel('height of water (m)');

%% Summary:
% we now have an accurate function to model the shape of this particular
% whirlpool, and we have the theta-directional velocity as a function of location.
% Future steps would be to delve into the z- and r- direction velocities of 
% the water, and one possible appliction would be finding the whirlpool's 
% "event horizon" for a ship with a defined top-speed. This would reveal 
% just how great a hazard charybdis' whirlpool would have been.