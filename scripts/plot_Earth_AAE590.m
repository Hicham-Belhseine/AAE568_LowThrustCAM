function [h] = plot_Earth_AAE590(npts_eq, Rotation, Texturemap)
% Plots the Earth as an ellipsoidal surface for a visual reference for AAE 590
% STM homework.
% Author: Nathan Houtz (2023-10-19)
% 
% Inputs:
%   npts_eq: [](1x1)<integer/double> position vector
%   Rotation: [deg](1x1)<double> Greenwich Sidereal Time (angle to rotate about
%             the z-axis)
%             OR
%             [](3x3)<double> ITRF --> J2000/GCRS rotation matrix such that:
%                   [r_GCRS] = [Rotation] * [r_ITRS]
%   texturemap: [](1x*)<char> file name of a texturemap on the Matlab Path
% 
% Outputs:
%   h: [](1x1)<graphics handle> handle for the surface object, which can be
%      referenced by legend()
% 
% ------------------------------- Example usage --------------------------------
% % Define some things to plot as an example:
% >> npts_eq = 4; %[] number of points around the equator
% >> GMST = 198.723095; %[deg] Greenwich Sidereal Time @ 2023-10-01 12:34:56
% >> C_ITRS2GCRS = [-0.94876977,    0.31595992,     0.002293967; ...
%                   -0.31596067,   -0.94877229,     3.5655734e-05; ...
%                    0.0021877177, -0.00069097413   0.999997368];
% >> % Texturemap = 'peppers.png'; %[] texturemap or image name
% >> Texturemap = 'kobi.png'; %[] Image that comes with Matlab
% >> % Texturemap = 'colorCheckerTestImage.jpg'; %[] Image that comes with Matlab
% 
% % Plot the Earth with 50 points around the equator in the ITRF
% >> plot_Earth_AAE590()
% 
% % Plot the Earth with <npts_eq> points around the equator in ITRF
% >> plot_Earth_AAE590(npts_eq)
% 
% % Plot the Earth rotated by GAST for the date 2023-10-01 12:34:56, roughly
% % giving the orientation in J2000:
% >> plot_Earth_AAE590([], GAST)
% 
% % Plot the Earth with the ITRS --> GCRS rotation matrix for the date
% % 2023-10-01 12:34:56
% >> plot_Earth_AAE590([], C_ITRS2GCRS)
%
% % Plot an ellipsoid the size of the Earth with a different image in ITRF:
% >> plot_Earth_AAE590([], [], Texturemap)

% Check the inputs
if nargin<1 || isempty(npts_eq)
    npts_eq = 50;
end
if nargin<2 || isempty(Rotation)
    R_ITRS2GCRS = eye(3,3);
elseif isscalar(Rotation)
    c = cosd(Rotation);
    s = sind(Rotation);
    R_ITRS2GCRS = [ c,      -s,     0; ...
                    s,      c,      0; ...
                    0,      0,      1   ];
elseif all(size(Rotation)==[3,3])
    R_ITRS2GCRS = Rotation;
else
    [M, N] = size(Rotation);
    error(['<Rotation> must be either a scalar angle in degrees, or a ', ...
        '(3x3) rotation matrix that rotates from ITRS to GCRS. Received:', ...
        '\n\tsize(Rotation) = [%i, %i]'])
end
if nargin<3 || isempty(Texturemap)
    Texturemap = imread('2k_Earth_July_topography_bathymetry.png');
else
    Texturemap = imread(Texturemap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin function:

% Generate points representing Earth and rotate them correctly:
RE_eq = 6378.137; %[km] Earth equatorial radius (WGS 84)
RE_pol = 6356.752314245; %[km] Earth polar radius (WGS 84)

% Generate an ellipsoid and rotate it appropriately:
[Earthx, Earthy, Earthz] = ellipsoid(0, 0, 0, RE_eq, RE_eq, RE_pol, npts_eq);
Rotated_points = R_ITRS2GCRS * [Earthx(:).'; Earthy(:).'; Earthz(:).'];
Earthx = reshape(Rotated_points(1,:), npts_eq+1, npts_eq+1);
Earthy = reshape(Rotated_points(2,:), npts_eq+1, npts_eq+1);
Earthz = reshape(Rotated_points(3,:), npts_eq+1, npts_eq+1);

% Generate a surface plot of the Earth:
if nargout>=1
    h = surf(Earthx, Earthy, -Earthz, 'FaceColor', 'texturemap', ...
        'CData', Texturemap, 'EdgeColor', 'none');
else
    surf(Earthx, Earthy, -Earthz, 'FaceColor', 'texturemap', ...
        'CData', Texturemap, 'EdgeColor', 'none');
end
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
