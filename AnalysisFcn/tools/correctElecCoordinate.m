



a = [-3.33, -26.02, 25.66]; %first point coordinate
b = [23.92, -16.02, 30.67]; % another point coordinate

interval = 3.5; % contact intervals in mm
numPt = 15; % number of points 

linev = b-a;


dd = 0:interval:(numPt-1)*interval;
tinterval = dd./norm(linev);

interp_ptx=linev(1)*tinterval'+a(1);
interp_pty=linev(2)*tinterval'+a(2);
interp_ptz=linev(3)*tinterval'+a(3);
interp_pt = [interp_ptx interp_pty interp_ptz];
