classdef Car
   properties
       Position {mustBeNumeric}
       Velocity {mustBeNumeric}
       Mass {mustBeNumeric}
       Length {mustBeNumeric} % needed when we don't roughly approximate N
       Width {mustBeNumeric} % needed when we don't roughly approximate N
       % ALL the physical properties
       % spring info (elongation for each spring)
       
       % lengths of each spring
       L_body_fl {mustBeNumeric}
       L_tire_fl {mustBeNumeric}
       L_body_fr {mustBeNumeric}
       L_tire_fr {mustBeNumeric}
       L_body_rl {mustBeNumeric}
       L_tire_rl {mustBeNumeric}
       L_body_rr {mustBeNumeric}
       L_tire_rr {mustBeNumeric}
   end
   methods
        % constructor with / without parameters (too long to put in a picture)
        function obj = Car(pos, vel, mass, len, wid, ...
                    l_body_fl, l_tire_fl, l_body_fr, l_tire_fr, ...
                    l_body_rl, l_tire_rl, l_body_rr, l_tire_rr)
            vel_default = rand(3,1); 
            vel_default = 60/3.6*vel_default/norm(vel_default); % default 16.666 m/s (60 km/h)
            if nargin == 0
                obj.Position = zeros(3,1);
                obj.Velocity = vel_default;
                obj.Mass = 10;
                obj.Length = 1;
                obj.Width = 0.5;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 1
                obj.Position = zeros(3,1);
                obj.Velocity = vel;
                obj.Mass = 10;
                obj.Length = 1;
                obj.Width = 0.5;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 2
                obj.Position = pos;
                obj.Velocity = vel;
                obj.Mass = 10;
                obj.Length = 1;
                obj.Width = 0.5;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 3
                obj.Position = pos;
                obj.Velocity = vel;
                obj.Mass = mass;
                obj.Length = 1;
                obj.Width = 0.5;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 4
                obj.Position = pos;
                obj.Velocity = vel;
                obj.Mass = mass;
                obj.Length = len;
                obj.Width = 0.5;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 5
                obj.Position = pos;
                obj.Velocity = vel;
                obj.Mass = mass;
                obj.Length = len;
                obj.Width = wid;
                obj.L_body_fl = 1;
                obj.L_tire_fl = 1;
                obj.L_body_fr = 1;
                obj.L_tire_fr = 1;
                obj.L_body_rl = 1;
                obj.L_tire_rl = 1;
                obj.L_body_rr = 1;
                obj.L_tire_rr = 1;
            elseif nargin == 13
                obj.Position = pos;
                obj.Velocity = vel;
                obj.Mass = mass;
                obj.Length = len;
                obj.Width = wid;
                obj.L_body_fl = l_body_fl;
                obj.L_tire_fl = l_tire_fl;
                obj.L_body_fr = l_body_fr;
                obj.L_tire_fr = l_tire_fr;
                obj.L_body_rl = l_body_rl;
                obj.L_tire_rl = l_tire_rl;
                obj.L_body_rr = l_body_rr;
                obj.L_tire_rr = l_tire_rr;
            end
        end
        function pos = get_position(obj)
            pos = obj.Position;
        end
        function vel = get_velocity(obj)
            vel = obj.Velocity;
        end
        function mass = get_mass(obj)
            mass = obj.Mass;
        end
        function dist_vector = get_dist_vector(obj, Point_P)
           % vector from the center of mass of the car to a specific
           % point from the global reference frame
           % Inputs: obj, Point_P - 3D vector
           dist_vector = Point_P - obj.Position;
        end
   end
end