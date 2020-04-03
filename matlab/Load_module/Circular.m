classdef Circular < Profiles
   properties
       Position % center of geometric figure
       Radius % in case of circular motion
   end
   methods
       function obj = Circular(pos, radius)
           obj.Name = 'Circular';
           obj.Name = '0'; % velodrome default
           if nargin == 0
                obj.Position = zeros(3,1);
                obj.Radius = 100;
            elseif nargin == 1
                obj.Position = pos;
                obj.Radius = 100;
            elseif nargin == 2
                obj.Position = pos;
                obj.Radius = radius;
            end
       end
       
       function normal_force = update_normal_force(obj, car_obj, time, force_field)
           % input: force vector
           % return normal forces on each wheel (Ni=N/4)
           g = 9.81; % gravitational acceleration (on y direction)
           % normal_force = Profile_motion.update_normal_force(time,force_field);
           normal_force_global = -0.25 * ((car_obj.Mass * sum(car_obj.Velocity.^2) / obj.Radius) * ...
                                            car_obj.get_dist_vector(obj.Position) ... % centripet force
                                           -car_obj.Mass * [0;g;0] + force_field); 
           normal_force = [normal_force_global; ...
                           normal_force_global; ...
                           normal_force_global; ...
                           normal_force_global];
       end
   end
end