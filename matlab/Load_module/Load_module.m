classdef Load_module
  % branch from 11DOF.cpp
   properties
       Active_Profile
       Car_obj
   end
   methods
       function obj = Load_module(Profile_type, car1)
           pos_vec_rand = rand(3,1); % random position generator
           if nargin == 2
               obj.Active_Profile = Profile_type;
               obj.Car_obj = car1;
           elseif nargin == 1
               obj.Active_Profile = Profile_type;
               obj.Car_obj = Car();
               obj.Car_obj.Position = obj.Active_Profile.Position + ...
                                      obj.Active_Profile.Radius * pos_vec_rand / norm(pos_vec_rand);
           elseif nargin == 0
               obj.Active_Profile = Circular();
               obj.Car_obj = Car();
               obj.Car_obj.Position = obj.Active_Profile.Position + ...
                                      obj.Active_Profile.Radius * pos_vec_rand / norm(pos_vec_rand);
           end
           if abs(norm(obj.Active_Profile.Position-obj.Car_obj.Position)-obj.Active_Profile.Radius) > 1e-12
             % !!! Be careful to not violate this condition - assert in CPP
               disp('ERROR!!!! Distance between car and center of curve is different then radius');               
               obj.Active_Profile.Position = obj.Active_Profile.Position - ...
                                             obj.Active_Profile.Radius * pos_vec_rand / norm(pos_vec_rand);
               disp(abs(norm(obj.Active_Profile.Position-obj.Car_obj.Position)));
           end    
       end
       
       function obj = set_Profile(obj, Profile_type)
           obj.Active_Profile = Profile_type;
       end
       
       function normal_force = update_normal_force(obj, time, force_field)
           % call Profile_type update force
           normal_force = obj.Active_Profile.update_normal_force(obj.Car_obj, time, force_field);
       end
       
   end
end