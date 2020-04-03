classdef Profile_description
    
%     Base_class -> Profiles
%     Derived_class -> Circular
%       Load_module uses Circular profile
%        - methods: 
%           - update_normal_force (
%           - get_position
%       Load_module receives an object of Circular class (all the methods
%       from load_module - position, radius, etc.)
%       
%       
   
%     construct with initial profile type given
%     member: active profile
%     get_position
%     update_force
%     switcher (profile_type) only called from outside
%   
   properties
       Name % name to describe the motion type (e.g. 'Velodrome' or '0')
       Position % center of geometric figure
       Radius % in case of circular motion
   end
   methods
        function obj = Profile_description(name, pos, radius)
            if nargin == 0
                obj.Name = '0'; % velodrome default
                obj.Position = zeros(3,1);
                obj.Radius = 100;
            elseif nargin == 1
                obj.Name = name;
                obj.Position = zeros(3,1);
                obj.Radius = 100;
            elseif nargin == 2
                obj.Name = name;
                obj.Position = pos;
                obj.Radius = 100;
            elseif nargin == 3
                obj.Name = name;
                obj.Position = pos;
                obj.Radius = radius;
            end
        end
        function pos = get_position(obj)
            pos = obj.Position;
        end
   end
end
