classdef solver < handle
    % This class is the mother of all solver implementations
    properties (GetAccess = 'public', SetAccess='public')  
    tEnd;     %time at simulaton end 
    h;        %time step
    M;        %mass matrix
    D;        %damping matrix
    K;        %stiffness matrix
    u_init;   %vector of initial displacements 
    ud_init;  %vector of initial velocities 
    u_n_p_1;  %u_n+1
    u_n;      %u_n
    u_n_m_1;  %u_n-1 
    step;     %internal counter
    end
    methods (Access=public)
       %Constructor
       function obj=solver(h_,tEnd_,M_,D_,K_,u_init_,ud_init_)
       disp('Call constructor of abstract super-class');
       %Initialize with user input
       obj.h=h_;  
       obj.tEnd=tEnd_; 
       obj.M=M_; 
       obj.D=D_; 
       obj.K=K_; 
       obj.u_init=u_init_;
       obj.ud_init=ud_init_;
       
       obj.step=0;
       end
       function incStepCounter(obj)
       %Increment time step counter
       obj.step = obj.step+1;
       %Update vectors
       obj.u_n_m_1=obj.u_n;
       obj.u_n=obj.u_n_p_1;
       end 
    %methods public
    end
    
    methods(Abstract = true)
       [output,doutput,ddoutput]=doSolve(obj,input);
       [interfaceJacobian]=getInterfaceJacobian(obj,i,j);
       [interfaceJacobian]=getInterfaceJacobianVelocities(obj,i,j);
       [interfaceJacobian]=getInterfaceJacobianAccelerations(obj,i,j);
    %methods abstract
    end
    
    
%end classdef
end