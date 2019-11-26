classdef BEsolver < solver
    % Backward Euler time integrator
    properties (GetAccess = 'public', SetAccess='public') 
        A; %system matrix
        B; %system matrix
    end
    methods (Access=public)
        %Constructor
       function obj=BEsolver(h_,tEnd_,M_,D_,K_,u_init_,ud_init_)
            disp('Generating Object BEsolver');
            obj = obj@solver(h_,tEnd_,M_,D_,K_,u_init_,ud_init_); 
            obj.setInitialConditions();
       end       
       function [interfaceDisplacement, interfaceVelocity, interfaceAcceleration]=doSolve(obj,interfaceForce)
            obj.u_n_p_1=obj.A\(obj.B*obj.u_n-((1/(obj.h*obj.h))*obj.M)*obj.u_n_m_1+interfaceForce);
            interfaceDisplacement=obj.u_n_p_1;
            interfaceVelocity= (1/obj.h)*(obj.u_n_p_1 - obj.u_n) ;
            interfaceAcceleration= (1/(obj.h*obj.h))*(obj.u_n_p_1 - 2* obj.u_n + obj.u_n_m_1) ;
       end 
       function [interfaceJacobian]=getInterfaceJacobian(obj,i,j)
            invA=inv(obj.A);
            interfaceJacobian=invA(i,j);
       end 
       function [interfaceJacobian]=getInterfaceJacobianVelocities(obj,i,j)
            invA=inv(obj.A);
            interfaceJacobian=(1/obj.h)*invA(i,j);
       end
       function [interfaceJacobian]=getInterfaceJacobianAccelerations(obj,i,j)
            invA=inv(obj.A);
            interfaceJacobian=(1/(obj.h*obj.h))*invA(i,j);
       end        
       % For debugging
       function [interfaceForce]=getForce(obj,interfaceDisp)
            interfaceForce=obj.A*interfaceDisp-(obj.B*obj.u_n-((1/obj.h)*obj.M)*obj.u_n_m_1);
       end 
    %end methods public
    end
    %%   PROTECTED FUNCTIONS
    methods (Access=protected)
        function setInitialConditions(obj)
            obj.u_n=obj.u_init;
            obj.u_n_m_1=obj.u_init-obj.h*obj.ud_init;
            obj.A=((1/(obj.h*obj.h))*obj.M+(1/obj.h)*obj.D+obj.K);
            obj.B=((2/(obj.h*obj.h))*obj.M+(1/obj.h)*obj.D);
        end
    %end methods private
    end
%end classdef
end