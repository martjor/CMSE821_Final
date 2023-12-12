classdef Crystal_growth
    %CRYSTAL_GROWTH This class simulates the dendritic crystal growth
    %   Detailed explanation goes here
    
    properties
        % Allen-Cahn equation
        M
        alpha
        epsilon
        gamma
        Tm

        % Heat Equation
        K

        % Discretized domain
        x
        y
        h
        sz

        % Function handles
        phi_initial_condition
        T_initial_condition

        % Stencils
        stencil_phi
        stencil_T
        stencil_grad_x 
        stencil_grad_y 

        % Indices
        p

        % Boundary conditions
        phi_bc
        T_bc
    end

    properties (Dependent)
        T_naught
        phi_naught
    end
    
    methods
        function obj = Crystal_growth(params,...
                                      phi_naught,...
                                      T_naught,...
                                      phi_bc,...
                                      T_bc)
            arguments
                params
                phi_naught
                T_naught
                phi_bc = @(phi) Crystal_growth.boundary_conditions(phi)
                T_bc = @(T,x) Crystal_growth.boundary_conditions(T)
            end

            %CRYSTAL_GROWTH Constructs an instance of the problem given the
            %   set of parameters contained in params.

            %   Set Allen-Cahn parameters
            obj.M = params.M;
            obj.alpha = params.alpha;
            obj.epsilon = params.epsilon;
            obj.gamma = params.gamma;
            obj.Tm = params.Tm;

            % Set heat equation parameters
            obj.K = params.K;

            % Set discretization
            obj.h = params.h;
            obj.x = params.xlim(1):obj.h:params.xlim(2);
            obj.y = params.ylim(1):obj.h:params.ylim(2);
            [obj.x,obj.y] = meshgrid(obj.x,obj.y);

            % Determine points accesible by stencils
            m = length(obj.y);
            n = length(obj.x);
            obj.sz = [m,n];

            obj.p = reshape((2:m-1)' + (1:n-2) * m,1,1,[]);

            % Determine stencils
            obj.stencil_phi = params.stencil_phi * (1/obj.h^2);
            obj.stencil_T   = params.stencil_T   * (1/obj.h^2);
            obj.stencil_grad_x = [0 0  0; -1 0 1; 0 0 0] * (1/(2*obj.h));
            obj.stencil_grad_y = [0 -1 0;  0 0 0; 0 1 0] * (1/(2*obj.h));

            % Function handles
            obj.phi_initial_condition = phi_naught;
            obj.T_initial_condition = T_naught;
            obj.phi_bc = phi_bc;
            obj.T_bc = T_bc;
        end

        function val = get.phi_naught(obj)
            val = obj.phi_initial_condition(obj.x,obj.y);
        end

        function val = get.T_naught(obj)
            val = obj.T_initial_condition(obj.x,obj.y);
        end
        
        function [phi,T] = step(obj,phi,T,k)
            %METHOD1 Summary of this method goes here
            % Create copy of previous phi
            phi_prev = phi; 

            % Take one step in phi field
            phi = phi + k * -obj.M * (obj.df_dphi(phi) - ...
                                      obj.laplacian(phi,...
                                                    obj.stencil_phi,...
                                                    epsilon=obj.epsilon) +...
                                      obj.m(phi,T));

            % Take one step in Temperature field
            T = T + k * (obj.laplacian(T,obj.stencil_T) + obj.K * (phi-phi_prev)/k);

            % Boundary conditions
            phi = obj.phi_bc(phi);
            T = obj.T_bc(T,obj.x);
        end

        function res = df_dphi(obj,phi)
            %DF_DPHI Calculates the first term in the Allen-Cahn Equation
            res = 1/2 * phi .* (1-phi) .* (1-2*phi);
        end

        function res = m(obj,phi,T)
            %M Calculates the third term in the Allen-Cahn Equation
            res = obj.alpha/pi * phi .* (1-phi) .* ...
                  atan(obj.gamma * (T-obj.Tm));
        end

        function res = laplacian(obj,f,stencil_arr,args)
            arguments
                obj
                f
                stencil_arr
                args.epsilon = 1;
            end
            
            % Check for coefficient type
            if isnumeric(args.epsilon)
                % Constant coefficient case. The laplacian can be evaluated
                % without consideration of an additional scalar field
                res =  args.epsilon^2 * stencil(f,obj.p,stencil_arr);

            elseif isa(args.epsilon,'cell')
                % Nonconstant coefficient (scalar field). 

            end
        end

        function [fx,fy,mag,idx] = grad(obj,f)
            fx = stencil(f,obj.p,obj.stencil_grad_x);
            fy = stencil(f,obj.p,obj.stencil_grad_y);

            mag = hypot(fx,fy);
            idx = mag > 0.01 * obj.h;
        end

        function [nx,ny,theta] = inward_vector(obj,f)
            % Allocate memory
            nx = zeros(obj.sz);
            ny = zeros(obj.sz);
            theta = zeros(obj.sz);

            % Calculate gradient
            [fx,fy,mag,idx] = obj.grad(f);

            % Calculate inward vector
            nx(idx) = fx(idx) ./ mag(idx);
            ny(idx) = fy(idx) ./ mag(idx);
            
            % Calculate angle
            theta(idx) = atan2(ny(idx),nx(idx));
        end
    end

    methods (Static)
        function U = boundary_conditions(U)
            U(1,2:end-1) = U(2,2:end-1);          % North row
            U(end,2:end-1) = U(end-1,2:end-1);    % South row 
            U(1:end,1) = U(1:end,2);              % West row
            U(1:end,end) = U(1:end,end-1);        % East row
        end
    end
end

