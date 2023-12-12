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

        
    end

    properties (Dependent)
        T_naught
        phi_naught
    end

    properties (Hidden)
        % Function handles
        handle

        % Stencils
        mask
        stencil
    end
    
    methods
        function obj = Crystal_growth(params,...
                                      phi_naught,...
                                      T_naught,...
                                      mask,...
                                      phi_bc,...
                                      T_bc)
            arguments
                params
                phi_naught      
                T_naught 
                mask                = struct;
                phi_bc              = @(phi,x) Crystal_growth.bc(phi);
                T_bc                = @(T,x)   Crystal_growth.bc(T);
            end

            %CRYSTAL_GROWTH Constructs an instance of the problem given the
            %   set of parameters contained in params.

            %   Set Allen-Cahn parameters
            obj.M                   = params.M;
            obj.alpha               = params.alpha;
            obj.epsilon             = params.epsilon;
            obj.gamma               = params.gamma;
            obj.Tm                  = params.Tm;

            % Set heat equation parameters
            obj.K                   = params.K;

            % Set discretization
            obj.h           = params.h;
            obj.x           = params.xlim(1):obj.h:params.xlim(2);
            obj.y           = params.ylim(1):obj.h:params.ylim(2);
            obj.sz          = [numel(obj.y) numel(obj.x)];
            [obj.x,obj.y]   = meshgrid(obj.x,obj.y);
            
            
            % Determine stencils
            obj.stencil = Stencil(size(obj.x));
            obj.mask = obj.generate_mask(obj.h,mask);

            % Function handles
            obj.handle.phi.bc           = phi_bc;
            obj.handle.phi.initial      = phi_naught;
            obj.handle.T.bc             = T_bc;
            obj.handle.T.initial        = T_naught;
        end

        function val = get.phi_naught(obj)
            val = obj.handle.phi.initial(obj.x,obj.y);
        end

        function val = get.T_naught(obj)
            val = obj.handle.T.initial(obj.x,obj.y);
        end
        
        function [phi,T] = step(obj,phi,T,k)
            %METHOD1 Summary of this method goes here
            % Create copy of previous phi
            phi_prev = phi; 

            % Take one step in phi field
            phi = phi + k * -obj.M * (obj.df_dphi(phi) - ...
                                      obj.diffusion(phi,...
                                                    coeff=obj.epsilon) +...
                                      obj.m(phi,T));

            % Take one step in Temperature field
            T = T + k * (obj.diffusion(T) + obj.K * (phi-phi_prev)/k);

            % Boundary conditions
            phi     = obj.handle.phi.bc(phi);
            T       = obj.handle.T.bc(T,obj.x);
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

        function res = diffusion(obj,f,args)
            %DIFFUSION calculates the diffusion term of a scalar field. It
            %can either accept scalar fields or scalars as coefficients
            arguments
                obj
                f
                args.coeff = 1;
            end
            
            % Check for coefficient type
            if isnumeric(args.coeff)
                % Constant coefficient case. The diffusion term simply
                % corresponds to the laplacian of the scalar field
                res =  args.coeff^2 * obj.stencil.apply(f,obj.mask.lap);

            elseif isa(args.coeff,'function_handle')
                % Scalar field as coefficient case. The scalar field that
                % will serve as the coefficient must be evaluated first.

                % Evaluate coefficient and its derivative
                [g,g_prime] = args.coeff(f);

                % Evaluate laplacian term
                lap = obj.stencil.apply(f,obj.mask.lap,g.^2);

                % Evaluate cross terms
                %   Evaluate product
                prod    = g .* g_prime;

                term1   = obj.stencil.apply(prod,obj.mask.dx) .*...
                          obj.stencil.apply(f,obj.mask.dy);

                term2   = obj.stencil.apply(prod,obj.mask.dy) .*...
                          obj.stencil.apply(f,obj.mask.dx);
                            
                % Evaluate result
                res = lap - term1 + term2;
            end
        end
    end


    methods (Static)
        function U = bc(U)
            U(1,2:end-1) = U(2,2:end-1);          % North row
            U(end,2:end-1) = U(end-1,2:end-1);    % South row 
            U(1:end,1) = U(1:end,2);              % West row
            U(1:end,end) = U(1:end,end-1);        % East row
        end

        function mask = generate_mask(h,mask)
            %GENERATE_MASK Generates the masks for the stencils to
            %calculate the derivatives of the fields

            % Central difference for the laplacian
            if ~isfield(mask,'lap')
                mask.lap = [0  1  0;
                            1 -4  1;
                            0  1  0];
            end

            % Central difference for the derivatives
            if ~isfield(mask,'dx')
                mask.dx = [0  0  0;
                           -1 0  1;
                           0  0  0];
            end

            if ~isfield(mask,'dy')
                mask.dy = [0 -1  0;
                           0  0  0;
                           0  1  0];
            end

            mask.lap        = mask.lap * (1/h^2); 
            mask.dx         = mask.dx  * (2*h);
            mask.dy         = mask.dy  * (2*h);
        end
    end
end

