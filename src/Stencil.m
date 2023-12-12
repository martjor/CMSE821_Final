classdef Stencil
    %STENCIL performs stencil-like operations on arrays
    properties
        idx_p
        idx_neigh
    end
    methods
        function obj = Stencil(sz)
            obj.idx_p       = Stencil.index_points(sz);
            obj.idx_neigh   = Stencil.index_neighbors(obj.idx_p,sz(1));
        end

        function U_new = apply(obj,U,A,kappa,idx_p,idx_neigh)
            arguments
                obj
                U
                A
                kappa      = [];
                idx_p      = obj.idx_p;
                idx_neigh  = obj.idx_neigh;
            end

            % Allocate memory
            U_new = zeros(size(U));

            % Calculate coefficients
            if isempty(kappa)
                kappa = A;
            else
                % Create mask to surrounding points
                mask = A;
                mask(2,2) = 0;
        
                kappa = mask .* (kappa(idx_neigh) + kappa(idx_p))/2;
                kappa(2,2,:) = -sum(kappa,[1 2]);
            end
        
            % Update points
            U_new(idx_p) = sum(kappa .* U(idx_neigh),[1 2]);
        end
    end

    methods (Static)
        function idx_p = index_points(sz)
            % INDEX_POINTS calculate the linear indices of the points that
            % are accessible by the stencil

            % Generate indices to accessible points
            idx_p = reshape((2:sz(1)-1)' + (1:sz(2)-2) * sz(1),1,1,[]);
        end

        function idx_neigh = index_neighbors(idx_p,m)
            % INDEX_NEIGHBORS calculates the linear indices of the
            % neighboring points 

            % Calculate linear indices of neighboring points
            idx_neigh = idx_p + (-1:1)' + (-1:1) * m;
        end
    end
end

