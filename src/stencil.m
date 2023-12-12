function U_new = stencil(U,p,A,args)
arguments
    U
    p
    A
    args.kappa = [];
end

%STENCIL This function calculates a centered stencil specified by the array
%A
    % Allocate memory
    U_new = zeros(size(U));
    % Calculate linear indices of neighboring points
    idx = p + (-1:1)' + (-1:1) * size(U,1);

    % Calculate coefficients
    if isempty(args.kappa)
        kappa = A;
    else
        % Create mask to surrounding points
        mask = A;
        mask(2,2) = 0;

        kappa = mask .* (args.kappa(idx) + args.kappa(p))/2;
        kappa(2,2,:) = -sum(kappa,[1 2]);
    end

    % Update points
    U_new(p) = sum(kappa .* U(idx),[1 2]);
end

