function fdraw(X,Y,U,args)
arguments
    X
    Y
    U
    args.title = []
end
%DRAW Draws the solution
    surf(X,Y,U,'EdgeColor','none')
    colorbar
    view(2)
    axis equal
    xlim([0 9])
    ylim([0 9])

    if ~isempty(args.title)
        title(args.title)
    end
end

