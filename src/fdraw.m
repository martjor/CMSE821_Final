function fdraw(X,Y,U,args)
arguments
    X
    Y
    U
    args.title = []
    args.color_label = []
end
%DRAW Draws the solution
    surf(X,Y,U,'EdgeColor','none')
    c = colorbar;
    view(2)
    axis equal
    xlim([0 9])
    ylim([0 9])
    xlabel('x')
    ylabel('y')
    clim([0 1])

    if ~isempty(args.color_label)
        c.Label.String = args.color_label;
    end

    if ~isempty(args.title)
        title(args.title)
    end
end

