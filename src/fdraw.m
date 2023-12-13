function fdraw(X,Y,U,args)
arguments
    X
    Y
    U
    args.title = []
    args.bar_lim = [0 1]
    args.bar_ticks = [];
end
%DRAW Draws the solution
    surf(X,Y,U,'EdgeColor','none')
    
    view(2)
    axis equal
    c = colorbar;
    xlim([0 9])
    ylim([0 9])
    xlabel('x')
    ylabel('y')
    
    if ~isempty(args.bar_lim)
        clim(args.bar_lim)
    end

    if ~isempty(args.bar_ticks)
        ticks = linspace(c.Limits(1),c.Limits(2),length(args.bar_ticks));
        c.Ticks = ticks;
        c.TickLabels = args.bar_ticks;
    end

    if ~isempty(args.title)
        title(args.title)
    end
end

