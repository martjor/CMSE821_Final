function [phi,T] = animate(simulation,t_span,k,options)
arguments
    simulation
    t_span
    k
    options.draw_rate = 0.05;
end
%ANIMATE Creates an animation of the evolving crystal
    % Discretize time
    t_span = t_span(1):k:t_span(2);
    n = length(t_span);
    draw_rate = floor(n * options.draw_rate);

    for i = 1:length(t_span)
        if i == 1
            phi = simulation.phi_naught;
            T = simulation.T_naught;
        else
            [phi,T] = simulation.step(phi,T,k);
        end
        
        if mod(i-1,draw_rate) == 0   
            % Create tiled layout
            if i == 1
                tiledlayout('horizontal');
                
                ax_phi = nexttile;
                fdraw(simulation.x,simulation.y,phi,"title","\phi")
                colormap(ax_phi,"parula")
     
                ax_T = nexttile;
                fdraw(simulation.x,simulation.y,T,"title","T")
                ylabel("")
                yticklabels([])
                colormap(ax_T,"hot")

                t = text(ax_phi,ax_phi.XLim(2) * 0.10,...
                                ax_phi.YLim(2) * 0.90,...
                                sprintf("t=%.2f",t_span(i)),...
                                'Color','m',...
                                'FontSize',12);
            end

            % Update phi
            set(ax_phi.Children(2),"ZData",phi);

            % Update T
            set(ax_T.Children,"ZData",T);

            % Update label
            t.String = sprintf("t=%.2f",t_span(i));

            drawnow limitrate         
        end
    end
end

