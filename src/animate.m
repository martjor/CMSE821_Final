function [phi,T] = animate(simulation,t_span,k,options)
arguments
    simulation
    t_span
    k
    options.draw_rate = 0.05;
end
%ANIMATE Creates an animation of the evolving crystal
    % Discretize time
    t = t_span(1):k:t_span(2);
    n = length(t);
    draw_rate = floor(n * options.draw_rate);

    for i = 1:length(t)
        if i == 1
            phi = simulation.phi_naught;
            T = simulation.T_naught;
        else
            [phi,T] = simulation.step(phi,T,k);
        end
        
        if mod(i-1,draw_rate) == 0
            subplot(1,2,1)
            fdraw(simulation.x,simulation.y,phi,"title",sprintf("t=%f",t(i)))

            subplot(1,2,2)
            fdraw(simulation.x,simulation.y,T,"title",sprintf("t=%f",t(i)))
            pause(0.1)
            drawnow
        end
    end
end

