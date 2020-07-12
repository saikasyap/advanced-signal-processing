
       
        x = [8 8 8 8 8 8 8 8];    % Convert logical to double
        h = [-8 -8 -8 -8 -8 -8 -8 -8];
        y = conv(x,h);
        figure(1)
        plot(y) 