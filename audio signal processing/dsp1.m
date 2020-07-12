
span = 2.5; npts = 100;
x = linspace(0,span,npts); % creates npts values in the interval [0,span]
y = besseli(0,x); % evaluates Bessell of order 0, at x values
plot(x,y); % plot results
X = linspace(0,span,npts);
ylim('manual')
ylim([0,3.5])
grid on
title('Zeroth Order Bessel')
xlabel('Argument')
ylabel('Value')