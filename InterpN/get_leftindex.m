function x_leftindex = get_leftindex(X,x)
x_leftindex = floor( (x-X(1))/(X(2)-X(1)) ) + 1;