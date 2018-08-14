function [x_value,x_volumn] = get_leftright(X,x)
    x_leftindex = get_leftindex(X,x);
    step = X(2) - X(1);
    if x == X(1)
        x_value(1) = X(1);
        x_value(2) = X(1);
        x_volumn(1) = X(1);
        x_volumn(2) = X(1) + step;
    elseif x == X(end)
        x_value(1) = X(end);
        x_value(2) = X(end);
        x_volumn(1) = x(end);
        x_volumn(2) = x(end) + step;
    else
        x_value(1) = X(x_leftindex);
        x_value(2) = X(x_leftindex) + step;
        x_volumn(1) = X(x_leftindex);
        x_volumn(2) = X(x_leftindex) + step;
    end