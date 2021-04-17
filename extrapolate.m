function fun = extrapolate(x1,x2,y1,y2, x_or_y)
    m = (y2 - y1) / (x2 - x1);
    b =  (y1 * x2 - y2 * x1) / (x2 - x1);
    if x_or_y == 'y'
        fun = @(x) (m * x + b); 
    else
        fun = @(y) ( (y - b) / m );
    end
 
end