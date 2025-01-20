% index 1-3 are the original system
% index 4-12 are the linearized evolution equation for the reference framework
function yp=lorenz(t,y)
    sigma=16;r=45.92;b=4;
    yp=zeros(12,1);
    yp(1) = sigma*(y(2)-y(1));
    yp(2) = (r-y(3))*y(1)-y(2);
    yp(3) = y(2)*y(1)-b*y(3);
    for i=0:2
        yp(4+i)=sigma*(y(7+i)-y(4+i));
        yp(7+i)=(r-y(3))*y(4+i)-y(7+i)-y(1)*y(10+i);
        yp(10+i)=y(2)*y(4+i)+y(1)*y(7+i)-b*y(10+i);
    end
end
