%% Input Data
x0 = 10;
y0 = 10;
x_0 = 0;
y_0 = 0;
xf = 50;
yf = 50;
x_f = 0;
y_f = 0;
K0 = 0;
K_0 = 1;
Kf = 0;
K_f = 1;
t0 = 0;
tf = 10;
tc = 5;
xc = 25;
yc = 25;

%% Under Constrained

Dx = [x0; x_0; xf; x_f; xc];
Ax = [getB(t0,0) getB(t0,1) getB(t0,2) getB(t0,3) getB(t0,4) getB(t0,5);
    getBder(t0,0) getBder(t0,1) getBder(t0,2) getBder(t0,3) getBder(t0,4) getBder(t0,5);
    getB(tf,0) getB(tf,1) getB(tf,2) getB(tf,3) getB(tf,4) getB(tf,5);
    getBder(tf,0) getBder(tf,1) getBder(tf,2) getBder(tf,3) getBder(tf,4) getBder(tf,5);
    getB(tc,0) getB(tc,1) getB(tc,2) getB(tc,3) getB(tc,4) getB(tc,5)];
Wx = pinv(Ax)*Dx;

[F00, F01, F02, F03, F04, F05] = get_coeff(Wx(1),Wx(2),Wx(3),Wx(4),Wx(5),Wx(6),t0,t0,tf);
[Ff0, Ff1, Ff2, Ff3, Ff4, Ff5] = get_coeff(Wx(1),Wx(2),Wx(3),Wx(4),Wx(5),Wx(6),t0,tf,tf);
[Fc0, Fc1, Fc2, Fc3, Fc4, Fc5] = get_coeff(Wx(1),Wx(2),Wx(3),Wx(4),Wx(5),Wx(6),t0,tc,tf);

Dy = [K0; yf-y0; Kf; yc-y0 ; K_0 ; K_f];
Ay = [getB(t0,0) getB(t0,1) getB(t0,2) getB(t0,3) getB(t0,4) getB(t0,5);
    Ff0 Ff1 Ff2 Ff3 Ff4 Ff5;
    getB(tf,0) getB(tf,1) getB(tf,2) getB(tf,3) getB(tf,4) getB(tf,5);
    Fc0 Fc1 Fc2 Fc3 Fc4 Fc5;
    getBder(t0,0) getBder(t0,1) getBder(t0,2) getBder(t0,3) getBder(t0,4) getBder(t0,5);
    getBder(tf,0) getBder(tf,1) getBder(tf,2) getBder(tf,3) getBder(tf,4) getBder(tf,5)];
Wk = (eye(6)/Ay)*Dy;

inc = 0.1;
sz = ((tf-t0)/inc) +1;

time = zeros(1,sz);
X = zeros(1,sz);
Y = zeros(1,sz);
theta = zeros(1,sz);
Xdot = zeros(1,sz);
Ydot = zeros(1,sz);

for i=1:sz
    t=(i-1)*inc;
    B = [getB(t,0) getB(t,1) getB(t,2) getB(t,3) getB(t,4) getB(t,5)];
    [Ft0, Ft1, Ft2, Ft3, Ft4, Ft5] = get_coeff(Wx(1),Wx(2),Wx(3),Wx(4),Wx(5),Wx(6),t0,t,tf);
    Bdev = [getBder(t,0) getBder(t,1) getBder(t,2) getBder(t,3) getBder(t,4) getBder(t,5)];
    time(1,i) = t;
    X(1,i) = B(1)*Wx(1)+B(2)*Wx(2)+B(3)*Wx(3)+B(4)*Wx(4)+B(5)*Wx(5)+B(6)*Wx(6);
    Y(1,i) = y0 + Ft0*Wk(1)+Ft1*Wk(2)+Ft2*Wk(3)+Ft3*Wk(4)+Ft4*Wk(5)+Ft5*Wk(6);
    theta(1,i) = atan( B(1)*Wk(1)+B(2)*Wk(2)+B(3)*Wk(3)+B(4)*Wk(4)+B(5)*Wk(5)+B(6)*Wk(6) );
    Xdot(1,i) = Bdev(1)*Wx(1)+Bdev(2)*Wx(2)+Bdev(3)*Wx(3)+Bdev(4)*Wx(4)+Bdev(5)*Wx(5)+Bdev(6)*Wx(6);
    Ydot(1,i) = Xdot(1,i)*tan(theta(1,i));
end

figure;
plot(time,X);
title('x vs t for Under Constrained');
xlabel('time');
ylabel('X');

figure;
plot(time,Y);
title('y vs t for Under Constrained');
xlabel('time');
ylabel('Y');

figure;
plot(time,theta);
title('theta vs t for Under Constrained');
xlabel('time');
ylabel('theta');

figure;
plot(time,Xdot);
title('Xdot vs t for Under Constrained');
xlabel('time');
ylabel('Xdot');

figure;
plot(time,Ydot);
title('Ydot vs t for Under Constrained');
xlabel('time');
ylabel('Ydot');

figure;
plot(X,Y);
title('Y vs X for Under Constrained');
xlabel('X');
ylabel('Y');