syms x1 real
syms x2 real
syms x3 real
syms y1 real
syms y2 real
syms y3 real
syms u1 real
syms u2 real
syms u3 real
syms x real
syms y real
syms v1 real
syms v2 real
syms v3 real
syms v4 real

A = [1, x1, y1; 1, x2, y2; 1, x3, y3];
A1 = A;
A2 = A;
A3 = A;
A1(:,1) = [u1, u2, u3];
A2(:,2) = [u1, u2, u3];
A3(:,3) = [u1, u2, u3];

detA = det(A);

a = det(A1)/detA;
b = det(A2)/detA;
c = det(A3)/detA;

disp("a = ")
pretty(a)
disp("b = ")
pretty(b)
disp("c = ")
pretty(c)

u = simplify(a + b * x + c * y);

disp("u = ")
pretty(u)

alpha1 = (y3 * x2 - y2 * x3 + (y2 - y3) * x + (x3 - x2) * y) / detA;
alpha2 = (y3 * x1 - y1 * x3 + (y3 - y1) * x + (x1 - x3) * y) / detA;

S_conj = [1, -0.5, 0, -0.5; -0.5, 1, -0.5, 0; 0, -0.5, 1, -0.5; -0.5, 0, -0.5, 1];
U = [v1; v2; v3; v4];
W_square = expand(U' * S_conj * U);
disp("W_square =")
pretty(W_square)
