% Works to obtain the Euclidean Distance of the passed vector

% Executes on being called, with input vectors X and Y.

function value = euclideanDistance(X, Y)

[r, c] = size(X);       % The length of the vector...

e = [];

% Euclidean Distance = sqrt [ (x1-y1)^2 + (x2-y2)^2 + (x3-y3)^2 ... soon on]

for i = 1:c
    e(i) = (X(i)-Y(i))^2;
end

Euclid = sqrt(sum(e));

%Obtain energyLevel...
value = Euclid;


