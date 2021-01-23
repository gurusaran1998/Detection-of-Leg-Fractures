

m = zeros(300);
[row,col] = size(m);
r = 5;
m = zeros(round(ceil(max(row,col)/2/(r+1))*3*(r+1)));
siz = size(m,1);
sx = round(siz/2);
i = 1:round(siz/2/(r+1));
j = 1:round(0.9*siz/2/(r+1));
j = j-round(median(j));
m(sx+2*j*(r+1),(2*i-1)*(r+1)) = 1;
se = strel('disk',r);
m = imdilate(m,se);
m = m(round(siz/2-row/2-6):round(siz/2-row/2-6)+row-1,round(siz/2-col/2-6):round(siz/2-col/2-6)+col-1);