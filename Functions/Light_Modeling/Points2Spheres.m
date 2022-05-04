function OutVol=Points2Spheres(InVol,radius)

% This function take a binary volume that contains points of ones and
% outputs a binary volume that contains spheres centered at the point
% locations each with radius radius.  If the image is not binary, it is
% made binary with no threshold.

InVol(InVol~=0)=1;
dims=size(InVol);
[Xs,Ys,Zs]=ind2sub(size(InVol),find(InVol==1));
OutVol=InVol;

for s=1:length(Xs)
for x=(max([floor(Xs(s)-radius),1])):min([dims(1),ceil(Xs(s)+radius)])
    for y=(max([floor(Ys(s)-radius),1])):min([dims(2),ceil(Ys(s)+radius)])
        for z=(max([floor(Zs(s)-radius),1])):min([dims(3),ceil(Zs(s)+radius)])
            if (((x-Xs(s))^2 + (y-Ys(s))^2 + (z-Zs(s))^2)<=radius^2)
                OutVol(x,y,z)=1;
            end
        end
    end
end
end