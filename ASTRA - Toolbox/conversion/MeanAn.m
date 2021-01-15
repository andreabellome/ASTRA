function [MPL] = MeanAn(rrPL)

xPL = rrPL(1);
yPL = rrPL(2);

MPL = atan2(yPL, xPL);
MPL = wrapTo2Pi(MPL);

end