function [sequence] = generateOutputSequenceTXT(path)

sequence = [];
s = path(:,7);
for indi = 1:length(s)
    pl = s(indi);
    if pl == 1
        plSeq = '-M-';
    elseif pl == 2
        plSeq = '-V-';
    elseif pl == 3
        plSeq = '-E-';
    elseif pl == 4
        plSeq = '-M-';
    elseif pl == 5
        plSeq = '-J-';
    elseif pl == 6
        plSeq = '-S-';
    elseif pl == 7
        plSeq = '-U-';
    elseif pl == 8
        plSeq = '-N-';
    end
    sequence = [sequence plSeq];
end

end