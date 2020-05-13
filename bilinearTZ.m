function [zz, pp, wzm] = bilinearTZ(z, p , wzm, fpr)
    % Transform analog filter to discrete filter
    pp = []; zz = [];
    for k=1:length(z)
       zz = [ zz (2*fpr+z(k))/(2*fpr-z(k)) ];
       wzm = wzm*(2*fpr-z(k));
    end
    for k=1:length(p)
        pp = [ pp (2*fpr+p(k))/(2*fpr-p(k)) ];
        wzm = wzm/(2*fpr-p(k));
    end
    l1 = length(p) - length(z);
    l2 = length(z) - length(p);
    if (l1 > 0) zz = [ zz -1*ones(1, l1) ]; end
    if (l2 > 0) pp = [ pp -1*ones(1, l2) ]; end
end

