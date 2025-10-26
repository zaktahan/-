function [fwhm, x_interp, field_cut_interp, index1, index2] = interp_FWHM(intensity, x)
% Calculates FWHM of the intensity along central row

N = size(intensity,1);
field_cut = intensity(:,N/2);

x_interp = min(x):0.0001:max(x);
field_cut_interp = interpn(x,field_cut, x_interp);

len = length(field_cut_interp);
field_cut_positive = field_cut_interp(ceil(len/2):end);
field_cut_negative = field_cut_interp(1:ceil(len/2));

width = max(field_cut_interp);

if field_cut_interp(ceil(len/2)) >= width/2
    index1 = find(field_cut_negative <= width/2,1,'last');
    index2 = find(field_cut_positive >= width/2,1,'last');
    fwhm = x_interp(index2+ceil(len/2)-1) - x_interp(index1);
else
    index1 = find(field_cut_interp(ceil(len/2):end) >= width/2,1,'first') + ceil(len/2)-1;
    index2 = find(field_cut_interp(index1+1:end) <= width/2,1,'first') + index1;
    fwhm = x_interp(index2) - x_interp(index1);
end

end
