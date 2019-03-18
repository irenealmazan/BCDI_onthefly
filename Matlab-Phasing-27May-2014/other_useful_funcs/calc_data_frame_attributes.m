function [attributes ] = calc_data_frame_attributes(data)
%jclark
%calc com, ma, sum etc from a data frame

attributes.max_val=max(data(:)); 

[yy,xx]=ind2sub(size(data),find( squeeze(abs(data)) == ((max(abs(data(:)))))));

attributes.max_x=xx;
attributes.max_y=yy;

com=center_of_mass(data);
attributes.cent_of_mass_x=com(1);
attributes.cent_of_mass_y=com(2); 

attributes.sum_data=sum(data(:));

end

