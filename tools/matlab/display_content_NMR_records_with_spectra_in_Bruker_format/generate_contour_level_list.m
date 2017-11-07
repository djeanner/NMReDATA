function [cont_level_list]=generate_contour_level_list(from_level,to_level,factor_level)

nb_level=log(to_level/from_level)/(log(factor_level));
cont_level_list=from_level*power(factor_level,0:1:nb_level);

end