function [vec_x,vec_y] = rotate_vector(vec_1,vec_2,theta)
%		[vec_x,vec_y] = rotate_vector(vec_1,vec_2,theta)	
% 		Represent vec_1 and vec_2 in a coordinate system that is 
%        rotated COUNTER-CLOCKWISE FOR THETA DEGREES
%		vec_y,vec_x, two orthoganal components of vector after the rotation
%		vec_1 and vec_2: arbitary orthoganal component of OBS rotated
%		theta counter-clockwise from H1 to E.
%      theta is in degree
%   THIS IS A ROTATION OF COORDINATE SYSTEM NOT VECTOR ITSELF
r_theta = theta*pi/180;    
vec_x = cos(r_theta) * vec_1 + sin(r_theta) * vec_2;
vec_y = -sin(r_theta) * vec_1 + cos(r_theta) * vec_2;


return
end
