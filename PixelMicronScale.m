function [pixelPerMicronScale, rotation_matrix]=PixelMicronScale(masked_image_file)

% 
% 
% 
% 
% 
% 
% 


pixelPerMicronX = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_x');
pixelPerMicronY = 1/h5readatt(masked_image_file, '/mask', 'pixels2microns_y');

normScale = sqrt((pixelPerMicronX ^ 2 + pixelPerMicronX ^ 2) / 2);
pixelPerMicronScale =  normScale * [sign(pixelPerMicronX) sign(pixelPerMicronY)];

% Compute the rotation matrix.
%rotation = 1;
angle = atan(pixelPerMicronY / pixelPerMicronX);
if angle > 0
    angle = pi / 4 - angle;
else
    angle = pi / 4 + angle;
end
cosAngle = cos(angle);
sinAngle = sin(angle);
rotation_matrix = [cosAngle, -sinAngle; sinAngle, cosAngle];