function area_distortion(F, V, uv)
Area = FaceArea(F, V);
uvArea = FaceArea(F, uv);
AreaQuality = (uvArea/sum(uvArea)) ./ (Area/sum(Area));
fprintf('----------------------------------\n')
fprintf('mean of area distortion = %.2e\n', mean(abs(AreaQuality - 1)) );
fprintf('SD   of area distortion = %.2e\n',  std(abs(AreaQuality - 1)) );
fprintf('----------------------------------\n')
end


function Area = FaceArea(F, V)
E12 = V(F(:,2),:) - V(F(:,1),:);
E23 = V(F(:,3),:) - V(F(:,2),:);
E31 = V(F(:,1),:) - V(F(:,3),:);
Len_E12 = sqrt(sum(E12.^2, 2));
Len_E23 = sqrt(sum(E23.^2, 2));
Len_E31 = sqrt(sum(E31.^2, 2));
s = 0.5*(Len_E12+Len_E23+Len_E31);
Area = sqrt(s.*(s-Len_E12).*(s-Len_E23).*(s-Len_E31));
end

