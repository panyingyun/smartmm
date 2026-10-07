function folding_plot_mesh(F, S)
Vno = size(S, 1);
e = ones(Vno, 1);
RGB = [129/255, 159/255, 247/255];
Vrgb = e*RGB;
figure
patch('Faces', F, 'Vertices', S, 'FaceVertexCData', Vrgb,...
    'EdgeColor', 'k', 'EdgeAlpha', 0.25, 'FaceColor', 'interp');
hold on

[FoldingNum, FoldingInd] = FoldingSphere(F, S);
if (FoldingNum > 0)
    patch('Faces', F(FoldingInd,:), 'Vertices', S, 'FaceColor', 'r',...
    'EdgeColor', 'r', 'EdgeAlpha', 1, 'FaceAlpha', 1);
end
view([-1, 1, 0.5])
camlight
axis equal off
end




function [FoldingNum, FoldingInd] = FoldingSphere(F, S)
CF = (S(F(:,1),:) + S(F(:,2),:) + S(F(:,3),:) ) / 3;
CF = CF ./ vecnorm(CF,2,2);
NF = FaceNormal(F, S);
IP = dot(CF, NF, 2);
FoldingInd = IP<0;
FoldingNum = sum(FoldingInd);
Fno = size(F,1);
if FoldingNum > Fno/2
    FoldingNum = Fno - FoldingNum;
    FoldingInd = ~FoldingInd;
end
FoldingInd = find(FoldingInd);
end

function NF = FaceNormal(F, V)
E12 = V(F(:,2),:) - V(F(:,1),:);
E13 = V(F(:,3),:) - V(F(:,1),:);
NF = cross(E12, E13);
NF = NF ./ vecnorm(NF,2,2);
end