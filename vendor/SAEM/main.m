%% Spherical Authalic Energy Minimization
%
% Main program:
% [S] = SphericalAEM(F, V)
% 
% Required Input:
% F: #F x 3 triangulations of a closed triangle mesh
% V: #V x 3 vertex coordinates of a closed triangle mesh
%
% Output:
% S: #V x 3 vertex coordinates of the spherical area-preserving map
%
% Optional Input:
% SphericalAEM( __, "MaxIter", Value): the maximum iterative number (default: 100)
% SphericalAEM( __, "Tol", Value): the tolerance of stoping criteria (default: 1e-5)
%
%
%% Riemannian Bijective Correction
%
% Main prograom:
% [S] = RiemanianBijectiveCorrection(F, S)
%
% Required Input:
% F: #F x 3 triangulations of a closed triangle mesh
% S: #V x 3 vertex coordinates of the spherical map
%
% Output:
% S: #V x 3 vertex coordinates of the spherical bijective map
%
%
%% Shape Description by Spherical Harmonic
%
% Main prograom:
% [SH, coef, ReV] = ShapeDescription(V, S, L)
%
% Required Input:
% V: #V x 3 vertex coordinates of a closed triangle mesh
% S: #V x 3 vertex coordinates of the spherical map
% L: the maximum degree of spherical harmonics
%
% Output:
% SH: #V x (L+1)^2 spherical harmonics of spherical map
% coef: (L+1)^2 x 3 coefficient of spherical harmonics
% Rev: #V x 3 reconstructed vertex coordinates by spherical coordinate
%
%
%% Remark:
% If you use this code in your own work, please cite the following paper:
% [1] S.-Y. Liu, and M.-H. Yueh, "Spherical Area-Preserving Parameterization 
%       via Energy Minimization"
% doi: 10.1137/25M1736979
%
% License:
% This software is released for academic and research purposes only.
% Commercial use is not permitted without prior written permission from the authors.
% Copyright (c), Shu-Yung Liu and Mei-Heng Yueh
%
%
%% Example: Area-preserving map for DavidHead
clear; clc

load('DavidHead.mat');

S = SphericalAEM(F, V);
area_distortion(F, V, S)

plot_mesh(F, V);
title('Original Surface');

plot_mesh(F, S);
title('Area-Preserving Map');



%% Example: Unfold the RightHand
load('RightHand_Folding.mat');

folding_plot_mesh(F, S);
title('Folding Map')

S = RiemanianBijectiveCorrection(F, S);

folding_plot_mesh(F, S);
title('Bijective Map')



%% Example: Shape Description
load('DavidHead_map.mat');

d = 30;
[SH, coef, ReV] = ShapeDescription(V, S, d);

plot_mesh(F, ReV);
title('Reconstructed Surface');