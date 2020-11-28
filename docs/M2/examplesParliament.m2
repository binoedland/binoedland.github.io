loadPackage ("Parliament")

-- The following constructs example 2.14 from Di Rocco, Jabbusch, Smith using the packacke ToricVectorBundles
fa = projectiveSpaceFan 2;
A = transpose( matrix (QQ,{{1,0,0},{0,1,0},{0,0,1}}));
B = transpose(  matrix (QQ,{{0,0,1},{0,1,0},{1,0,0}}));
C = transpose( matrix (QQ,{{1,-1,0},{0,-1,1},{1,0,0}}));
L1={C,B,A};
D = matrix{{-4,0,1}};
E = matrix {{-3,0,2}};
F = matrix {{-3,-2,1}};
L2 = {F,E,D};
V = toricVectorBundle(3,fa,L1,L2);
details V
-- We construct the parliament of polytopes
P = parliament V;
-- View the characters at each maximal cone
getCharacters P
-- View the matroid of the parliament
getMatroid P
-- View the polytopes of the parliament
getPolytopes P
-- Check if the vector bundle is globally generated
globallyGenerated P



-- Example 4.7
fa = projectiveSpaceFan 2;
A = transpose( matrix (QQ,{{1,0,0},{0,1,0},{0,0,1}}));
B = transpose(  matrix (QQ,{{0,0,1},{0,1,0},{1,0,0}}));
C = transpose( matrix (QQ,{{1,-1,0},{0,-1,1},{1,0,0}}));
L1={C,B,A};
D = matrix{{-2,1,2}};
E = matrix {{-2,0,2}};
F = matrix {{-4,-3,-1}};
L2 = {F,E,D};
V = toricVectorBundle(3,fa,L1,L2);
details V;
-- We construct the parliament of polytopes
P = parliament V;
-- View the characters at each maximal cone
getCharacters P
-- View the matroid of the parliament
getMatroid P
-- View the polytopes of the parliament
getPolytopes P
-- Check if the vector bundle is globally generated
globallyGenerated P


 -- Example 2.9
PP = pp1ProductFan 2
A = matrix (QQ,{{0,0,1},{1,0,0},{0,1,0}});
B = matrix (QQ,{{0,0,1},{1,0,0},{0,1,0}});
C = matrix (QQ,{{1,1,0},{0,0,1},{1,0,0}});
D = matrix (QQ,{{1,0,0},{1,1,0},{0,0,1}});
L1={A,B,C,D};
E = matrix{{-2,-1,0}};
F = matrix {{-2,-1,0}};
G = matrix {{-1,0,1}};
H = matrix {{-1,0,1}};
L2 = {E,F,G,H};
V = toricVectorBundle(3,PP,L1,L2);
details V
-- We construct the parliament of polytopes
P = parliament V;
-- View the characters at each maximal cone
getCharacters P
-- View the matroid of the parliament
getMatroid P
-- View the polytopes of the parliament
getPolytopes P
-- Check if the vector bundle is globally generated
globallyGenerated P


---Example 2.16
Sigma = hirzebruchFan 1;
A = transpose( matrix (QQ,{{1,0},{0,1}}));
B = transpose (matrix (QQ,{{1,0},{0,1}}));
C = transpose (matrix (QQ,{{0,1},{1,0}}));
D = transpose (matrix (QQ,{{1,1},{1,0}}));
L1={D,C,B,A};
E = matrix{{-4,2}};
F = matrix {{-3,-2}};
G = matrix {{-5,0}};
H = matrix {{-3,1}};
L2 = {H,G,F,E};
V = toricVectorBundle(2,Sigma,L1,L2);
details V
-- We construct the parliament of polytopes
P = parliament V;
-- View the characters at each maximal cone
getCharacters P
-- View the matroid of the parliament
getMatroid P
-- View the polytopes of the parliament
getPolytopes P
-- Check if the vector bundle is globally generated
globallyGenerated P
