%------------------------------------------------------------
% Programmer(s):  Daniel R. Reynolds @ SMU
%------------------------------------------------------------
% Copyright (c) 2018, Southern Methodist University.
% All rights reserved.
% For details, see the LICENSE file.
%------------------------------------------------------------
% driver to check analytical properties of various RK methods
% in butcher.m

clear; close all
!\rm butcher_tests.txt
diary butcher_tests.txt

% test parameters
tol = 1e-9;      % tolerance to use in tests
plotdir = 'stability_regions';    % folder to use for stability region plots
plotregions = false;

% method names, filenames, and stability region bounding boxes to test
tests = { {'DBM-5-3-ERK',               'dbm53_erk',           [-4,1,-5,5]},
          {'Ascher(2,3,3)-ERK',         'a233_erk',            [-3,0.5,-3,3]},
          {'Ascher(2,3,2)-ERK',         'a232_erk',            [-3,0.5,-3,3]},
          {'Ascher(2,2,2)-ERK',         'a222_erk',            [-2.5,0.5,-3,3]},
          {'Ascher(3,4,3)-ERK',         'a343_erk',            [-3,0.5,-4,4]},
          {'Ascher(4,4,3)-ERK',         'a443_erk',            [-3,0.5,-3,3]},
          {'Knoth-Wolke-ERK',           'kw_erk',              [-3,0.5,-3,3]},
          {'3/8-Rule-ERK',              '38rule_erk',          [-3,0.5,-4,4]},
          {'Cooper4-ERK',               'cooper4_erk',         [-3,0.5,-3,3]},
          {'Cooper6-ERK',               'cooper6_erk',         [-3,0.5,-4,4]},
          {'ERK-1-1',                   'erk11',               [-2.5,0.5,-1.5,1.5]},
          {'ERK-2-2',                   'erk22',               [-3,0.5,-3,3]},
          {'ERK-4-4',                   'erk44',               [-3,0.5,-3.5,3.5]},
          {'ARK(2,3,2)-ERK',            'ark232_erk',          [-3,0.5,-3,3]},
          {'SSP2(2,2,2)-ERK',           'ssp2222_erk',         [-2.5,0.5,-3,3]},
          {'SSP2(3,3,2)-lpm1-ERK',      'ssp2332lpm1_erk',     [-5,0.5,-3,3]},
          {'SSP2(3,3,2)-lpm2-ERK',      'ssp2332lpm2_erk',     [-5,0.5,-3,3]},
          {'SSP2(3,3,2)-lpum-ERK',      'ssp2332lpum_erk',     [-5,0.5,-3,3]},
          {'SSP2(3,3,2)-lspum-ERK',     'ssp2332lspum_erk',    [-3.5,0.5,-3,3]},
          {'SSP2(3,3,2)-a-ERK',         'ssp2332a_erk',        [-5,0.5,-3,3]},
          {'SSP2(3,3,2)-b-ERK',         'ssp2332b_erk',        [-5,0.5,-3,3]},
          {'SSP3(3,3,2)-ERK',           'ssp3332_erk',         [-3,0.5,-3,3]},
          {'SSP3(3,3,3)-ERK',           'ssp3333_erk',         [-3,1,-3,3]},
          {'SSP3(4,3,3)-ERK',           'ssp3433_erk',         [-3,0.5,-3,3]},
          {'SSPRK(3,3)-Shu-Osher-ERK',  'ssprk33so_erk',       [-3,0.5,-3,3]},
          {'ARK3(2)4L[2]SA-ERK',        'ark324_erk',          [-4,0.5,-4,4]},
          {'ARK4(3)6L[2]SA-ERK',        'ark436_erk',          [-5,1,-5,5]},
          {'ARK4(3)7L[2]SA-ERK',        'ark437_erk',          [-8,1,-6,6]},
          {'ARK5(4)8L[2]SA-ERK',        'ark548_erk',          [-10,1,-5,5]},
          {'ARK5(4)8L[2]SAb-ERK',       'ark548b_erk',         [-5,1,-4,4]},
          {'Sayfy-Aburub-4-3-ERK',      'sa43_erk',            [-3,0.5,-2.5,2.5]},
          {'Heun-Euler-ERK',            'he_erk',              [-2.5,0.5,-2,2]},
          {'Bogacki-Shampine-ERK',      'bs_erk',              [-3.5,0.5,-3,3]},
          {'Fehlberg-ERK',              'fehlberg_erk',        [-4,1,-4,4]},
          {'Cash-Karp-ERK',             'cashkarp_erk',        [-15,5,-15,15]},
          {'Dormand-Prince-ERK',        'dp_erk',              [-15,5,-15,15]},
          {'ERK-3-3',                   'erk33',               [-3,0.5,-3,3]},
          {'Merson-4-3-ERK',            'm43_erk',             [-4,1,-4,4]},
          {'Zonneveld-4-3-ERK',         'z43_erk',             [-3,0.5,-3,3]},
          {'Verner-6-5-ERK',            'v65_erk',             [-5,1,-4,4]},
          {'Fehlberg-8-7-ERK',          'f87_erk',             [-10,5,-9,9]},
          {'DBM-5-3-ESDIRK',            'dbm53_esdirk',        [-1,9,-6,6]},
          {'Ascher(2,3,3)-SDIRK',       'a233_sdirk',          [-5,15,-10,10]},
          {'Ascher(2,3,2)-SDIRK',       'a232_sdirk',          [-5,15,-10,10]},
          {'Ascher(2,2,2)-SDIRK',       'a222_sdirk',          [-5,15,-10,10]},
          {'Ascher(3,4,3)-SDIRK',       'a343_sdirk',          [-5,10,-6,6]},
          {'Ascher(4,4,3)-SDIRK',       'a443_sdirk',          [-2,10,-5,5]},
          {'Cooper4-ESDIRK',            'cooper4_esdirk',      [-2,15,-7,7]},
          {'Cooper6-ESDIRK',            'cooper6_esdirk',      [-2,10,-6,6]},
          {'SDIRK-2-2',                 'sdirk22',             [-5,15,-8,8]},
          {'ARK(2,3,2)-SDIRK',          'ark232_sdirk',        [-2,15,-10,10]},
          {'SSP2(2,2,2)-SDIRK',         'ssp2222_sdirk',       [-2,15,-10,10]},
          {'SSP2(3,3,2)-lpm1-SDIRK',    'ssp2332lpm1_sdirk',   [-2,30,-20,20]},
          {'SSP2(3,3,2)-lpm2-SDIRK',    'ssp2332lpm2_sdirk',   [-2,30,-20,20]},
          {'SSP2(3,3,2)-lpum-SDIRK',    'ssp2332lpum_sdirk',   [-2,30,-20,20]},
          {'SSP2(3,3,2)-lspum-SDIRK',   'ssp2332lspum_sdirk',  [-2,30,-20,20]},
          {'SSP2(3,3,2)-a-DIRK',        'ssp2332a_dirk',       [-2,15,-10,10]},
          {'SSP2(3,3,2)-b-DIRK',        'ssp2332b_dirk',       [-2,12,-10,10]},
          {'SSP3(3,3,2)-SDIRK',         'ssp3332_sdirk',       [-2,15,-10,10]},
          {'SSP3(3,3,3)-ESDIRK',        'ssp3333_esdirk',      [-4,1,-3,3]},
          {'SSP3(4,3,3)-SDIRK',         'ssp3433_sdirk',       [-5,25,-15,15]},
          {'EDIRK-3-3',                 'edirk33',             [-2,15,-8,8]},
          {'ESDIRK-3-3',                'esdirk33',            [-2,15,-8,8]},
          {'ARK3(2)4L[2]SA-ESDIRK',     'ark324_esdirk',       [0,9,-5,5]},
          {'ARK4(3)6L[2]SA-ESDIRK',     'ark436_esdirk',       [0,30,-20,20]},
          {'ARK4(3)7L[2]SA-ESDIRK',     'ark437_esdirk',       [-5,80,-50,50]},
          {'ARK5(4)8L[2]SA-ESDIRK',     'ark548_esdirk',       [-5,45,-25,25]},
          {'ARK5(4)8L[2]SAb-ESDIRK',    'ark548b_esdirk',      [-5,30,-15,15]},
          {'TRBDF2-ESDIRK',             'trbdf2_esdirk',       [-10,15,-8,8]},
          {'TRX2-ESDIRK',               'trx2_esdirk',         [-10,2,-8,8]},
          {'Billington-SDIRK',          'b_sdirk',             [-30,15,-15,15]},
          {'Cash(5,2,4)-SDIRK',         'c524_sdirk',          [0,90,-50,50]},
          {'Cash(5,3,4)-SDIRK',         'c534_sdirk',          [0,12,-5,5]},
          {'Kvaerno(4,2,3)-ESDIRK',     'k423_esdirk',         [0,90,-50,50]},
          {'Kvaerno(5,3,4)-ESDIRK',     'k534_esdirk',         [0,40,-20,20]},
          {'Kvaerno(7,4,5)-ESDIRK',     'k745_esdirk',         [0,120,-60,60]},
          {'SDIRK-2-1',                 'sdirk21',             [0,6,-3,3]},
          {'SDIRK-5-4',                 'sdirk54',             [-15,25,-15,15]},
          {'IRK-1-1',                   'irk11',               [-0.5,2.5,-1.5,1.5]},
          {'Crank-Nicolson-2-2-IRK',    'cn22_irk',            [-5,5,-5,5]},
          {'SIRK-2-2',                  'sirk22',              [-5,15,-8,8]},
          {'Gauss-2-4-IRK',             'g24_irk',             [-5,5,-5,5]},
          {'RadauIIA-2-3-IRK',          'riia23_irk',          [0,7,-5,5]},
          {'LobattoIII-2-2-IRK',        'liii22_irk',          [-3,1,-3,3]},
          {'LobattoIIIA-2-2-IRK',       'liiia22_irk',         [-5,5,-5,5]},
          {'LobattoIIIC-2-2-IRK',       'liiic22_irk',         [-1,3,-3,3]},
          {'Gauss-3-6-IRK',             'g36_irk',             [-5,5,-5,5]},
          {'RadauI-3-5-IRK',            'ri35_irk',            [-20,5,-10,10]},
          {'RadauIA-3-5-IRK',           'ria35_irk',           [-2,15,-10,10]},
          {'RadauII-3-5-IRK',           'rii35_irk',           [-15,2,-9,9]},
          {'RadauIIA-3-5-IRK',          'riia35_irk',          [-2,15,-10,10]},
          {'LobattoIII-3-4-IRK',        'liii34_irk',          [-8,2,-5,5]},
          {'LobattoIIIA-3-4-IRK',       'liiia34_irk',         [-5,5,-5,5]},
          {'LobattoIIIB-3-4-IRK',       'liiib34_irk',         [-5,5,-5,5]},
          {'LobattoIIIC-3-4-IRK',       'liiic34_irk',         [-2,8,-6,6]},
          {'LobattoIII-4-6-IRK',        'liii46_irk',          [-15,5,-8,8]},
          {'LobattoIIIA-4-6-IRK',       'liiia46_irk',         [-5,5,-5,5]},
          {'LobattoIIIB-4-6-IRK',       'liiib46_irk',         [-5,5,-5,5]},
          {'LobattoIIIC-4-6-IRK',       'liiic46_irk',         [-2,15,-10,10]},
          {'RadauIIA-5-9-IRK',          'riia59_irk',          [-5,30,-20,20]},
          {'LobattoIII-5-8-IRK',        'liii58_irk',          [-20,2,-12,12]},
          {'LobattoIIIA-5-8-IRK',       'liiia58_irk',         [-5,5,-5,5]},
          {'LobattoIIIB-5-8-IRK',       'liiib58_irk',         [-5,5,-5,5]},
          {'LobattoIIIC-5-8-IRK',       'liiic58_irk',         [-2,20,-15,15]},
          {'Gauss-6-12-IRK',            'g612_irk',            [-5,5,-5,5]} };

fprintf('                             |         |       Method        |     Embedding\n');
fprintf('                             |         |                     |\n');
fprintf('            Name             |  s  qs  |  q  lq   A   B   L  |  p  lp   A   B   L\n');
fprintf('-----------------------------------------------------------------------------------\n');
for i = 1:length(tests)
   mname = tests{i}{1};
   fname = [ plotdir, '/', tests{i}{2} ];
   box = tests{i}{3};
   B = butcher(mname);
   s = size(B,2)-1;
   [q,p,qs,lq,lp,Bs,As,Ls,BsE,AsE,LsE] = check_rk(B,tol,0,plotregions,box,mname,fname);
   Bs_ = ' ';
   As_ = ' ';
   Ls_ = ' ';
   BsE_ = ' ';
   AsE_ = ' ';
   LsE_ = ' ';
   if (Bs==1)
      Bs_ = 'Y';
   end
   if (As==1)
      As_ = 'Y';
   end
   if (Ls==1)
      Ls_ = 'Y';
   end
   if (BsE==1)
      BsE_ = 'Y';
   end
   if (AsE==1)
      AsE_ = 'Y';
   end
   if (LsE==1)
      LsE_ = 'Y';
   end
   if (p > 0)
      fprintf(' %26s  | %2i  %2i  | %2i  %2i   %s   %s   %s  | %2i  %2i   %s   %s   %s\n',...
              mname,s,qs,q,lq,As_,Bs_,Ls_,p,lp,AsE_,BsE_,LsE_);
   else
      fprintf(' %26s  | %2i  %2i  | %2i  %2i   %s   %s   %s\n',mname,s,qs,q,lq,As_,Bs_,Ls_);
   end
end
fprintf('---------------------------------------------------------------------------------------\n');



diary off
% end of script
