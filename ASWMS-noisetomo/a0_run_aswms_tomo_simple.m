% Run ASWMS eikonal tomography using delay times from Bessel fits as the
% input dataset
%

% Convert matnoise Bessel fits to the CSmeasure data structure
a1_matnoise_to_aswms_struct

% Run eikonal tomography treating each virtual source as an event
a6_a_eikonal_eq

% Stack eikonal maps
a7_a_stack_phv

