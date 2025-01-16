#------------------------------------------------------------------------------------------------------------------------------------------
#ATOMBOND, ATOMBONDS
#Defining unusual atom bonds
#==========================================================================================
#Iron, IRON, FE, Fe, fe
unbond /////Fe,/////C
distance Bond_FeC = (elem Fe), (elem C), 2.5
color Brown, Bond_FeC
#==========================================================================================
#Unbond In, Ga, Zn, Al, etc with O for amorphous solids
select InL25, elem O within 2.5 of name In
unbond /////In,InL25
delete InL25
select SiL25, elem O within 2.5 of name Si
unbond /////Si,SiL25
delete SiL25
select SnL25, elem O within 2.5 of name Sn
unbond /////Sn,SnL25
delete SnL25
select ZnL25, elem O within 2.5 of name Zn
unbond /////Zn,ZnL25
delete ZnL25
select GaL25, elem O within 2.5 of name Ga
unbond /////Ga,GaL25
delete GaL25
#select AlL25, elem O within 2.5 of name Al
#unbond /////Al,AlL25
#delete AlL25
#==========================================================================================
#Rhodium, RHODIUM, RH, Rh, rh
bond /////Rh,/////P
unbond /////Rh,/////C
distance Bond_RhH = (/////Rh), (elem H), 2.5; color Grape, Bond_RhH
distance Bond_RhC = (/////Rh), (elem C), 2.6; color Grape, Bond_RhC
distance Bond_RhN = (/////Rh), (elem N), 2.5; color Grape, Bond_RhN
distance Bond_RhO = (/////Rh), (elem O), 2.5; color Grape, Bond_RhO
distance Bond_RhF = (elem Rh), (elem F), 3.0; color Grape, Bond_RhF
distance Bond_RhCl = (elem Rh), (elem Cl), 3.0; color Grape, Bond_RhCl
distance Bond_RhBr = (elem Rh), (elem Br), 3.0; color Grape, Bond_RhBr
#select RhL26, (elem C or elem Cl or elem P) within 2.4 of name Rh
#bond /////Rh,RhL26
#delete RhL26
#select all-original, *_high
#select Rh, elem Rh and all-original
#unbond Rh, all
#delete Rh
#delete all-original
#------------------------------------------------------------------------------------------------------------------------------------------
#==========================================================================================
#Palladium, PD, Pd, pd
bond /////Pd,/////P
unbond /////Pd,/////C
distance Bond_PdH = (/////Pd), (elem H), 2.2; color Grape, Bond_PdH
distance Bond_PdC = (/////Pd), (elem C), 2.4; color Grape, Bond_PdC
distance Bond_PdN = (/////Pd), (elem N), 2.4; color Grape, Bond_PdN
distance Bond_PdO = (/////Pd), (elem O), 2.4; color Grape, Bond_PdO
distance Bond_PdF = (elem Pd), (elem F), 3.0; color Grape, Bond_PdF
distance Bond_PdCl = (elem Pd), (elem Cl), 3.0; color Grape, Bond_PdCl
distance Bond_PdBr = (elem Pd), (elem Br), 3.0; color Grape, Bond_PdBr
select PdL26, (elem C) within 2.2 of name Pd
bond /////Pd,PdL26
delete PdL26
select PdL26, (elem Br) within 3.0 of name Pd
bond /////Pd,PdL26
delete PdL26
#==========================================================================================
#Aluminum, Aluminium, AL, Al, al
#select AlL30, (elem O or elem Cl or elem N) within 1.7 of name Al
#bond /////Al,AlL30
distance Bond_AlO = (name Al), (/////O), 2.0; color Eggplant, Bond_AlO
distance Bond_AlCl = (name Al), (/////Cl), 1.7; color Eggplant, Bond_AlCl
delete AlL30
#==========================================================================================
#Gallium, Gauminium, GA, Ga, ga
#select GaL30, (elem O or elem Cl or elem N) within 2.2 of name Ga
#bond /////Ga,GaL30
#distance Bond_GaO = (name Ga), (/////O), 2.2; color Eggplant, Bond_GaO
#distance Bond_GaCl = (name Ga), (/////Cl), 2.2; color Eggplant, Bond_GaCl
#delete GaL30

#select OH12, elem H within 1.2 of name O
#bond /////O,OH12
#select OO2, elem O within 2.0 of name O
#unbond /////O,/////O

#------------------------------------------------------------------------------------------------------------------------------------------
#SPHERE, SPHERES, STICK, STICKS, LABEL, LABELS, GAP, GAPS, TEXTURE, TEXTURES
#Setting molecular representations (bonds, atoms, representation, label, dash)
#==========================================================================================
#sphere (spheres)
show sphere, *_med
show sphere, *_high
show sphere, *_orig
set sphere_scale, 0.25
set sphere_scale, 1.00, *_orig
set sphere_transparency, 0.5, *_orig
disable *_orig
show stick, *_med
show stick, *_high
set stick_radius=0.15, *_high
#===============================================================================
#Low Level settings
set lines, *_low
#===============================================================================
#dat Level settings
show sticks, *_dat
hide lines, *_dat
set stick_radius, 0.05, *_dat
set stick_color, black, *_dat
set stick_transparency, 0.5, *_dat
#==========================================================================================
#line (lines)
set line_smooth, 1
#==========================================================================================
#label (labels)
set label_font_id, 13
#set label_position (0,0,5)
#set texture_fonts=1
set label_size= 26
set label_color, Black
set label_outline_color, White
set label_distance_digits, 2
set label_angle_digits, 1
set label_dihedral_digits, 1
#set float_labels, 1
set label_shadow_mode, 0
#------------------------------------------------------------------------------------------------------------------------------------------
#unactivating, disactivating objects and groups
disable Distance_*
disable ALL_*
disable Steric_*
disable ESP_*
disable TMBond_*
disable BOND_*
disable Classic_ESP_*
delete XH_32
