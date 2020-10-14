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
