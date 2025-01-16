#CONFIG
 
TITLE!
TMOL!
 
desy 1.d-6
*
no
b "c" 3-21g hondo
b "h" 3-21g hondo
b "o" 3-21g hondo
b "p" 3-21g hondo
b "au" LANL2DZ ECP
*
eht
 
 CHARGE!
 
 
 
 
scf
iter
300

ri
on
m 2000
jbas
b "c" universal
b "h" universal
b "o" universal
b "p" universal
b "au" universal-ecp-60
*
*

marij

 
dft
on
func
b3-lyp
grid
m3
*


*
#END
