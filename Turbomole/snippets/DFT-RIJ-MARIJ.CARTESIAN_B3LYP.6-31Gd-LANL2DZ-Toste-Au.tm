#CONFIG
 
TITLE!
TMOL!
 
desy 1.d-6
*
no
b "c" 6-31G*
b "h" 6-31G*
b "o" 6-31G*
b "p" 6-31G*
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
