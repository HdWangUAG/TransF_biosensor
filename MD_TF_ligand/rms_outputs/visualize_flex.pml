# PyMOL visualization script (RMSF in B-factors)
load annotated_rmsf_bfactor.pdb, prot
hide everything, prot
show cartoon, prot
spectrum b, blue_white_red, prot
cartoon putty, prot
set putty_scale_min, 0
set putty_scale_max, 4.556
set putty_radius, 0.3
# Flexible cutoff (mean+SD) = 4.126 Å
select flexible, prot and b > 4.126
color red, flexible
show sticks, flexible
set cartoon_highlight_color, red
zoom prot
bg_color white
set ray_opaque_background, off
# ray 1600,1200; png flexible_regions.png, dpi=300
