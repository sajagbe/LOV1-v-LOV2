
### To visualize Electrostatic Projection Map mol2 Output:

select high_charge, partial_charge > 0.25
color blue, high_charge
select low_charge, partial_charge < -0.25
color red, low_charge
select neutral_charge, partial_charge > -0.25 and partial_charge < 0.25
color white, neutral_charge


### To visualize RMSD for full protein residues:

spectrum b, blue_white_red, minimum=0.32, maximum=1.32



### To visualize RMSD for residues 4Å from LF heavy atoms:

spectrum b, blue_white_red, minimum=0.32, maximum=0.73

