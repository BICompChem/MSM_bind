use strict;
my $i=1;

while ($i<151) 
{
	my $generate_pdb="echo Protein Protein_MOL | gmx trjconv -f sys$i/traj_comp.xtc -s sys$i/$i.tpr -n index.ndx -pbc mol -skip 5 -o final_trajs/start.pdb -nice 0 -center -e 0\n";
	`$generate_pdb`;

	my $traj1 = "echo Protein Protein_MOL | gmx trjconv -f sys$i/traj_comp.xtc -s sys$i/$i.tpr -n index.ndx -pbc mol -skip 5 -o sys$i/center$i.xtc -nice 0 -center\n";
	`$traj1`;

	my $traj2 = "echo backbone system| gmx trjconv -f sys$i/center$i.xtc -s final_trajs/start.pdb -fit rot+trans -o final_trajs/$i.xtc\n";
	`$traj2`;

	$i++;
}
