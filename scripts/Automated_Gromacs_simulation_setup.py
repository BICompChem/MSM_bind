# Automated Gromacs Setup #
# Author: Aniket Magarkar #
# Version: 1.0 #
# email: aniket_suresh.magarkar@boehringer-ingelheim.com #

############################################################################################################################
# This script has following dependencies 
# 1. Gromacs 5.1.* : https://doi.org/10.1016/j.softx.2015.06.001
# 2. Acpype : https://doi.org/10.1186/1756-0500-5-367
# 3. VMD : https://doi.org/10.1016/0263-7855(96)00018-5
# 4. RDKIT : http://www.rdkit.org/
# 5. Accesory files (MDP files, and tcl files, provided here)
# 6. The output generated from this script can be used as input for MSM_bind_script 
############################################################################################################################
#
#	This script performs following tasks 
#	1. Extracts protein and crystallized water from ptn.pdb file and generates topology
#	2. Calculates charge of ligand molecule and generates GAFF topology with Acpype
#	3. Combines protein and ligand topologies, the system is hydrated, saline salt concentration is added, additional system charges are neutralized
# 	4. 150 replicates of the system are created, all the systems are then minimized, equilibrated and submitted for production run of 200ns
#	Note: Ligand is placed randomly around protein independantly in in 150 replicates
#
############################################################################################################################
#
# The script takes following 3 inputs 
# 1. Working directory : Path to directory where all files will be generated
# 2. Input protein file path: PDB file, named ptn.pdb 
# 3. Input ligand file path:  mol2 file, lig.mol2
#
############################################################################################################################

print '############################################################################################################################\n'
print '# 	This script has following dependencies '
print '# 	1. Gromacs 5.1.* : https://doi.org/10.1016/j.softx.2015.06.001'
print '# 	2. Acpype : https://doi.org/10.1186/1756-0500-5-367'
print '# 	3. VMD : https://doi.org/10.1016/0263-7855(96)00018-5'
print '# 	4. RDKIT : http://www.rdkit.org/'
print '# 	5. Accesory files (MDP files, and tcl files, provided here)'
print '# 	6. The output generated from this script can be used as input for MSM_bind_script '
print '############################################################################################################################'
print '#'
print '#	This script performs following tasks '
print '#	1. Extracts protein and crystallized water from ptn.pdb file and generates topology'
print '#	2. Calculates charge of ligand molecule and generates GAFF topology with Acpype'
print '#	3. Combines protein and ligand topologies, the system is hydrated, saline salt concentration is added, additional system charges are neutralized'
print '# 	4. 150 replicates of the system are created, all the systems are then minimized, equilibrated and submitted for production run of 200ns'
print '#	Note: Ligand is placed randomly around protein independantly in in 150 replicates'
print '#'
print '############################################################################################################################'

import os
import commands
import string
from rdkit import Chem
from random import *
import re
import time

my_working_dir =input("Enter the working directory path: ")
my_protein_path=input("Enter path for protein input file 'ptn.pdb': ")
my_ligand_path =input("Enter path for ligand input file 'lig.mol2': ")
job_submit_command =input("Enter command for submitting jobs to your cluster :")
# User input dir
#my_working_dir='/home/magarkar/Desktop/negative_data_simulation/sim3/all_input2'
# User input protein
#my_protein_path='/home/magarkar/Desktop/negative_data_simulation/sim3/all_input2/protein/ptn.pdb'
# User input ligand
#my_ligand_path='/home/magarkar/Desktop/negative_data_simulation/sim3/all_input2/ligand/lig.mol2'

def clean_up(my_working_dir):
    print "cleaning up"
    os.chdir(my_working_dir)
    #print called
    commands.getoutput("rm system/*#*")
    commands.getoutput("rm work_protein/*#*")
    commands.getoutput("rm work_ligand/*#*")

def combine_topologies():
    
    commands.getoutput('cp ../work_ligand/MOL.itp .')
    commands.getoutput('cp ../work_ligand/MOL_GMX.gro .')

    fr=open("topol2.top","w")
    with open("topol.top","r") as all_lines:
        for line in all_lines:
            a=re.search('/forcefield',line)
            b=re.search('ions.itp',line)
            c=re.search('Protein_chain_X     1',line)
            if a:
                fr.write(line)
                fr.write('#include "ffMOL.itp"')
                fr.write('\n')
            elif b:
                fr.write(line)
                fr.write('#include "MOL.itp"')
            elif c:
                fr.write(line)
                fr.write('MOL                 1\n\n')
            else:
                fr.write(line)
    commands.getoutput('mv topol.top backup_topol.top')
    commands.getoutput('mv topol2.top topol.top')


def equilibrate(my_working_dir,my_protein_path):
    os.chdir(my_working_dir)
    commands.getoutput('rm -fr eq')
    os.mkdir('eq')
    os.chdir('eq')
    
    command_string='cp '+ my_protein_path + ' ' + 'user_ptn.pdb'
    commands.getoutput(command_string)
    
    commands.getoutput('cp ../scripts/process_protein.tcl .')
    commands.getoutput('cp ../scripts/em.mdp .')
    commands.getoutput('cp ../scripts/eq_vsite_md.mdp .')
        
    extract_protein="vmd -dispdev text user_ptn.pdb -e process_protein.tcl"
    commands.getoutput(extract_protein)
    
    generate_protein_topology='gmx pdb2gmx -f %s -ff amber99sb-ildn -water spce' % ('ptn2.pdb')
    commands.getoutput(generate_protein_topology)
    
    box = 'gmx editconf -f %s -o %s -d 0.8 -bt cubic' % ('conf.gro','box.gro')
    commands.getoutput(box)
    
    # Protein charge calculation
    commands.getoutput("rm all_charges")
    commands.getoutput("egrep \" q \" topol.top > all_charges")
    total_charge=0.0
    with open("all_charges","r") as all_lines:
        for line in all_lines:
                a=line
                b=a.split()
                total_charge = total_charge + float(b[-1])
    
    initial_na_ions=0
    initial_cl_ions=0
    
    if (total_charge>0):
        initial_cl_ions=int(total_charge)
    elif (total_charge>0):
        initial_na_ions=abs(int(total_charge))
    
    solvate='gmx solvate -cp %s -cs -p -o conf.gro' % ('box.gro')
    commands.getoutput(solvate)
    
    tpr_gen1='gmx grompp -f em.mdp'
    commands.getoutput(tpr_gen1)
    
    ionize='echo SOL| gmx genion -s -np %s -nn %s -p -o ions.gro'%(initial_na_ions,initial_cl_ions)
    commands.getoutput(ionize)
    
    tpr_gen2='gmx grompp -f em.mdp -c ions.gro'
    commands.getoutput(tpr_gen2)
    commands.getoutput("gmx mdrun")
    
    commands.getoutput("rm *#*")

    print "Submitting a Gromacs 5 job on the GPU queue\n"
    print "Running equilibration for 100 ps (default)..............."
    
    tpr_gen3='gmx grompp -f eq_vsite_md.mdp -c confout.gro -maxwarn 1'
    commands.getoutput(tpr_gen3)
    
    current_dir=os.getcwd()
    my_tpr=current_dir+'/topol.tpr'
    
    commands.getoutput("rm -fr *.gro")
    
    run_on_gpu="job_submit_command -s %s " % (my_tpr)
    commands.getoutput(run_on_gpu)
    
    check=1
    
    while check>0:
        #print "looking for confout.gro"
        out=commands.getoutput("ls *.gro")
        look=re.search('confout.gro',out)

        if look:
            check=0
        else:
            check=1
            time.sleep(10)
            print "....."
    
    
    print "##########"
    print os.getcwd()
    print "##########"
    extract_protein="vmd -dispdev text confout.gro -e process_protein.tcl"
    commands.getoutput(extract_protein)
    
    commands.getoutput("*#*")
    print "Equilibration done\n"
    

def process_protein(my_working_dir,my_protein_path):
    print "Processing Protein........"
#    print "This will:"
#    print "1. Genearate topology"
#    print "2. Calculate charge on the protein"
#    print "3. Calculate number of water molecules needed to solvate"

    os.chdir(my_working_dir)
    commands.getoutput('rm -fr work_protein')
    os.mkdir('work_protein')
    #command_string='cp '+ my_protein_path + ' ' + 'work_protein/ptn.pdb'
    command_string='cp '+ my_working_dir + '/eq/ptn2.pdb ' + 'work_protein/ptn.pdb'
    commands.getoutput(command_string)

    extract_protein="vmd -dispdev text work_protein/ptn.pdb -e process_protein.tcl"
    commands.getoutput(extract_protein)
    commands.getoutput("cp work_protein/ptn2.pdb work_protein/ptn.pdb")
    
    os.chdir('work_protein')

    generate_protein_topology='gmx pdb2gmx -f %s -ff amber99sb-ildn -water spce -vsite hydrogens' % ('ptn.pdb')
    commands.getoutput(generate_protein_topology)
    box = 'gmx editconf -f %s -o %s -d 1.4 -bt cubic' % ('conf.gro','box.gro')
    commands.getoutput(box)

    # Protein charge calculation
    commands.getoutput("rm all_charges")
    commands.getoutput("egrep \" q \" topol.top > all_charges")
    total_charge=0.0
    with open("all_charges","r") as all_lines:
        for line in all_lines:
                a=line
                b=a.split()
                total_charge = total_charge + float(b[-1])

    solvate = 'gmx solvate -cp %s -cs -o %s' % ('box.gro','solvated.gro')
    commands.getoutput(solvate)
    
    number_of_water_molecules=commands.getoutput("egrep SOL solvated.gro | wc -l")
    number_of_water_molecules=int(number_of_water_molecules)/3
    reduced_number_of_water_molecules=number_of_water_molecules-int(number_of_water_molecules/5)
    resolvate='gmx solvate -cp %s -cs -o %s -maxsol %s' % ('box.gro','solvated.gro',reduced_number_of_water_molecules)
    commands.getoutput(resolvate)
    combine_topologies()
    print "Processing Protein Done"
    
    return total_charge,reduced_number_of_water_molecules


def charge_calculation(my_working_dir,total_charge_protein,total_charge_ligand,reduced_number_of_water_molecules):
    total_system_charge=total_charge_protein+total_charge_ligand

    ions=int(reduced_number_of_water_molecules*0.00227)
    na_ions=ions
    cl_ions=ions

    if total_system_charge<0:
        na_ions=na_ions+abs(total_system_charge)
    elif total_system_charge>0:
        cl_ions=cl_ions+abs(total_system_charge)
    return na_ions,cl_ions


def ligand_charge(my_working_dir,my_ligand_path):
    print "Processing Ligand........"
    os.chdir(my_working_dir)
    commands.getoutput('rm -fr work_ligand')
    os.mkdir('work_ligand')
    command_string='cp '+ my_ligand_path + ' ' + 'work_ligand/lig.mol2'
    commands.getoutput(command_string)
    os.chdir('work_ligand')
    mol=Chem.MolFromMol2File('lig.mol2')
    ligand_charge=Chem.GetFormalCharge(mol)
    build_topology='acpype -i lig.mol2 -n %s' % (ligand_charge)
    commands.getoutput(build_topology)
    print "Processing Ligand Done"
    return ligand_charge


def make_multiple_copies(my_working_dir,reduced_number_of_water_molecules,na_ions,cl_ions):
    print "Making multiple replicas................"
    os.chdir(my_working_dir)
    
    commands.getoutput("rm -fr system")
    os.mkdir('system')
    
    os.chdir('system')
    
    na_ions=int(na_ions)
    cl_ions=int(cl_ions)
    
    dir_initial='sys'
    
    for i in range (1,151):
        commands.getoutput("rm *#*")
        print "Processing replica no",i
        random_number=randint(0,100)
        random_number=str(random_number)
        dir_name=dir_initial+str(i)
        os.mkdir(dir_name)

        copy_topology='cp ../work_protein/topol.top '+dir_name+'/topol.top'
        copy_topology2='cp ../work_ligand/*.itp '+dir_name+'/'
        commands.getoutput(copy_topology)
        commands.getoutput(copy_topology2)
        topology=dir_name+'/topol.top'
        
        inserted_file_name=dir_name+'/'+str(i)+'.gro'
        solvated_file_name=dir_name+'/'+'solv'+str(i)+'.gro'
        tpr_file_name=dir_name+'/'+str(i)+'.tpr'
        ion_file_name=dir_name+'/'+'ions_'+str(i)+'.gro'

        gro_file_name=dir_name+'/ions.gro'
        tpr_file_name=dir_name+'/ions.tpr'
        tpr_file_name2=dir_name+'/'+str(i)+'.tpr'
        
        put_ligand_in_ptn='gmx insert-molecules -f %s -ci %s -nmol 1 -seed %s -o %s' % ('../work_protein/box.gro','../work_ligand/MOL_GMX.gro',random_number,inserted_file_name)
        commands.getoutput(put_ligand_in_ptn)
        
        solvate='gmx solvate -cp %s -cs -maxsol %s -o %s -p %s' % (inserted_file_name,reduced_number_of_water_molecules,solvated_file_name,topology)
        commands.getoutput(solvate)
        
        tpr_gen1='gmx grompp -f ../em.mdp -c %s -o %s -p %s' % (solvated_file_name,tpr_file_name,topology)
        commands.getoutput(tpr_gen1)
        
        ionize='echo SOL| gmx genion -s %s -np %s -nn %s -o %s -p %s'%(tpr_file_name,na_ions,cl_ions,ion_file_name,topology)
        commands.getoutput(ionize)
        
        tpr_gen2='gmx grompp -f ../em.mdp -c %s -o %s -p %s' % (ion_file_name,tpr_file_name,topology)
        commands.getoutput(tpr_gen2)
        
        prevdir=os.getcwd()
        os.chdir(dir_name)
        commands.getoutput("gmx mdrun -deffnm ions")
        commands.getoutput("rm -fr *#*")
        os.chdir(prevdir)
        
        
        tpr_gen3='gmx grompp -f ../vsite_md_gpu1.mdp -c %s -o %s -p %s -maxwarn 2' % (gro_file_name,tpr_file_name2,topology)
        commands.getoutput(tpr_gen3)
        
        run_tpr=os.getcwd()+'/'+tpr_file_name2
        d_path=os.getcwd()+'/'+dir_name

        delete_gro_files='rm '+dir_name+'/*.gro'
        commands.getoutput(delete_gro_files)
        
        
        run_on_gpu="job_submit_command -s %s" % (run_tpr,d_path)
        print "## 1 ##"
        print run_on_gpu
        commands.getoutput(run_on_gpu)
        
        check_gro_files='ls '+dir_name+'/*.gro'
        out=commands.getoutput(check_gro_files)

        print "Setting up the run for %s" % (dir_name) 

        check=1
        while check>0:
            out=commands.getoutput(check_gro_files)
            look=re.search('confout.gro',out)
            if look:
                check=0
            else:
                check=1
                print "....."
                time.sleep(30)

                
        for_production_input_gro=os.getcwd()+'/'+dir_name+'/confout.gro'
        for_production_input_top=os.getcwd()+'/'+dir_name+'/topol.top'
        for_production_output_tpr=os.getcwd()+'/'+dir_name+'/p'+str(i)+'.tpr'
        for_run_path=os.getcwd()+'/'+dir_name
        
        #print for_production_input_gro
        #print for_production_output_tpr
        #print for_run_path
        
        production_tpr="gmx grompp -f ../production.mdp -c %s -p %s -o %s -maxwarn 2" % (for_production_input_gro,for_production_input_top,for_production_output_tpr)
        #print production_tpr
        commands.getoutput(production_tpr)
        
        submit_production_tpr="job_submit_command -s %s" % (for_production_output_tpr)
        print submit_production_tpr
        commands.getoutput(submit_production_tpr)

def submit_jobs_to_gpu(my_working_dir,how_many):
    
    os.chdir(my_working_dir)
    print os.getcwd()
    
    for i in range (1,how_many+1):
        subdir_name="system/sys"+str(i)
        print i
        print subdir_name
        
        os.chdir(subdir_name)
        print os.getcwd()
        
        submit_production_tpr="job_submit_command -s p1.tpr" % (os.getcwd())
        print submit_production_tpr
        
        commands.getoutput(submit_production_tpr)
        
        os.chdir(my_working_dir)
        print os.getcwd()

def post_process():
	os.chdir(my_working_dir)
	os.mkdir("final_trajs")
	os.getoutput("cp ../scripts/post_process.pl")
	commands.getoutput("perl post_process.pl")

#============================ Main =============================================================

equilibrate(my_working_dir,my_protein_path)
total_charge_ligand=ligand_charge(my_working_dir,my_ligand_path)
total_charge_protein,reduced_number_of_water_molecules=process_protein(my_working_dir,my_protein_path)
na_ions,cl_ions=charge_calculation(my_working_dir,total_charge_protein,total_charge_ligand,reduced_number_of_water_molecules)
make_multiple_copies(my_working_dir,reduced_number_of_water_molecules,na_ions,cl_ions)
clean_up(my_working_dir)

