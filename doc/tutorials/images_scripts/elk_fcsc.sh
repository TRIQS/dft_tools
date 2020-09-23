#!/bin/bash

#(Crude) shell script for FCSC calculations using the Elk-TRIQS setup.

#You may put job scheduler (PBS, Slurm or other) information here
#Note mpirun command lines may have to be altered for the 
#scheduler used


#
#the Elk executable (may have to update path and version used)
export elk_ex="${HOME}/elk-6.2.8/src/elk"
#this variable is already linked in the $PATH in .bashrc
export triqs_ex="python"

#initial variables which can be tweeked
#max loop index - Input variable
export imax=20
#python file prefix - user's script
export triqs_name="dft_dmft_cthyb_elk"


#DFT+DMFT self-consistent cycle shell output file
export out_file="FCSC.OUT"
#elk interface tasks
export e2t_task=805
export t2e_task=808
#Master input file for the dft interface calculation.
export elk_master_file="elk_master.in"
#input elk file
export elk_infile="elk.in"
#output elk file
export elk_outfile="elk.out"
#output triqs file
export triqs_outfile="SC_DMFT.OUT"
#this variable will be used to replace the string TASK in elk_master.in
export task="TASK"
#check if master file exists
if [ ! -e "$elk_master_file" ]; then
  echo "No master input file" > $out_file
  exit
fi

echo " "> $elk_outfile
echo " "> $triqs_outfile

echo "Maximum number of iterations = "$imax > $out_file
echo "Number of processors used = "$NUMPROC >> $out_file
for (( i=1; i<=$imax; i++)); do 
  echo "SC loop "$i" of "$imax >> $out_file
  echo " ">> $out_file

  #copy master input to elk.in
  cp $elk_master_file $elk_infile
  #replace TASK string to interface task
  sed -i "s/$task*/$e2t_task/" $elk_infile
  #executes elk (generate projectors)
  echo "SC loop "$i" of "$imax >> $elk_outfile
  mpirun $elk_ex >> $elk_outfile
  echo " ">> $elk_outfile
  echo "Interfaced Elk to TRIQS via task "$e2t_task >> $out_file

  #executes triqs
  echo "SC loop "$i" of "$imax >> $triqs_outfile
  mpirun $triqs_ex $triqs_name".py" >> $triqs_outfile
  echo " " >> $triqs_outfile
  echo "TRIQS DMFT loop(s) done " >> $out_file

  #copy master input to elk.in
  cp $elk_master_file $elk_infile
  #replace TASK string to interface task
  sed -i "s/$task*/$t2e_task/" $elk_infile 
  #executes elk (density update)
  mpirun $elk_ex >> $elk_outfile
  echo " ">> $elk_outfile
  echo "Interfaced TRIQS dmft outputs into elk via task "$t2e_task >> $out_file
  echo " ">> $out_file
done
#end


