#! /bin/env csh

set maindir = /gpfs02/eic/ztu/ATHENA/exclusiveGroup_ATHENA/dd4hep_condor/physics_benchmarks
cd $maindir
set Exec = `PWD`/benchmarks/diffractive_phi/gen.sh

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${Exec}" >> CondorFile
echo "getenv = true" >>CondorFile

cd $maindir
# split into chunks
for j in {1..3}
    # arguments
    set ArgName = --ebeam 18 --pbeam 110 --config test_${j}
    set Files      = `basename ${input} | sed 's/.hepmc//g'`
    echo "input name " ${Files}    
    # Logfiles.
    set LogFile    = $maindir/logs/${Files}.out
    set ErrFile    = $maindir/logs/${Files}.err

    ### hand to condor
    set Args = (${ArgName})
    echo "" >> CondorFile
    echo "Output       = ${LogFile}" >> CondorFile
    echo "Error        = ${ErrFile}" >> CondorFile
    echo "Arguments    = ${Args}" >> CondorFile
    echo "Queue" >> CondorFile   

    echo Submitting:
    echo $Exec $Args
    echo "Logging output to " $LogFile
    echo "Logging errors to " $ErrFile
    echo

end
# condor_submit CondorFile

