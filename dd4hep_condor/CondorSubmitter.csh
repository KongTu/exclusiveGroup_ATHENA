#! /bin/env csh

set maindir = /gpfs02/eic/ztu/Analysis/BeAGLE/eAu_diffractive_VM/18x110_Q2_1_20
cd $maindir
set Exec = `pwd`/myexec.csh

####### Initialize condor file
echo ""  > CondorFile
echo "Universe     = vanilla" >> CondorFile
echo "Executable   = ${Exec}" >> CondorFile
echo "getenv = true" >>CondorFile

set j = 1
while ( $j <= 800 )
	cp $maindir/S3VJn003 $maindir/py_${j}
	sed -i "1d" py_${j}
	sed -i "1i\'outForPythiaMode/input_temp_${j}.txt' ! output file name" py_${j}
	cp $maindir/inputFiles/eAu_18x110_Q2_1_20_0.01_0.95_40k_Shd3_ShdFac=1.32_AllVec9193_Jpsinodecay_tune3_newGID.inp $maindir/inputFiles/input_temp_${j}.inp
    sed -i "43d" $maindir/inputFiles/input_temp_${j}.inp
	sed -i "43i\PY-INPUT                                                              py_${j}" $maindir/inputFiles/input_temp_${j}.inp
	@ j++
end

set v = `cat beam_spread_electron.txt`
set v1 = `cat beam_spread_gold.txt`
@ i = 1
while ( $i <= $#v )
	sed -i "26d" $maindir/inputFiles/input_temp_${i}.inp
	sed -i "26i\MOMENTUM  $v[$i]    $v1[$i]"  $maindir/inputFiles/input_temp_${i}.inp
	@ i = $i + 1
end

cd $maindir
# split into chunks
set base = $maindir/inputFiles/input_temp_*.inp

foreach input ( ${base}* )
    # arguments
    set ArgName = ${input}
    set Files      = `basename ${input} | sed 's/.inp//g'`
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
condor_submit CondorFile

