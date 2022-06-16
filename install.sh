#!/bin/bash


projdir=`pwd`

bindir=$projdir/src/step-bin/
mkdir -p $bindir
mkdir -p $projdir/src/log/

errs=0

##### Download and install BWA ######

echo "Downloading and installing BWA"
if [[ ! -s $bindir/bwa ]]; then

    if [[ ! -d $projdir/src/bwa ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/bwa.git &> $projdir/src/log/bwa_cloning.log
    fi

    if [[ ! -s $projdir/src/bwa/bwa ]]; then
	cd $projdir/src/bwa
	make &> $projdir/src/log/bwa_installation.log
    fi

    cp bwa $bindir
fi

if  [[ ! -s $bindir/bwa ]]; then
    echo " !! Error: bwa not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if bwa was downloaded properly:" $projdir/src/log/bwa_cloning.log 
    echo "   Check if the bwa was compiled properly:" $projdir/src/log/bwa_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/bwa/bwa $bindir/bwa 
    
    errs=$(($errs+1))
else
    echo " BWA succesfully installed!"
    rm -rf $projdir/src/bwa/
fi

##### Download and install minimap2 ######

echo "Downloading and installing minimap2"
if [[ ! -s $bindir/minimap2 ]]; then

    if [[ ! -d $projdir/src/minimap2 ]]; then
	cd $projdir/src/
	git clone https://github.com/lh3/minimap2 &> $projdir/src/log/minimap2_cloning.log
    fi

    if [[ ! -s $projdir/src/minimap2/minimap2 ]]; then
	cd $projdir/src/minimap2
	make &> $projdir/src/log/minimap2_installation.log
    fi

    cp minimap2 $bindir
fi

if  [[ ! -s $bindir/minimap2 ]]; then
    echo " !! Error: minimap2 not installed properly!"; 
    echo "   Please check the log files:" 
    echo "   Check if minimap2 was downloaded properly:" $projdir/src/log/minimap2_cloning.log 
    echo "   Check if the minimap2 was compiled properly:" $projdir/src/log/minimap2_installation.log

    # Cleaning up
    cd $projdir/src
    rm -rf $projdir/src/minimap2/minimap2 $bindir/minimap2 
    
    errs=$(($errs+1))
else
    echo " minimap2 succesfully installed!"
    rm -rf $projdir/src/minimap2/
fi


###### Compile steppingStone sources ######

echo; echo "Compiling steppingStone sources"

srcs=( step_fastq step_linkStones step_breakSort step_number step_breakProcess step_cleanProcess stepBreakPoint step_processStones step_shortReads step_edgeStones step_sortStones steppingStone)

cd $projdir/src
make &> $projdir/src/log/sources_compilation.log

echo; echo "Checking installation:"
for src in "${srcs[@]}"; do
    if [[ ! -s $bindir/$src ]]; then 
        echo " !! Error: executable $src missing in $bindir"
	echo "    Please check for errors the log file:" $projdir/src/log/sources_*	
        errs=$(($errs+1))
    fi
done

if [  $errs -gt 0 ]; then echo; echo " ****  Errors occurred! **** "; echo; exit; 
else echo " Congrats: installation successful!"; fi




