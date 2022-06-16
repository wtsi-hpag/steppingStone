# steppingStone v1.0
A pipeline to identify chromothripsis breakpoints and trace cancer rearrangements

### Download and Compile:

    $ git clone  https://github.com/wtsi-hpag/steppingStone.git 
    $ cd steppingStone 
    $ bash install.sh
		
If everything compiled successfully you must see the final comment: 
		"Congrats: installation successful!"		

#### External packages
The genome aligner BWA (https://github.com/lh3/bwa) and minimap2 (https://github.com/lh3/minimap2) are downloaded and compiled by steppingStone.

### Run the pipelines

           $ /full/path/to/steppingStone/src/steppingStone -nodes <nodes>  \
                 -reads /lustre/team117/zn1/project/cancer/PacBio-HiFi-sample.fastq.gz \
                 reference.fasta breakpoints.dat  

       Note: "-reads " expects a full path or a file in the working directory.

