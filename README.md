#Indriidae genome analysis

##1. Trimming the raw trees with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```bash
java -jar trimmomatic-0.36.jar PE\
	-phred33\  
	indrifilea1.fastq indrifilea2.fastq\  
	indrifilea1_paired.fastq indrifilea1_unpaired.fastq\  
	indrifilea2_paired.fastq indrifilea2_unpaired.fastq\  
	ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10\  
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36\  
```
##2. Extracting CDS from Genbank with the python script--[gbgenecdsextract.py](https://github.com/leirhyh/Indriidae-genome-analysis/blob/master/gbgenecdsextract.py) and select 1400 CDS.
```bash 
python gbgenecdsextract.py -input\_path sifakamRNA -output\_path geneselected
```

##3. Using [cGAP](https://github.com/TheCulliganMan/cgap) pipeline to extract the above CDS
```bash
docker pull theculliganman/cgap:latest  

docker run -itv indridir:/work theculliganman/cgap /bin/bash

python run\_cgap.py -refs\_path indri_ref\ 
	-forward indri1r1.fastq indri2tr1.fastq indri3r1.fastq\  	indri4tr1.fastq .... indri22tr1.fastq\ 
	-reverse indri1r2.fastq  indri2tr2.fastq  indri3r2.fastq\
	indri4tr2.fastq  .... indri22tr2.fastq\
	-c 8 -format\_db\
```

##4. Selecting genes based on 50% missing data by python scripts
* sorting command by ordered\_phy\_sortindriidae.py
```bash
python ordered_phy_sortindriidae.py \
	-input_path $STARTING_PHY_DIR \
	-output_path $SORTED_PHY_DIR
```
* convert phylip-relaxed file to standard phylip file by [phyluce_align_convert_one_align_to_another] (https://github.com/faircloth-lab/phyluce/blob/master/bin/align/phyluce_align_get_informative_sites)
```bash
python /home/genetics/anaconda2/pkgs/phyluce-1.5.0-py27_0/bin/phyluce_align_convert_one_align_to_another \
	--alignments $SORTED_PHY_DIR \
	--output $FORMAL_PHY_DIR \
	--input-format phylip-relaxed \
	--output-format phylip
```
* select phylip files based on missing data by [phyluce_align_get_informative_sites_reworked_ryan.py](https://github.com/TheCulliganMan/phyluce/blob/master/bin/align/phyluce_align_get_informative_sites_reworked_ryan.py)
```bash
 python phyluce\_align\_get\_informative\_sites\_reworked\_ryan \
	--input $FORMAL_PHY_DIR \
	--cores 60 \
	--output informative-file \
	--informative_sites_path informative-folder \
	--input-format phylip \
	--min_all_seqs_parsimony 50 \
	--output_phy_dir selectedphyfiles\
```

##5. Running [raxml](http://sco.h-its.org/exelixis/software.html) tree with 500 bootstrap replicates with python scripts.
```bash
python parallel\_process\_run\_rxml.py -input_path genefolder\
	-output genetreefolder -c 32\  
```

##6. Running species tree for all the genes with [ASTRAL](https://github.com/smirarab/ASTRAL)
```bash
 java -jar astral.4.8.o.jar -i combinedtree.tre\
	-a map.txt -o combinedtreeastra.tre\
	2> combinedtreeastral.log\  

 java -jar astral.4.0.11.jar -q combinedtreeastra.tre \  
	-i combinedtree.tre -a map.txt -o combinedtreeastrala.tre\
	2> combinedtreeastrala.log\

 java -jar astral.4.10.11.jar -q combinedtreeastra.tre\
	-i combinedtree.tre -a map.txt -t 4 
	-o combinedtreeastralat4.tre\
	2> combinedtreeastralat4.log\

 java -jar astral.4.10.11.jar -i combinedtree.tre\
	-a map.txt -b bootstraplist_files\
	-o combinedtreeatralbp\
	2> combinedtreeatralbp.log\
```

##7. Generating species trees for different randomly selected gene sets.
* Randomly selected gene sets by python scripts--[generandomselect.py](https://gist.github.com/leirhyh/932136d61ac9fd40d11e0ee79b151637)  
```bash	
python generandomselect.py -input\_path totalgenefile\
	-output\_path geneselectedfile -sample_num 200\
```
* As described in Step 6
* [MP-EST online version](http://bioinformatics.publichealth.uga.edu/SpeciesTreeAnalysis/index.php) 

##8. Comparing the above species trees by RF distance in [R3.2.3](https://www.r-project.org/) and select gene  set presenting the whole data set for following analyses

* setwd("C:\\Users\\Documents\\Genomic sequence sequencing")  

* library([ape](https://cran.r-project.org/web/packages/ape/index.html))  

* library([phangorn](http://cran.fhcrc.org/web/packages/phangorn/index.html))  
 
* tree1 <- read.tree( "basetreefile.tre")  

* tree2 <- read.tree( "filetocompare.tre")  

* RF.dist(tree1,tree2)  

##9. Data partition for the selected gene set by kmeans in [Partitionfiner 2.0](https://github.com/brettc/partitionfinder/releases)
```bash
docker pull theculliganman/partitionfinder
docker run -itv /partitiontest:/partitionfinder_data\ 
	theculliganman/partitionfinder /bin/bash\

python /partitionfinder/PartitionFinder.py\ 
	/partitionfinder_data/indripartition.phy
```
##10. Running [MrBayes v3.2.5](http://mrbayes.sourceforge.net/) with data partition from step 9.
```bash
mpirun -np 8 mb genefile.nex
```

##11. Running [BEASTv1.8.2](http://beast.bio.ed.ac.uk/downloads) with data partition from step 9 with the selected clock model.
```bash
java -jar beast.jar beagle bealge_cpu genefile.xml
```
