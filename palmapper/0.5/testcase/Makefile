dataDir = data
resultDir ?= $(dataDir)
ext = 
binDir = ..
#nrThreads ?= $(shell grep processor /proc/cpuinfo | wc -l)
nrThreads=1

Unmapped = split_1m.000.unmapped$(ext)
MappedSam = split_1m.000.sam$(ext)
MappedBed = split_1m.000.bedx$(ext)
Coverage = split_1m.000.wiggle$(ext)
#Spliced = split_1m.000.spliced$(ext)
LogFile = split_1m.000.log$(ext)

#Output file names for each testcase
ResultFilesTest = $(MappedSam) $(Unmapped) $(Coverage)
ResultFilesUnspliced = $(MappedSam) $(Unmapped) $(Coverage)
ResultFilesPolytrim = $(MappedSam) $(Unmapped) $(Coverage)
ResultFilesRtrim =  $(MappedSam) $(Unmapped) $(Coverage)
ResultFilesNoss =  $(MappedSam) $(Unmapped) $(Coverage)
ResultFilesBed =  $(MappedBed) $(Unmapped) $(Coverage)

#Output directory for each testcase
resultDirTest = $(resultDir)/map
resultDirUnspliced = $(resultDir)/mapUnspliced
resultDirPolytrim = $(resultDir)/mapPolytrim
resultDirRtrim = $(resultDir)/mapRtrim
resultDirNoss = $(resultDir)/mapNoSpliceSitePred
resultDirBed = $(resultDir)/mapBed

#Output subdirectory according to index
resultSubDirArray=array-index
resultSubDirBwa=bwt-index


#Checksum file for each testcase
result = result.md5
resultUnspliced =resultUnspliced.md5
resultPolytrim = resultPolytrim.md5
resultRtrim = resultRtrim.md5
resultNoss = resultNoss.md5
resultBed = resultBed.md5

DnaFile = $(dataDir)/c_elegans.WS209.dna.fa
SpliceAccFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/acc_pred.bspf
SpliceDonFile = $(dataDir)/C_elegans_SpliceSitePred_WS209/don_pred.bspf
# Since make is unable to handle tools creating multiple products, we use only one representative
SpliceFiles = $(SpliceAccFile) # $(SpliceDonFile)
# Since make is unable to handle tools creating multiple products, we use only one representative
IndexFiles := $(DnaFile).cid # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta
BwaIndexFiles := $(DnaFile).bwt # $(DnaFile).mfd $(DnaFile)/mrc $(DnaFile).mta

CurrentDir := $(shell pwd)

#PALMapper parameters
inFiles = -i $(dataDir)/c_elegans.WS209.dna.fa -q $(dataDir)/split_1m.000
accDon = -acc $(SpliceAccFile)/contig_%i%c -don $(SpliceDonFile)/contig_%i%c 
outSam= -f sam 
outBed= -f bedx
qpalma =  -L 35 -K 12 -C 55 -S -qpalma $(dataDir)/parameters.qpalma 
other = -M 4 -G 1 -E 4 -z 1 -l 20 -seed-hit-cancel-threshold 1000
polytrim = -polytrim 40
rtrim = -rtrim 25
report =
threads = -threads $(nrThreads)
bwa = -bwa 12

#Parameter base for all testcases
allFlags=$(inFiles) $(other) $(threads) $(report) -rlim 10000
MapCmd = $(binDir)/palmapper $(allFlags) 

#outFiles =  -u $(resultDir)/$(Unmapped) #-H $(resultDir)/$(Spliced)
#-o $(resultDir)/$(MappedSam)
#outCoverage =  -report-coverage-wig $(resultDir)/$(Coverage)

numtest=12
numtestarray=6
numtestbwa=6

all: mapall1 bwamapall1
	@ (cd $(CurrentDir); bash  -c "./test_output.sh $(numtest)")

test: map

mapall1: map mapUnspliced mapPolytrim mapRtrim mapNoss mapBedx

bwamapall1: bwamap bwamapUnspliced bwamapPolytrim bwamapRtrim bwamapNoss bwamapBedx

mapall: mapall1
	@ (cd $(CurrentDir); bash  -c "./test_output.sh $(numtestarray)")

bwamapall: bwamapall1
	@ (cd $(CurrentDir); bash  -c "./test_output.sh $(numtestbwa)")

$(dataDir)/split_1m.000: 
	bzip2 -d $(dataDir)/split_1m.000.bz2

$(DnaFile): $(SpliceFiles)

$(IndexFiles): $(DnaFile) # pmindex
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/pmindex -i $(DnaFile) -v"

$(BwaIndexFiles): $(DnaFile) # bwa
	mkdir -p $(dataDir)
	bash -c "time $(binDir)/bwa index $(DnaFile) "

$(SpliceFiles):
	mkdir -p $(dataDir)
	wget -c ftp://ftp.tuebingen.mpg.de/pub/fml/raetsch-lab/predictions/splice/C_elegans/C_elegans_SpliceSitePred_WS209.tar.bz2
	(cd $(dataDir) && tar xvf ../C_elegans_SpliceSitePred_WS209.tar.bz2 && rm ../C_elegans_SpliceSitePred_WS209.tar.bz2) 
	bzip2 -d $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta.bz2
	mv $(dataDir)/C_elegans_SpliceSitePred_WS209/genome/C_elegans_WS209_genome.fasta $(DnaFile)

valgrind: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "valgrind $(MapCmd)"

cachegrind: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
	mkdir -p $(dataDir)
	bash -c "valgrind --tool=cachegrind $(MapCmd)"

clean:
	rm -rf $(dataDir)/C_elegans_SpliceSitePred_WS209/ $(dataDir)/c_elegans.WS209.dna.fa* $(resultDirTest) $(resultDirUnspliced) $(resultDirRtrim) $(resultDirPolytrim) $(resultDirNoss) $(resultDirBed) $(CurrentDir)/*.test_passed $(CurrentDir)/*.test_failed

allmd5: result.md5 resultUnspliced.md5 resultPolytrim.md5 resultRtrim.md5 resultNoss.md5 resultBed.md5

#Checksum file creation rules
result.md5: $(resultDirTest)/$(resultSubDirArray)/$(MappedSam) $(resultDirTest)/$(resultSubDirArray)/$(Unmapped) $(resultDirTest)/$(resultSubDirArray)/$(Coverage)
	mkdir -p $(resultDirTest)/$(resultSubDirArray)
	(cd $(resultDirTest)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesTest))) > $(result)

resultUnspliced.md5: $(resultDirUnspliced)/$(resultSubDirArray)/$(MappedSam) $(resultDirUnspliced)/$(resultSubDirArray)/$(Unmapped) $(resultDirUnspliced)/$(resultSubDirArray)/$(Coverage)
	echo $(resultDirUnspliced)/$(resultSubDirArray)
	mkdir -p $(resultDirUnspliced)/$(resultSubDirArray)
	(cd $(resultDirUnspliced)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesUnspliced))) > $(resultUnspliced)

resultRtrim.md5: $(resultDirRtrim)/$(resultSubDirArray)/$(MappedSam) $(resultDirRtrim)/$(resultSubDirArray)/$(Unmapped) $(resultDirRtrim)/$(resultSubDirArray)/$(Coverage)
	mkdir -p $(resultDirRtrim)/$(resultSubDirArray)
	(cd $(resultDirRtrim)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesRtrim))) > $(resultRtrim) 

resultPolytrim.md5: $(resultDirPolytrim)/$(resultSubDirArray)/$(MappedSam) $(resultDirPolytrim)/$(resultSubDirArray)/$(Unmapped) $(resultDirPolytrim)/$(resultSubDirArray)/$(Coverage)
	mkdir -p $(resultDirPolytrim)/$(resultSubDirArray)
	(cd $(resultDirPolytrim)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesPolytrim))) >  $(resultPolytrim) 

resultNoss.md5: $(resultDirNoss)/$(resultSubDirArray)/$(MappedSam) $(resultDirNoss)/$(resultSubDirArray)/$(Unmapped) $(resultDirNoss)/$(resultSubDirArray)/$(Coverage)
	mkdir -p $(resultDirNoss)/$(resultSubDirArray)
	(cd $(resultDirNoss)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesNoss))) >  $(resultNoss) 

resultBed.md5: $(resultDirBed)/$(resultSubDirArray)/$(MappedBed) $(resultDirBed)/$(resultSubDirArray)/$(Unmapped) $(resultDirBed)/$(resultSubDirArray)/$(Coverage)
	mkdir -p $(resultDirBed)/$(resultSubDirArray)
	(cd $(resultDirBed)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesBed))) >  $(resultBed) 


#Testcases with array-based index
map: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000	
#@ echo "Array-based index: Basic Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirTest)/$(resultSubDirArray)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) -u $(resultDirTest)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirTest)/$(resultSubDirArray)/$(MappedSam) -report-coverage-wig $(resultDirTest)/$(resultSubDirArray)/$(Coverage)> $(resultDirTest)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirTest)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesTest)) > /dev/null)
	@ (cd $(resultDirTest)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(result) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed


mapUnspliced: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: Unspliced Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirUnspliced)/$(resultSubDirArray)
	bash -c "time $(MapCmd) $(outSam) -u $(resultDirUnspliced)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirUnspliced)/$(resultSubDirArray)/$(MappedSam) -report-coverage-wig $(resultDirUnspliced)/$(resultSubDirArray)/$(Coverage) > $(resultDirUnspliced)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirUnspliced)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesUnspliced))  > /dev/null)
	@ (cd $(resultDirUnspliced)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(resultUnspliced) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

mapRtrim: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: Rtrim Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirRtrim)/$(resultSubDirArray)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim) -u $(resultDirRtrim)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirRtrim)/$(resultSubDirArray)/$(MappedSam) -report-coverage-wig $(resultDirRtrim)/$(resultSubDirArray)/$(Coverage) > $(resultDirRtrim)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirRtrim)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesRtrim))  > /dev/null)
	@ (cd $(resultDirRtrim)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(resultRtrim) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

mapPolytrim: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: Polytrim Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirPolytrim)/$(resultSubDirArray)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim) -u $(resultDirPolytrim)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirPolytrim)/$(resultSubDirArray)/$(MappedSam) -report-coverage-wig $(resultDirPolytrim)/$(resultSubDirArray)/$(Coverage) > $(resultDirPolytrim)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirPolytrim)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesPolytrim))  > /dev/null)
	@ (cd $(resultDirPolytrim)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(resultPolytrim) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

mapNoss: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: No splice predictions Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirNoss)/$(resultSubDirArray)
	bash -c "time $(MapCmd) -no-ss-pred $(outSam) $(qpalma) -u $(resultDirNoss)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirNoss)/$(resultSubDirArray)/$(MappedSam) -report-coverage-wig $(resultDirNoss)/$(resultSubDirArray)/$(Coverage) > $(resultDirNoss)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirNoss)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesNoss))  > /dev/null)
	@ (cd $(resultDirNoss)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(resultNoss) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

mapBedx: $(binDir)/palmapper $(DnaFile) $(IndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: Bedx Testcase"	
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirBed)/$(resultSubDirArray)
	bash -c "time $(MapCmd) $(accDon)  $(outBed) $(qpalma) -u $(resultDirBed)/$(resultSubDirArray)/$(Unmapped) -o $(resultDirBed)/$(resultSubDirArray)/$(MappedBed) -report-coverage-wig $(resultDirBed)/$(resultSubDirArray)/$(Coverage) > $(resultDirBed)/$(resultSubDirArray)/$(LogFile)"
	@ (cd $(resultDirBed)/$(resultSubDirArray) ; md5sum $(notdir $(ResultFilesBed))  > /dev/null)
	@ (cd $(resultDirBed)/$(resultSubDirArray) ; md5sum -c $(CurrentDir)/$(resultBed) && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

#Testcases with bwa-based index
bwamap: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Bwt-based index: Basic Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirTest)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) -u $(resultDirTest)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirTest)/$(resultSubDirBwa)/$(MappedSam) -report-coverage-wig $(resultDirTest)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirTest)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirTest)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesTest))  > /dev/null)
	@ (cd $(resultDirTest)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(result)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

bwamapUnspliced: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Bwt-based index: Unspliced Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirUnspliced)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) $(outSam) -u $(resultDirUnspliced)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirUnspliced)/$(resultSubDirBwa)/$(MappedSam) -report-coverage-wig $(resultDirUnspliced)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirUnspliced)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirUnspliced)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesUnspliced))  > /dev/null)
	@ (cd $(resultDirUnspliced)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(resultUnspliced)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

bwamapRtrim: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Bwt-based index: Rtrim Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirRtrim)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(rtrim) -u $(resultDirRtrim)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirRtrim)/$(resultSubDirBwa)/$(MappedSam) -report-coverage-wig $(resultDirRtrim)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirRtrim)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirRtrim)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesRtrim))  > /dev/null)
	@ (cd $(resultDirRtrim)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(resultRtrim)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

bwamapPolytrim: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Bwt-based index: Polytrim Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirPolytrim)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) $(accDon)  $(outSam) $(qpalma) $(polytrim) -u $(resultDirPolytrim)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirPolytrim)/$(resultSubDirBwa)/$(MappedSam) -report-coverage-wig $(resultDirPolytrim)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirPolytrim)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirPolytrim)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesPolytrim))  > /dev/null)
	@ (cd $(resultDirPolytrim)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(resultPolytrim)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

bwamapNoss: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: No splice predictions Testcase"
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirNoss)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) -no-ss-pred $(outSam) $(qpalma) -u $(resultDirNoss)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirNoss)/$(resultSubDirBwa)/$(MappedSam) -report-coverage-wig $(resultDirNoss)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirNoss)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirNoss)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesNoss))  > /dev/null)
	@ (cd $(resultDirNoss)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(resultNoss)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed

bwamapBedx: $(binDir)/palmapper $(DnaFile) $(BwaIndexFiles) $(dataDir)/split_1m.000
#@ echo "Array-based index: Bedx Testcase"	
#@ echo "-------------------------------------------------"
	@ rm -f $(CurrentDir)/$@.test_passed $(CurrentDir)/$@.test_failed
	@ mkdir -p $(resultDirBed)/$(resultSubDirBwa)
	bash -c "time $(MapCmd) $(accDon)  $(outBed) $(qpalma) -u $(resultDirBed)/$(resultSubDirBwa)/$(Unmapped) -o $(resultDirBed)/$(resultSubDirBwa)/$(MappedBed) -report-coverage-wig $(resultDirBed)/$(resultSubDirBwa)/$(Coverage) $(bwa) > $(resultDirBed)/$(resultSubDirBwa)/$(LogFile)"
	@ (cd $(resultDirBed)/$(resultSubDirBwa) ; md5sum $(notdir $(ResultFilesBed))  > /dev/null)
	@ (cd $(resultDirBed)/$(resultSubDirBwa) ; md5sum -c $(CurrentDir)/$(resultBed)  && touch $(CurrentDir)/$@.test_passed) || touch $(CurrentDir)/$@.test_failed


