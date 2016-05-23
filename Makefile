all: partitionByPhasedSNVs partitionBySNV readToSNVList testSeqan testReadSpeed 

SEQAN=/net/eichler/vol5/home/mchaisso/software/seqan/include
BOOSTLIB=/net/eichler/vol5/home/mchaisso/software/lib
BOOST=/net/eichler/vol5/home/mchaisso/software/include
BLASR=/net/eichler/vol5/home/mchaisso/projects/blasr/cpp/common
VCFLIB=/net/eichler/vol5/home/mchaisso/software/vcflib/vcflib
HTSLIB=$(VCFLIB)/tabixpp/htslib
CPPOPTS= -O3 -g
CPP=g++ 

partitionByPhasedSNVs: PartitionByPhasedSNVs.cpp FastaIndex.h
	$(CPP) -static $(CPPOPTS) $< \
     -I $(SEQAN) \
     -I $(BLASR) \
     -I $(BOOST) \
     -I $(VCFLIB)/include -I $(HTSLIB) \
     -L $(BOOSTLIB) -l boost_program_options \
     -L $(VCFLIB)/lib -l vcflib  -L $(HTSLIB) -l hts  -lpthread -lz \
     -o $@ 

