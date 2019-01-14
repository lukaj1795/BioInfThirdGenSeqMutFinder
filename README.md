# BioInfThirdGenSeqMutFinder
Third generation sequencing - mutation finder

project for [Bioinformatics course](https://www.fer.unizg.hr/predmet/bio)

[Faculty of Electrical Engineering and Computing, Zagreb, Croatia](https://www.fer.unizg.hr/en)

[Luka Jukić](https://github.com/lukaj1795)

[Ivan Moštak](https://github.com/IvanMostak)

[Data used for testing](https://www.dropbox.com/s/xxj53t44rehlc77/Bioinfo_18_19_train_data.tar.gz?dl=0)

## Installation and running
```
Download or copy from git, then go to source folder

# compile
g++ *.cpp -O3 -m64 -o main.exe -std=c++11

# running
./main.exe ecoli.fasta ecoli_simulated_reads.fasta

# running with cutom w and k
./main.exe ecoli.fasta ecoli_simulated_reads.fasta 7 14

```
