qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.pfa.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.PFA.txt ## Step2
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.lscp.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.LSCp.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.gbm.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.GBM.txt ## Step2
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.pcsc1.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.PCSC1.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.hf.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.HF.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.brain.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.BRAIN.txt ## 
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.hESCs.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.hESC.txt ## Running 
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.lscn.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.LSCn.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.lscbulk.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.LSCbulk.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.HSC.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.HSC.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.hematdiff.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.HEMATDIFF.txt ## Running
qsub -cwd -q lupiengroup -e stdout/repeat/ -o stdout/repeat/ -N step1.hematprog.pipeline scripts/repeat/pipeline/Step1.bg.gen.part2.sh scripts/repeat/pipeline/config.HEMATPROG.txt ## Running
