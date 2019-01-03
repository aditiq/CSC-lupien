
#!/bin/bash

## Run KCCA with different cluster numbers

#for f in {3..100} ;
#do
    f=100  ; qsub -cwd -N KCCA.Enhancer."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.Enhancer.sh $f ;
#done

#for f in {3..100} ;
#do
    f=100  ; qsub -cwd -N KCCA.Promoter."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.Promoter.sh $f ;
#done

f=100
qsub -cwd -N KCCA.LSCnegvspos."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/LSCposandneg.Consensus.Catalogue.Binarymat.txt" "LSCnegvspos" $f ;
qsub -cwd -N KCCA.GBMvsHF."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/GBM.HF.Consensus.Catalogue.Binarymat.txt" "GBMvsHF" $f ;
qsub -cwd -N KCCA.PFA."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/PFA.Consensus.Catalogue.Binarymat.txt" "PFA" $f ;
qsub -cwd -N KCCA.GBM."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/GBM.Consensus.Catalogue.Binarymat.txt" "GBM" $f ; ##Did not finish running.Couldnt converge
qsub -cwd -N KCCA.LSCp."$f" -q all.q -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/LSCp.Consensus.Catalogue.Binarymat.txt" "LSCp" $f ;
qsub -cwd -N KCCA.CSC.ESC."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.other.sh "data/ConsensusSet/PCSC1/CSC.ESC.Consensus.Catalogue.Binarymat.txt" "CSC.ESC" $f ;


num="5 25 50 75 100"
for f in $num ;
do
    qsub -cwd -N KCCA.KitchenSink.Promoter."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.kitchensink.promoter.sh $f ;
done


num="5 25 50 75 100"
for f in $num ;
do
    qsub -cwd -N KCCA.KitchenSink.Enhancer."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.kitchensink.enhancer.sh $f ;
done

f=100
qsub -cwd -N KCCA.KitchenSink2.Enhancer."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.kitchensink2.enhancer.sh $f ;
qsub -cwd -N KCCA.KitchenSink2.Promoter."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06_KCCA.kitchensink2.promoter.sh $f ;


num="5 25 50 75 100"
for f in $num ;
do
    qsub -cwd -N KCCA.Brain."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06b.2_KCCA.Brain.sh $f ;
    qsub -cwd -N KCCA.Blood."$f" -q lupiengroup -e stdout/kcca -o stdout/kcca scripts/run_06b.2_KCCA.Blood.sh $f ;

done
