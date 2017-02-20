#Script for running analysis code through signal files

LAST=$1 #Last index of samples
signalPath='/home/aflorez/HN_VBF_signal'
outputPath='/home/aflorez/HN_Analysis/TauChannel/signal'


for i in `seq 0 10 ${LAST}`;
do
  ./PhenoAnalyzer ${signalPath}/VBF_m_n3_1p5_tau_decay_vbfcuts_QCD0_seed${i}/Events/run_01/output_delphes.root ${outputPath}/output_signal_seed${i}.root
done

cd ${outputPath}

hadd output_signal.root *.root
