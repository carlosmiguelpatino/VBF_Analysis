#Script for running analysis code through signal files
make compile_ROOT_Delphes

signalPath='/home/aflorez/HN_VBF_signal'
outputPath='/home/aflorez/HN_Analysis/TauChannel'


for i in `seq 0 10 90`;
do
  ./PhenoAnalyzer ${signalPath}/VBF_m_n3_1p5_tau_decay_vbfcuts_QCD0_seed${i}/Events/run_01/output_delphes.root ${outputPath}/signal_1p5TeV/output_signal_seed${i}.root
done

hadd -f ${outputPath}/signal_1p5TeV/signal_1p5TeV_700mass_30MET_50pt1.root ${outputPath}/signal_1p5TeV/output_signal_seed*.root
rm ${outputPath}/signal_1p5TeV/output_signal_seed*.root

./PhenoAnalyzer ${signalPath}/VBF_m_n3_1p0_tau_decay_vbfcuts/Events/run_01/output_delphes.root ${outputPath}/signal_1p0TeV/signal_1p0TeV_700mass_30MET_50pt1.root 

./PhenoAnalyzer ${signalPath}/VBF_m_n3_1p25_tau_decay_vbfcuts/Events/run_01/output_delphes.root ${outputPath}/signal_1p25TeV/signal_1p25TeV_700mass_30MET_50pt1.root

./PhenoAnalyzer ${signalPath}/VBF_m_n3_2p0_tau_decay_vbfcuts/Events/run_01/output_delphes.root ${outputPath}/signal_2p0TeV/signal_2p0TeV_700mass_30MET_50pt1.root
