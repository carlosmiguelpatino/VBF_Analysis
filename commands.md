Generate Delphes File
==================

1. Copy Delphes paramaters card
   * now available in HN_VBF_signal/VBF_m_n3_1TeV/Events/run_01/delphes_card_CMS.tcl
2. Unzip .hep file
   * gunzip compressed_file.hep.gz
3. DelphesSTDHEP delphes_card.tlc output_file.root input_file.hep
   * DelphesSTDHEP delphes_card_CMS.tcl output_delphes.root tag_1_pythia_events.hep
   
Generate MadGraph5 Simulation
=======================

1. Open MadGraph5
   * mg5_mg5_aMC
2. Import model for simulation using `import`
   * import model SM_HeavyN_NLO
3. Generate process using `generate`
   * generate p p > n3 ta+ j j, n3 > ta+ j j QCD=0
4. If other processes neede, use `add process`
5. Set and ouput folder for all files generated during simulation using `output`
6. Launch simulation using `launch -m`

