Delphes Generation
==================

1. Copy Delphes paramaters card
   * now available in HN_VBF_signal/VBF_m_n3_1TeV/Events/run_01/delphes_card_CMS.tcl
2. Unzip .hep file
   * gunzip compressed_file.hep.gz
3. DelphesSTDHEP delphes_card.tlc output_file.root input_file.hep
   * DelphesSTDHEP delphes_card_CMS.tcl output_delphes.root tag_1_pythia_events.hep
