hadd -f hists_EWKWJetsWmn.root hists_EWKWJetsWmn_*.root  
hadd -f hists_EWKWJetsWen.root hists_EWKWJetsWen_*.root  
python ../../printYields.py -c EWKWJetsWmn -i hists_EWKWJetsWmn.root -s ../../systematics_ewk_Wmn.txt -b binStats_EWKWJetsWmn.txt --doEWK 1
python ../../printYields.py -c EWKWJetsWen -i hists_EWKWJetsWen.root -s ../../systematics_ewk_Wen.txt -b binStats_EWKWJetsWen.txt --doEWK 1
