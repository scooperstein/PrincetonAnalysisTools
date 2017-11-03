hadd -f hists_WmnHighPt.root hists_WmnHighPt_*.root
hadd -f hists_WenHighPt.root hists_WenHighPt_*.root
hadd -f hists_ttWmn.root hists_ttWmn_*.root
hadd -f hists_ttWen.root hists_ttWen_*.root
hadd -f hists_wlfWmn.root hists_wlfWmn_*.root
hadd -f hists_wlfWen.root hists_wlfWen_*.root
hadd -f hists_whfWmnLow.root hists_whfWmnLow_*.root
hadd -f hists_whfWenLow.root hists_whfWenLow_*.root
hadd -f hists_whfWmnHigh.root hists_whfWmnHigh_*.root
hadd -f hists_whfWenHigh.root hists_whfWenHigh_*.root
cat binStats_WmnHighPt_*.txt > binStats_WmnHighPt.txt
cat binStats_WenHighPt_*.txt > binStats_WenHighPt.txt
cat binStats_ttWmn_*.txt > binStats_ttWmn.txt
cat binStats_ttWen_*.txt > binStats_ttWen.txt
cat binStats_wlfWmn_*.txt > binStats_wlfWmn.txt
cat binStats_wlfWen_*.txt > binStats_wlfWen.txt
cat binStats_whfWmnLow_*.txt > binStats_whfWmnLow.txt
cat binStats_whfWenLow_*.txt > binStats_whfWenLow.txt
cat binStats_whfWmnHigh_*.txt > binStats_whfWmnHigh.txt
cat binStats_whfWenHigh_*.txt > binStats_whfWenHigh.txt

python ../../printYields.py -c WmnHighPt -i hists_WmnHighPt.root -s ../../systematics_Wmn_highsplitBT.txt -b binStats_WmnHighPt.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln    -bnt 30   
python ../../printYields.py -c WenHighPt -i hists_WenHighPt.root -s ../../systematics_Wen_highsplitBT.txt -b binStats_WenHighPt.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln    -bnt 30  
python ../../printYields.py -c ttWmn -i hists_ttWmn.root -s ../../systematics_Wmn_highsplitBT.txt -b binStats_ttWmn.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln    
python ../../printYields.py -c ttWen -i hists_ttWen.root -s ../../systematics_Wen_highsplitBT.txt -b binStats_ttWen.txt -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln   
python ../../printYields.py -c wlfWmn -i hists_wlfWmn.root -s ../../systematics_Wmn_highsplitBT.txt -b binStats_wlfWmn.txt   -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln     
python ../../printYields.py -c wlfWen -i hists_wlfWen.root -s ../../systematics_Wen_highsplitBT.txt -b binStats_wlfWen.txt   -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln   
python ../../printYields.py -c whfWmnLow -i hists_whfWmnLow.root -s ../../systematics_Wmn_highsplitBT.txt -b binStats_whfWmnLow.txt   -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln          
python ../../printYields.py -c whfWenLow -i hists_whfWenLow.root -s ../../systematics_Wen_highsplitBT.txt -b binStats_whfWenLow.txt     -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln   
python ../../printYields.py -c whfWmnHigh -i hists_whfWmnHigh.root -s ../../systematics_Wmn_highsplitBT.txt -b binStats_whfWmnHigh.txt   -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln       
python ../../printYields.py -c whfWenHigh -i hists_whfWenHigh.root -s ../../systematics_Wen_highsplitBT.txt -b binStats_whfWenHigh.txt     -r TT_Wln,Wj0b_Wln,Wj1b_Wln,Wj2b_Wln     
