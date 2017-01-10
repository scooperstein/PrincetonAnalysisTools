## Run this after all VHbbAnalysis jobs are finished

mkdir -p Data
mv Run* Data
mkdir -p TT_Powheg
mv TT_powheg TT_Powheg
mkdir -p Signal
mv WH* WplusH* WminusH* ZH* Signal
hadd -f output_data.root Data/Run*/*.root
hadd -f output_mc.root */*.root
hadd -f output_ttpowheg.root TT_Powheg/TT_powheg/*.root
hadd -f output_signal.root Signal/*/*.root
