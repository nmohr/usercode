./manualStack.py -C 8TeVconfig/configPlotManual -T log10
./manualStack.py -C 8TeVconfig/configPlotManual
./manualStack.py -C 8TeVconfig/configPlotManual -S True
python contourPlots.py
root ccc_VZ.C
