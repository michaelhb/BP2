rm @OUTPUT_DIR@/testplot.csv

delta_min=0.01
delta_max=0.4
delta_step=0.01
n_fields=3

@GEN_PLOT@ "N=25" $n_fields 25 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
@GEN_PLOT@ "N=50" $n_fields 50 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
@GEN_PLOT@ "N=100" $n_fields 100 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
@GEN_PLOT@ "N=200" $n_fields 200 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
@GEN_PLOT@ "N=300" $n_fields 300 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
cat @OUTPUT_DIR@/testplot.csv | python3 @SHOW_PLOT@ "varying number of grid points (N): $n_fields fields"