rm @OUTPUT_DIR@/testplot.csv

delta_min=0.01
delta_max=0.2
delta_step=0.005
n_fields=3

@GEN_PLOT@ "BP2, N=200" $n_fields 200 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
@GEN_PLOT_CT@ "CosmoTransitions" $n_fields $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv

# PLOT 1 - 3 FIELDS
# delta_min=0.01
# delta_max=0.2
# delta_step=0.005
# n_fields=3

# @GEN_PLOT@ "BP2, N=200" $n_fields 200 $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv
# @GEN_PLOT_CT@ "CosmoTransitions" $n_fields $delta_min $delta_max $delta_step | tee -a @OUTPUT_DIR@/testplot.csv

cat @OUTPUT_DIR@/testplot.csv | python3 @SHOW_PLOT@ "Comparing BP2 to CosmoTransitions: $n_fields fields"