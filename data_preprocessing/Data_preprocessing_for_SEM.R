# Code for new model: ELS_total<->NIH_reading, Ab_PA<->NIH_total
# Data Preprocessing

# 1. Make csv files that consists of:
# - subjectkey, age, sex, high.educ,  abcd_site,  marriage, BMI, income
# - ELS_total, Ab_PA
# - PC1~16
# - NIH_total, NIH_reading, NIH_fluid, NIH_crystal
# 
# 2. Make interaction term
# - ELS*EA
# - Ab_PA*EA
# - ELS*PC1
# - Ab_PA*PC1

demo_list <- c('subjectkey', 'age', 'sex', "race.ethnicity", "high.educ", "income", "married", "abcd_site", "BMI")
els_list <- c("ELS_total" , "Ab_PA")
