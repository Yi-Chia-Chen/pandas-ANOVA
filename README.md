# Python ANOVA with Pandas (pandas-ANOVA)

This is a python 3 script using pandas to perform between-group ANOVA by Yi-Chia Chen.

## Instructions
The data to be analyzed should be organized and saved in "input.txt" in the same folder of the script. It should be a chart with rows of subjects and columns of variables (separated by tabs and line breaks). The last columns will be treated as dependent variables and all other columns independent variables. Column titles should be in the first row and these will be printed as factor names in the output ANOVA table.

## Version History
- 0.1.0b (2020.02.05): First beta version
- 0.1.1b (2020.02.06): replace for loops with reduce() and comprehension

## Planned Improvements
- Add within-subject ANOVA
- Automatically analyze simple main effect for siginificant interaction effects
- Develop testing tool sets
