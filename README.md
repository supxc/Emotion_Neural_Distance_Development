# Emotion Generation-Regulation Neural Distance: Example Code and Interaction Results

Associated manuscript: *Developmental Changes in the Neural Distance Between Emotion Generation and Emotion Regulation*.

## Important Data and Code Notice

**All R scripts in this package are illustrative/example scripts. No real participant-level data, raw MRI images, individual neural-distance values, or individual questionnaire responses are distributed.**

## Package Layout

```text
repository/
  README.md
  example_code/
    _r_scripts/                  9 original example R scripts
    _outputs/
      400P7N_Subcortexlabel.csv   ROI index-to-label lookup, not participant data
  interaction_results/
    age_by_sex_interaction.xlsx
    depression_group_x_age_MI.xlsx
    anxiety_group_x_age_MI.xlsx
```

## Script Index

| Script | Purpose | Status in this release |
|---|---|---|
| `01_simulate_data.R` | Generate synthetic demographics, behavior, symptom scores and 432 ROI distance variables | Demonstration input generator; not real data |
| `02_age_effect_analysis.R` | Global and ROI-level age associations | Example of the main age-analysis workflow |
| `03_brain_visualization.R` | ROI labeling and cortical-network/subcortical summary plots from age results | Example plotting workflow; not raw-image analysis |
| `04_behavior_effect_analysis.R` | Neural-distance associations with emotion generation/regulation behavior | Example behavioral workflow |
| `05_behavior_sig_visualization.R` | Label and summarize significant behavior-related ROIs | Example plotting workflow |
| `06_anxiety_depression_analysis.R` | Continuous symptom-score associations with distance using MI | Additional exploratory example; not the age-by-group interaction analysis |
| `08_scatter_roi_behaviors.R` | Scatterplots for selected ROI-behavior associations | Example plotting workflow |
| `09_age_by_sex_interaction.R` | Age-by-sex interactions in global and ROI distance | Corresponding type of analysis for the archived sex workbook |
| `11_clinical_group_x_age_interaction_MI.R` | MI-pooled age-by-threshold-group interactions | Corresponding type of analysis for the primary depression/anxiety workbooks |
