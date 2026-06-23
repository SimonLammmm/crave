#### Shiny server: initialisation ####

defaults_gq_cutoff_choices = c(1e-30, 1e-20, 1e-15, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.25)
defaults_gq_cutoff_selected = 1e-4
defaults_corr_cutoff_min = 0.4
defaults_corr_cutoff_max = 1
defaults_corr_cutoff_step = 0.05
defaults_corr_cutoff_value = 0.7

# Initialise and reset explore customise
defaults_explore_customise_height = 900
defaults_explore_customise_width = 0
defaults_explore_customise_width_auto = T
defaults_explore_customise_point_colour = "#ff0087"
defaults_explore_customise_goi_colour = "#007cff"
defaults_explore_customise_x_axis_log = "Automatic"
defaults_explore_customise_y_axis_log = "Automatic"
defaults_explore_customise_y_equals_x = T
defaults_explore_customise_density = F
defaults_explore_customise_rug_x = F
defaults_explore_customise_rug_y = F
explore_customise_rug_x = T
explore_customise_rug_y = T


# Initialise and reset correlate customise
defaults_correlate_customise_height = 900
defaults_correlate_customise_width = 900
defaults_correlate_customise_width_auto = T
defaults_correlate_customise_point_colour = "#ff0087"
defaults_correlate_customise_genequery_fillscheme = "Treatment"
defaults_correlate_customise_heatmap_high = "#ff0087"
defaults_correlate_customise_heatmap_mid = "#ffffff"
defaults_correlate_customise_heatmap_low = "#00007e"