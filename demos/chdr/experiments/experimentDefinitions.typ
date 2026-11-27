#import "../docs/manual/preamble.typ": *
#import "./E15_homogeneous_long_run_stability/E15-S.typ" as e15s

#let chapter = [
  = Experiment definitions <experiment-definitions>

  // Every included file starts with a first-level heading for its experiment.
  // The offset makes each experiment a section of this chapter.
  #set heading(offset: 1)

  #include "./E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E01.typ"

  #include "./E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E02.typ"

  #include "./E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E03.typ"

  #include "./E01E02E03E04_wave_pulses_in_homogeneous_dielectric/E04.typ"

  #include "./E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E05.typ"

  #include "./E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E06.typ"

  #include "./E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E07.typ"

  #include "./E05E06E07E08_gaussian_wave_pulses_in_homogeneous_dielectric/E08.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E09.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E10.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E11.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E12.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E13.typ"

  #include "./E09E10E11E12E13E14_comparison_with_ippl_standard_solver/E14.typ"

  #include "./E16E17_normal_incidence_on_planar_dielectric_interface/E16.typ"

  #include "./E16E17_normal_incidence_on_planar_dielectric_interface/E17.typ"

  #include "./E18E19_oblique_incidence_on_planar_dielectric_interface/E18.typ"

  #include "./E18E19_oblique_incidence_on_planar_dielectric_interface/E19.typ"

  #include "./E20_field_free_potential_pulse_on_planar_dielectric_interface/E20.typ"

  #e15s.section

  #include "./E16E17_normal_incidence_on_planar_dielectric_interface/E16S.typ"

  #include "./E18E19_oblique_incidence_on_planar_dielectric_interface/E18S.typ"

  #include "./E18E19_oblique_incidence_on_planar_dielectric_interface/E19S.typ"
]
