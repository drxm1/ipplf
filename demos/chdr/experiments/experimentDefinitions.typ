#import "../docs/manual/preamble.typ": *

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
]
