#import "../../../docs/manual/preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)

  = Material formulations

  #include "./docs/scTablePossibleOptionsMaterialTreatment.typ"

  #include "./docs/scChewsGaugeDecouplingAndRewrite.typ"
]

#standalone(
  title: "Material formulations",
  chapter,
)
