#import "../../../docs/manual/preamble.typ": *

#let chapter = [
  #show: equate.with(sub-numbering: false)

  = Material formulations

  #include "./docs/scTablePossibleOptionsMaterialTreatment.typ"

  #include "./docs/scChewsGaugeDecouplingAndRewrite.typ"

  #include "./docs/scO1VSImplementation.typ"
]

#standalone(
  title: "Material formulations",
  chapter,
)
