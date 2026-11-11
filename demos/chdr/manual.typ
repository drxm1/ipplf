#import "./docs/manual/preamble.typ": *

#import "./docs/manual/Units.typ" as unitsChap
#import "./docs/manual/SimplifyingAssumptions.typ" as simplifyingAssumptions
#import "./docs/manual/maxwellsEqsForOurMaterialModel.typ" as maxwellsEqsForOurMaterialModel
#import "./docs/manual/theVacuumUpdateInIPPL.typ" as theVacuumUpdateInIPPL
#import "./core/base/levelContract.typ" as theLevelContract
#import "./core/solver/formulation/formulations.typ" as formulations
#import "./opt-in/diagnostics/diagnostics.typ" as diagnostics
#import "./experiments/experimentDefinitions.typ" as experimentDefinitions
#import "./docs/manual/application.typ" as application

#manual(contents: true)[
  _\*Note: Parts of this manual were written with the help of ai agents._
  #unitsChap.chapter
  #simplifyingAssumptions.chapter
  #maxwellsEqsForOurMaterialModel.chapter
  #theVacuumUpdateInIPPL.chapter
  #theLevelContract.chapter
  #formulations.chapter
  #diagnostics.chapter
  #experimentDefinitions.chapter
  #application.chapter
]
