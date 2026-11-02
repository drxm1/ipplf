#import "./docs/manual/preamble.typ": *

#import "./docs/manual/Units.typ" as unitsChap
#import "./docs/manual/SimplifyingAssumptions.typ" as simplifyingAssumptions
#import "./docs/manual/maxwellsEqsForOurMaterialModel.typ" as maxwellsEqsForOurMaterialModel
#import "./docs/manual/theVacuumUpdateInIPPL.typ" as theVacuumUpdateInIPPL
#import "./core/base/levelContract.typ" as theLevelContract
#import "./core/solver/formulation/formulations.typ" as formulations

#manual[
  #unitsChap.chapter
  #simplifyingAssumptions.chapter
  #maxwellsEqsForOurMaterialModel.chapter
  #theVacuumUpdateInIPPL.chapter
  #theLevelContract.chapter
  #formulations.chapter
]
