module Abstracts
	
export ScatteringLaw, AnalyticalScatteringLaw, PhaseFunction

abstract ScatteringLaw
abstract AnalyticalScatteringLaw <: ScatteringLaw
abstract PhaseFunction

end # module
