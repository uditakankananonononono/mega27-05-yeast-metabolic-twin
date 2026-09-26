"""Locked stress-to-constraint parameters (pre-registration Section 4).

Functional forms and nominal values are fixed before any outcome run.
Sources:
- Temperature: Cardinal Temperature Model with Inflection (Rosso et al. 1993),
  cardinal values anchored to S. cerevisiae estimates (Salvado et al. 2011,
  Food Microbiol., ethanol x temperature fitness; Salvado et al. 2013,
  BMC Evol. Biol. 13:106 / PMC3067424, genus cardinal temperatures).
- Ethanol: Levenspiel-type product inhibition, anchored to published
  S. cerevisiae inhibition kinetics (Biotechnol. Bioeng. 27:3,
  doi:10.1002/bit.260270311; yield decline doi:10.1002/bit.260400213).
- Osmotic: obligate glycerol-production drain for osmoregulation
  (Hohmann 2002, Microbiol. Mol. Biol. Rev. 66:300, osmotic stress
  physiology; Yoshikawa 2009 / GSE59659 context).
- Nitrogen: ammonium uptake bound scaling vs the base medium level.
Values marked ASSUMED have no direct published value; they are swept +-50%
per the pre-registered sensitivity rule and are never fitted to validation
outcomes.
"""
from __future__ import annotations
from dataclasses import dataclass, field

# --- reaction identifiers (yeast-GEM, pinned in-repo SBML sha256
# 30842b15eefb0ef7e36cbdea86a9efddfacf69a871c8b054165faa9af9f6c8eb) ---
GROWTH_RXN = "r_2111"        # growth
NGAM_RXN = "r_4046"          # non-growth associated maintenance (0.7 mmol ATP/gDW/h)
GLUCOSE_EX = "r_1714"
ETHANOL_EX = "r_1761"
GLYCEROL_EX = "r_1808"
AMMONIUM_EX = "r_1654"
OXYGEN_EX = "r_1992"
PHOSPHATE_EX = "r_2005"

BASE_AMMONIUM_LB = -1000.0    # complete_y7 base: ammonium is in the FREE set (lb -1000)
BASE_GLUCOSE_LB = -20.0       # complete_y7 default glucose uptake (media.complete_y7 glucose=20)
BASE_OXYGEN_LB = -0.5        # microaerobic fermentation base; documented modeling choice
NGAM_VALUE = 0.7             # model default maintenance (mmol ATP/gDW/h)
VIABILITY_FRACTION = 0.01    # growth/no-growth boundary: f_T x f_E below this = nonviable (Amendment 1)


@dataclass(frozen=True)
class TemperatureModel:
    """CTMI relative growth factor, normalized to 1.0 at the 30 C reference."""
    t_min: float = 5.0
    t_opt: float = 32.0
    t_max: float = 42.0      # sweep 40.0-44.0 per pre-reg sensitivity rule
    t_ref: float = 30.0

    def _ctmi(self, t: float) -> float:
        if t <= self.t_min or t >= self.t_max:
            return 0.0
        num = (t - self.t_max) * (t - self.t_min) ** 2
        den = (self.t_opt - self.t_min) * (
            (self.t_opt - self.t_min) * (t - self.t_opt)
            - (self.t_opt - self.t_max) * (self.t_opt + self.t_min - 2 * t)
        )
        return num / den

    def factor(self, t: float) -> float:
        ref = self._ctmi(self.t_ref)
        return self._ctmi(t) / ref if ref > 0 else 0.0


@dataclass(frozen=True)
class EthanolModel:
    """Levenspiel-type inhibition factor, 1.0 at 0% v/v, 0 at p_max."""
    p_max: float = 14.0      # % v/v; matches locked grid edge
    n: float = 1.5           # ASSUMED shape exponent; sweep 0.75-2.25

    def factor(self, p: float) -> float:
        if p >= self.p_max:
            return 0.0
        return (1.0 - p / self.p_max) ** self.n


@dataclass(frozen=True)
class OsmoticModel:
    """Obligatory glycerol drain per unit external osmolarity."""
    k_gly: float = 1.0       # ASSUMED mmol glycerol/gDW/h per M; sweep 0.5-2.0


@dataclass(frozen=True)
class StressParameters:
    temperature: TemperatureModel = field(default_factory=TemperatureModel)
    ethanol: EthanolModel = field(default_factory=EthanolModel)
    osmotic: OsmoticModel = field(default_factory=OsmoticModel)


LOCKED = StressParameters()
