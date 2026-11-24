import SlopeFinal36
import SlopeFinal37
import SlopeFinal38
import SlopeFinal39
import SlopeFinal40
import SlopeFinal41
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative36 : (⟨(-3867428778566563479506671781/79228162514264337593543950336 : ℚ), (-156134112705032962275424505/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope36_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope36_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative36
theorem fullCompleteDerivative37 : (⟨(-2096764672800833920137134853/79228162514264337593543950336 : ℚ), (-39836495974711687529892365/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope37_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope37_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative37
theorem fullCompleteDerivative38 : (⟨(-1148720844247690126484866051/79228162514264337593543950336 : ℚ), (-45208992629390433446783309/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope38_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope38_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative38
theorem fullCompleteDerivative39 : (⟨(-154547480502641251384756657/19807040628566084398385987584 : ℚ), (-12960031858928255623068437/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope39_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope39_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative39
theorem fullCompleteDerivative40 : (⟨(-167629067057481242527046635/39614081257132168796771975168 : ℚ), (-210586048607764725417017/1237940039285380274899124224 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((10000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope40_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope40_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative40
theorem fullCompleteDerivative41 : (⟨(-181606850630727096448972067/79228162514264337593543950336 : ℚ), (-6988023357494340537001129/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/2500000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope41_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope41_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative41
