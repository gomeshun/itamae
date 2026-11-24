import SlopeFinal72
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative72 : (⟨(-387289523595349633/79228162514264337593543950336 : ℚ), (-22300095943307683/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1000000000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope72_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope72_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative72
