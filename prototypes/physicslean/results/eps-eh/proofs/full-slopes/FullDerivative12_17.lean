import SlopeFinal12
import SlopeFinal13
import SlopeFinal14
import SlopeFinal15
import SlopeFinal16
import SlopeFinal17
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative12 : (⟨(-1820970956982400200598871156485727/19807040628566084398385987584 : ℚ), (-123966616178349439330834568173571/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/1000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope12_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope12_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative12
theorem fullCompleteDerivative13 : (⟨(-4043699738340248241156066781254815/79228162514264337593543950336 : ℚ), (-74495166237970761648533493473423/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope13_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope13_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative13
theorem fullCompleteDerivative14 : (⟨(-1104609577739563230722574295702327/39614081257132168796771975168 : ℚ), (-41861673457088707522757573896153/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope14_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope14_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative14
theorem fullCompleteDerivative15 : (⟨(-1214744236786182063504513679188433/79228162514264337593543950336 : ℚ), (-43187580622523472123893433651429/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope15_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope15_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative15
theorem fullCompleteDerivative16 : (⟨(-167891510628614244098992773665131/19807040628566084398385987584 : ℚ), (-11799040949461003146384636378375/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/100 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope16_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope16_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative16
theorem fullCompleteDerivative17 : (⟨(-368568616217814121209493204267609/79228162514264337593543950336 : ℚ), (-14100586391965294782744318948433/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/2500000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope17_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope17_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative17
