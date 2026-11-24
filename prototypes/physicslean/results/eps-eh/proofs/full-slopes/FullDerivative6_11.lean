import SlopeFinal6
import SlopeFinal7
import SlopeFinal8
import SlopeFinal9
import SlopeFinal10
import SlopeFinal11
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative6 : (⟨(-260686355117219712722313592446904857/79228162514264337593543950336 : ℚ), (-2400930176381863896367481206464479/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope6_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope6_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative6
theorem fullCompleteDerivative7 : (⟨(-143045367943476555581610852725302035/79228162514264337593543950336 : ℚ), (-5283148599717745982079689799146537/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope7_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope7_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative7
theorem fullCompleteDerivative8 : (⟨(-78797662605024299900739183465557965/79228162514264337593543950336 : ℚ), (-2713436884197173598552101869978965/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/10000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope8_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope8_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative8
theorem fullCompleteDerivative9 : (⟨(-10949877872248422664017903033778753/19807040628566084398385987584 : ℚ), (-1536470587130871942850490616495167/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/250000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope9_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope9_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative9
theorem fullCompleteDerivative10 : (⟨(-6000491827822257478944334774663967/19807040628566084398385987584 : ℚ), (-453058888515247799413503160029313/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope10_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope10_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative10
theorem fullCompleteDerivative11 : (⟨(-3308311277340972812610678957306771/19807040628566084398385987584 : ℚ), (-479525139039228009442769472746199/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/100000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope11_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope11_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative11
