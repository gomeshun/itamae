import SlopeFinal30
import SlopeFinal31
import SlopeFinal32
import SlopeFinal33
import SlopeFinal34
import SlopeFinal35
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative30 : (⟨(-147658759153268361915330744501/79228162514264337593543950336 : ℚ), (-671302654509513093282801327/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope30_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope30_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative30
theorem fullCompleteDerivative31 : (⟨(-81251943163059118880544757185/79228162514264337593543950336 : ℚ), (-3160152519096023060439016383/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope31_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope31_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative31
theorem fullCompleteDerivative32 : (⟨(-687290332424006292618994659/1237940039285380274899124224 : ℚ), (-889060823483182637271519315/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((100/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope32_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope32_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative32
theorem fullCompleteDerivative33 : (⟨(-23981953427205823073488530411/79228162514264337593543950336 : ℚ), (-913101894120595234970033623/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/250000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope33_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope33_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative33
theorem fullCompleteDerivative34 : (⟨(-13095048364453162858097943269/79228162514264337593543950336 : ℚ), (-490951053595079716218547571/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope34_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope34_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative34
theorem fullCompleteDerivative35 : (⟨(-7132492569854539923568233201/79228162514264337593543950336 : ℚ), (-289572962704405787703827095/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/100000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope35_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope35_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative35
