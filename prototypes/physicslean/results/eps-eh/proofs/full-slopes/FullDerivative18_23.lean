import SlopeFinal18
import SlopeFinal19
import SlopeFinal20
import SlopeFinal21
import SlopeFinal22
import SlopeFinal23
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative18 : (⟨(-50547932516580779974716737179221/19807040628566084398385987584 : ℚ), (-7605400780966687475491108728971/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope18_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope18_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative18
theorem fullCompleteDerivative19 : (⟨(-27738001146899809029177520923399/19807040628566084398385987584 : ℚ), (-1950653114499475826551784617777/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope19_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope19_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative19
theorem fullCompleteDerivative20 : (⟨(-30780651117700418585306285939207/39614081257132168796771975168 : ℚ), (-2255197023407650725191093466503/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/10 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope20_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope20_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative20
theorem fullCompleteDerivative21 : (⟨(-16734322195236404481681202769673/39614081257132168796771975168 : ℚ), (-1300178480229619570133023180023/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/250000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope21_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope21_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative21
theorem fullCompleteDerivative22 : (⟨(-9185819303995741191294732067661/39614081257132168796771975168 : ℚ), (-169276668316169120448922957463/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope22_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope22_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative22
theorem fullCompleteDerivative23 : (⟨(-5044112438299997374002643785555/39614081257132168796771975168 : ℚ), (-178864080579155642105919787285/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/100000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope23_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope23_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative23
