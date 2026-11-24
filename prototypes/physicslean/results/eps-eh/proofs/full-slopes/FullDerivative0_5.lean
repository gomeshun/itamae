import SlopeFinal0
import SlopeFinal1
import SlopeFinal2
import SlopeFinal3
import SlopeFinal4
import SlopeFinal5
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative0 : (⟨(-4556703251604488282305588955045945625/39614081257132168796771975168 : ℚ), (-40758692964341533892216653182305239/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/1000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope0_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope0_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative0
theorem fullCompleteDerivative1 : (⟨(-5036367016965769471527272399347166589/79228162514264337593543950336 : ℚ), (-168460260808784995183357775745935381/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope1_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope1_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative1
theorem fullCompleteDerivative2 : (⟨(-702315714254873830580310785792742137/19807040628566084398385987584 : ℚ), (-24762257546132651122875657171064989/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope2_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope2_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative2
theorem fullCompleteDerivative3 : (⟨(-769158087086870536334566803473290901/39614081257132168796771975168 : ℚ), (-3568648198087659363588471937244015/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10000000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope3_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope3_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative3
theorem fullCompleteDerivative4 : (⟨(-850301419689481947740798779735053333/79228162514264337593543950336 : ℚ), (-14833610717623405313936433055986203/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/100000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope4_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope4_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative4
theorem fullCompleteDerivative5 : (⟨(-470223872322443167941371110160226495/79228162514264337593543950336 : ℚ), (-7949432137871110405018097336673603/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/2500000000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope5_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope5_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative5
