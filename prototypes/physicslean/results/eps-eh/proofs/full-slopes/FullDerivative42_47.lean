import SlopeFinal42
import SlopeFinal43
import SlopeFinal44
import SlopeFinal45
import SlopeFinal46
import SlopeFinal47
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative42 : (⟨(-98878441167803402321101421/79228162514264337593543950336 : ℚ), (-513986372297542175913317/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope42_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope42_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative42
theorem fullCompleteDerivative43 : (⟨(-1656204979686862843624271/2475880078570760549798248448 : ℚ), (-1125354688627162676666239/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope43_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope43_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative43
theorem fullCompleteDerivative44 : (⟨(-7142371905353527883117693/19807040628566084398385987584 : ℚ), (-1149489561192952980984581/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((100000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope44_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope44_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative44
theorem fullCompleteDerivative45 : (⟨(-1938341350785635710617555/9903520314283042199192993792 : ℚ), (-312167010509616793364729/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/250000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope45_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope45_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative45
theorem fullCompleteDerivative46 : (⟨(-8308875682020793750052323/79228162514264337593543950336 : ℚ), (-361449581422058721251527/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope46_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope46_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative46
theorem fullCompleteDerivative47 : (⟨(-4467681912471290587263773/79228162514264337593543950336 : ℚ), (-95393795420075265482149/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/100000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope47_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope47_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative47
