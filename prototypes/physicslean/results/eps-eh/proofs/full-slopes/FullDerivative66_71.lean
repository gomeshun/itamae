import SlopeFinal66
import SlopeFinal67
import SlopeFinal68
import SlopeFinal69
import SlopeFinal70
import SlopeFinal71
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative66 : (⟨(-5810898387466210995/19807040628566084398385987584 : ℚ), (-1158627888227595517/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope66_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope66_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative66
theorem fullCompleteDerivative67 : (⟨(-187552137581039435/1237940039285380274899124224 : ℚ), (-310134249465447201/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope67_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope67_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative67
theorem fullCompleteDerivative68 : (⟨(-3030325447483204725/39614081257132168796771975168 : ℚ), (-332300172454528525/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((100000000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope68_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope68_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative68
theorem fullCompleteDerivative69 : (⟨(-768299616049148615/19807040628566084398385987584 : ℚ), (-82155681642366147/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((177827941004/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope69_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope69_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative69
theorem fullCompleteDerivative70 : (⟨(-387284329485972437/19807040628566084398385987584 : ℚ), (-40753647378019975/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope70_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope70_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative70
theorem fullCompleteDerivative71 : (⟨(-195032227169377767/19807040628566084398385987584 : ℚ), (-21900240674785895/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((562341325190/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope71_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope71_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative71
