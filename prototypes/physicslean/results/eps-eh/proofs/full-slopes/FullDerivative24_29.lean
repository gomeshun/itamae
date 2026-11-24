import SlopeFinal24
import SlopeFinal25
import SlopeFinal26
import SlopeFinal27
import SlopeFinal28
import SlopeFinal29
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative24 : (⟨(-2780390685821421016222718392909/39614081257132168796771975168 : ℚ), (-107116946525014331036758679781/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope24_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope24_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative24
theorem fullCompleteDerivative25 : (⟨(-1511848474891226255744000428877/39614081257132168796771975168 : ℚ), (-117762509939155055398684685397/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope25_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope25_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative25
theorem fullCompleteDerivative26 : (⟨(-413633567167999667746068597175/19807040628566084398385987584 : ℚ), (-60404381368554196231650216137/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope26_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope26_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative26
theorem fullCompleteDerivative27 : (⟨(-455542817347149281136292375871/39614081257132168796771975168 : ℚ), (-33557104641411963330308779303/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10000000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope27_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope27_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative27
theorem fullCompleteDerivative28 : (⟨(-123865044629770144131464379513/19807040628566084398385987584 : ℚ), (-19698605234921178921927243203/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((10/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope28_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope28_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative28
theorem fullCompleteDerivative29 : (⟨(-135483135130772969002119305213/39614081257132168796771975168 : ℚ), (-5222437225028622531123140473/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/2500000000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope29_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope29_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative29
