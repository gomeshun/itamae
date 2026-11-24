import SlopeFinal60
import SlopeFinal61
import SlopeFinal62
import SlopeFinal63
import SlopeFinal64
import SlopeFinal65
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative60 : (⟨(-150897406068353456279/9903520314283042199192993792 : ℚ), (-58147393416030116623/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1000000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope60_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope60_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative60
theorem fullCompleteDerivative61 : (⟨(-78326629833284137637/9903520314283042199192993792 : ℚ), (-1951604173223991593/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope61_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope61_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative61
theorem fullCompleteDerivative62 : (⟨(-326088462944924382663/79228162514264337593543950336 : ℚ), (-1959625314401947503/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope62_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope62_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative62
theorem fullCompleteDerivative63 : (⟨(-170388912740548873711/79228162514264337593543950336 : ℚ), (-8181964344890822421/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope63_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope63_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative63
theorem fullCompleteDerivative64 : (⟨(-43961586640033673273/39614081257132168796771975168 : ℚ), (-4550768158819424647/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((10000000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope64_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope64_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative64
theorem fullCompleteDerivative65 : (⟨(-45333458911567433547/79228162514264337593543950336 : ℚ), (-293375927912706619/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((88913970502/5 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope65_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope65_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative65
