import SlopeFinal54
import SlopeFinal55
import SlopeFinal56
import SlopeFinal57
import SlopeFinal58
import SlopeFinal59
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative54 : (⟨(-55590264458756555009437/79228162514264337593543950336 : ℚ), (-1275965854160257423581/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/10000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope54_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope54_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative54
theorem fullCompleteDerivative55 : (⟨(-29422486300140190193599/79228162514264337593543950336 : ℚ), (-1287486801814067688767/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/1000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope55_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope55_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative55
theorem fullCompleteDerivative56 : (⟨(-7853364893582653254821/39614081257132168796771975168 : ℚ), (-703986181311494079235/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((100000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope56_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope56_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative56
theorem fullCompleteDerivative57 : (⟨(-515332690363432329403/4951760157141521099596496896 : ℚ), (-395852697823339332143/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/250 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope57_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope57_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative57
theorem fullCompleteDerivative58 : (⟨(-1087390208668462832751/19807040628566084398385987584 : ℚ), (-203514733565663378333/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/1000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope58_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope58_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative58
theorem fullCompleteDerivative59 : (⟨(-2284578648344494580029/79228162514264337593543950336 : ℚ), (-102814097474542037491/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/100 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope59_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope59_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative59
