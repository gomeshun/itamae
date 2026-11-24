import SlopeFinal48
import SlopeFinal49
import SlopeFinal50
import SlopeFinal51
import SlopeFinal52
import SlopeFinal53
import SashimiEhFullDerivative
open FloatLib.Numerics SashimiEhFullDerivative
set_option maxRecDepth 100000
private theorem leftEdge : SashimiEhFullGrid.lowEdge = (1/100000000 : ℚ) := by decide +kernel
private theorem rightEdge : SashimiEhFullGrid.highEdge = (10000000000/1 : ℚ) := by decide +kernel
theorem fullCompleteDerivative48 : (⟨(-1197004858419641531847237/39614081257132168796771975168 : ℚ), (-97371058444523164862189/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((1000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope48_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope48_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative48
theorem fullCompleteDerivative49 : (⟨(-1294629179228375437236109/79228162514264337593543950336 : ℚ), (-3464053322956881162833/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/25000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope49_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope49_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative49
theorem fullCompleteDerivative50 : (⟨(-687815709528963625504011/79228162514264337593543950336 : ℚ), (-30888275972703618544837/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((316227766017/100000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope50_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope50_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative50
theorem fullCompleteDerivative51 : (⟨(-91961885795810781113267/19807040628566084398385987584 : ℚ), (-15797237628111197372397/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((56234132519/10000 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope51_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope51_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative51
theorem fullCompleteDerivative52 : (⟨(-196514179391829298134691/79228162514264337593543950336 : ℚ), (-8286799028622438012209/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((10000000/1 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope52_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope52_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative52
theorem fullCompleteDerivative53 : (⟨(-26247218352725396366995/19807040628566084398385987584 : ℚ), (-298213157456304121647/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (deriv fullS ((44456985251/2500 : ℚ) : ℝ)) := by
  apply replayed (by simpa only [leftEdge] using certifiedSpectrum) (by simpa only [leftEdge,rightEdge] using certifiedBandEnd) slope53_radiusCorrect (by decide +kernel) (by simpa only [leftEdge] using slope53_sum238) (by simpa only [leftEdge,rightEdge] using concrete.2.1)
  decide +kernel
#print axioms fullCompleteDerivative53
