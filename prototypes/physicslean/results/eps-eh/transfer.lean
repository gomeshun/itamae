import SashimiEisensteinHu
open SashimiEisensteinHu FloatLib.Numerics

set_option maxRecDepth 100000
set_option maxHeartbeats 20000000
set_option Elab.async false
private def bs : List RationalInterval := [⟨(26699890767307081769024311263/39614081257132168796771975168 : ℚ), (53399781534614163538048622527/79228162514264337593543950336 : ℚ)⟩, ⟨(12324380835552230292329058939/79228162514264337593543950336 : ℚ), (96284225277751799158820773/618970019642690137449562112 : ℚ)⟩, ⟨(33451890839356053650607445695/39614081257132168796771975168 : ℚ), (8362972709839013412651861425/9903520314283042199192993792 : ℚ)⟩, ⟨(207502781351668752835363457/19807040628566084398385987584 : ℚ), (415005562703337505670726915/39614081257132168796771975168 : ℚ)⟩, ⟨(5981163722631073034125647614581/39614081257132168796771975168 : ℚ), (1495290930657768258744693127971/9903520314283042199192993792 : ℚ)⟩, ⟨(3572739045732410219822274241/39614081257132168796771975168 : ℚ), (1786369522866205111104284633/19807040628566084398385987584 : ℚ)⟩, ⟨(30618818463129726333345229299/39614081257132168796771975168 : ℚ), (30618818463129726336224277561/39614081257132168796771975168 : ℚ)⟩, ⟨(46589924881685299294835871361/39614081257132168796771975168 : ℚ), (23294962440842649647624436579/19807040628566084398385987584 : ℚ)⟩, ⟨(21657275934125292625671335171/39614081257132168796771975168 : ℚ), (43314551868250590088361105333/79228162514264337593543950336 : ℚ)⟩, ⟨(308950870566318646605186128057/39614081257132168796771975168 : ℚ), (617901741132637293210372256281/79228162514264337593543950336 : ℚ)⟩, ⟨(35750720062004382165707804343/9903520314283042199192993792 : ℚ), (286005760496035057366599914977/79228162514264337593543950336 : ℚ)⟩]
private theorem certifiedBoxes : coefficientBoxes? 0.674 0.315 0.049 2.7255 96 40 = some bs := by decide +kernel
private noncomputable def model := prepareReal 0.674 0.315 0.049 2.7255
private def box0 : RationalInterval := ⟨(4951760111239050139707226961/4951760157141521099596496896 : ℚ), (79228161779824802333011877789/79228162514264337593543950336 : ℚ)⟩
private theorem point0 : box0.Contains (transferReal model (1/1000000 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1/1000000 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box1 : RationalInterval := ⟨(39610473143426511359326275703/39614081257132168796771975168 : ℚ), (19805236571713255697263257857/19807040628566084398385987584 : ℚ)⟩
private theorem point1 : box1.Contains (transferReal model (1/10000 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1/10000 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box2 : RationalInterval := ⟨(62020886429708618660636423183/79228162514264337593543950336 : ℚ), (62020886429708618698969895529/79228162514264337593543950336 : ℚ)⟩
private theorem point2 : box2.Contains (transferReal model (1/100 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1/100 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box3 : RationalInterval := ⟨(10600665514633286125121507369/79228162514264337593543950336 : ℚ), (5300332757316643072202670019/39614081257132168796771975168 : ℚ)⟩
private theorem point3 : box3.Contains (transferReal model (1/10 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1/10 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box4 : RationalInterval := ⟨(190491213230712980553604385/39614081257132168796771975168 : ℚ), (190491213230712980611533171/39614081257132168796771975168 : ℚ)⟩
private theorem point4 : box4.Contains (transferReal model (1/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box5 : RationalInterval := ⟨(7317981341574862724776283/79228162514264337593543950336 : ℚ), (7317981341574862727058985/79228162514264337593543950336 : ℚ)⟩
private theorem point5 : box5.Contains (transferReal model (10/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (10/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box6 : RationalInterval := ⟨(109536864391866324902947/79228162514264337593543950336 : ℚ), (109536864391866324938777/79228162514264337593543950336 : ℚ)⟩
private theorem point6 : box6.Contains (transferReal model (100/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (100/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box7 : RationalInterval := ⟨(18239726438218855487/79228162514264337593543950336 : ℚ), (18239726438218855497/79228162514264337593543950336 : ℚ)⟩
private theorem point7 : box7.Contains (transferReal model (10000/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (10000/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box8 : RationalInterval := ⟨(2552580111921813/79228162514264337593543950336 : ℚ), (1276290055960909/39614081257132168796771975168 : ℚ)⟩
private theorem point8 : box8.Contains (transferReal model (1000000/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (1000000/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
private def box9 : RationalInterval := ⟨(164059349899/39614081257132168796771975168 : ℚ), (164059349901/39614081257132168796771975168 : ℚ)⟩
private theorem point9 : box9.Contains (transferReal model (100000000/1 : ℚ)) :=
  rangeFast_sound (RationalInterval.contains_point (100000000/1 : ℚ))
    (fun _ _ => coefficientBoxes_sound certifiedBoxes) (bits := 96) (degree := 40)
    (by decide +kernel)
theorem concrete : box0.Contains (transferReal model (1/1000000 : ℚ)) ∧
  box1.Contains (transferReal model (1/10000 : ℚ)) ∧
  box2.Contains (transferReal model (1/100 : ℚ)) ∧
  box3.Contains (transferReal model (1/10 : ℚ)) ∧
  box4.Contains (transferReal model (1/1 : ℚ)) ∧
  box5.Contains (transferReal model (10/1 : ℚ)) ∧
  box6.Contains (transferReal model (100/1 : ℚ)) ∧
  box7.Contains (transferReal model (10000/1 : ℚ)) ∧
  box8.Contains (transferReal model (1000000/1 : ℚ)) ∧
  box9.Contains (transferReal model (100000000/1 : ℚ)) := ⟨point0, point1, point2, point3, point4, point5, point6, point7, point8, point9⟩
#print axioms concrete
