import SashimiEhHost
open SashimiEhHost SashimiEisensteinIntegral FloatLib.Numerics
set_option maxHeartbeats 0
set_option maxRecDepth 100000
private def cr := (⟨(4249493329856275766452151469309587/5192296858534827628530496329220096 : ℚ), (4249493336317676395081585996139171/5192296858534827628530496329220096 : ℚ)⟩ : RationalInterval)
private def pr := (⟨(3181581070066054162317212445488541/2596148429267413814265248164610048 : ℚ), (3181581076224908531977775316882823/2596148429267413814265248164610048 : ℚ)⟩ : RationalInterval)
private def sval := (⟨(-622293623953105633001/1180591620717411303424 : ℚ), (-2489174495812286061699/4722366482869645213696 : ℚ)⟩ : RationalInterval)
private def dsval := (⟨(130161345812930201545613451810539850515476413218770065485866510963807581575698122505262264775398815064533906258251141886930284725976591969359476886461972593377603234713266211316716135409379591727/104125088364876246092836878161019485005525250424424861892591118466515540704970652056659662052902809278250139161590038890680832783054336764502898808954649520763424917226667962959237355730638000000 : ℚ), (3645557269147406008416397093736193902632682412376861201827541986064216060072850894049023708532127451720039813594329285460547438704966425603479620221780096431840727300549687901368518340521007286281/1315878635199353295995644442564463977699538131263870967618648798109063400131620041598478489257496336817093364629893178453701021658477939648424134333562131927275452108005993068218231201171875000000 : ℚ)⟩ : RationalInterval)
private theorem hc : cr.Contains childRadius := radius_check_sound (bits:=128) (degree:=40) (by decide +kernel)
private def heathRecord : SashimiHeath.Row := ⟨(⟨(1/1 : ℚ), (1/1 : ℚ)⟩ : RationalInterval), sval, (⟨(7961909890118221830291/4722366482869645213696 : ℚ), (1990477472529555457573/1180591620717411303424 : ℚ)⟩ : RationalInterval)⟩
private theorem hSlope : sval.Contains (deriv (SashimiHeath.D matched.om) 0) := by
  have h := (SashimiHeath.check_sound (m:=matched.om) (z:=0) (n:=90) (out:=heathRecord) (by decide +kernel)).2.1
  simpa only [heathRecord, Rat.cast_zero] using h
theorem host_M200_0 (hDS : dsval.Contains deltaS) : (⟨(1000000000000/1 : ℚ), (1000000000000/1 : ℚ)⟩ : RationalInterval).Contains (target .m200 (0/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_0 (hDS : dsval.Contains deltaS) : (⟨(-33027195152470619607349041083106954995659/39614081257132168796771975168 : ℚ), (-31262750784409978817915265489929230988275/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (0/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_1 (hDS : dsval.Contains deltaS) : (⟨(65873185668821601330344506255209213891553/79228162514264337593543950336 : ℚ), (72254613449399543180931548596054528390231/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_1 (hDS : dsval.Contains deltaS) : (⟨(-63394309009501376158681264113966319518073/79228162514264337593543950336 : ℚ), (-13987933903822480229436258617148579558927/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_2 (hDS : dsval.Contains deltaS) : (⟨(54379977834007927494445301753806381968621/79228162514264337593543950336 : ℚ), (65198630445510496054739609926249404153309/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_2 (hDS : dsval.Contains deltaS) : (⟨(-59479478764224833719293337469725234617577/79228162514264337593543950336 : ℚ), (-23795945731771339489082361142797257613111/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_3 (hDS : dsval.Contains deltaS) : (⟨(44671546287515265199001655738406221427143/79228162514264337593543950336 : ℚ), (58402024494126410789157968807033068779751/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_3 (hDS : dsval.Contains deltaS) : (⟨(-13781320784398075874112836972821318314235/19807040628566084398385987584 : ℚ), (-19810850780574530112326035884790115858055/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_4 (hDS : dsval.Contains deltaS) : (⟨(18282171423201848083387024554517573624975/39614081257132168796771975168 : ℚ), (1626062601264377643744076526704476955751/2475880078570760549798248448 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_4 (hDS : dsval.Contains deltaS) : (⟨(-1582909323967220262685671654236323945567/2475880078570760549798248448 : ℚ), (-8156817995962736305841412476857313626721/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_5 (hDS : dsval.Contains deltaS) : (⟨(14923102928334329471799446621380730775991/39614081257132168796771975168 : ℚ), (11542687447934451773164333795012839413357/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_5 (hDS : dsval.Contains deltaS) : (⟨(-46222755762600225373828161421072959972989/79228162514264337593543950336 : ℚ), (-6673790310444165240758730236617950892181/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_6 (hDS : dsval.Contains deltaS) : (⟨(24309583429390215545719527742946165422651/79228162514264337593543950336 : ℚ), (10208990777251094639341728173480425696841/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_6 (hDS : dsval.Contains deltaS) : (⟨(-20965737745426661452774689672913040613699/39614081257132168796771975168 : ℚ), (-10875555915456899583386534041410089510559/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_7 (hDS : dsval.Contains deltaS) : (⟨(4941320338324178479280590226384084228791/19807040628566084398385987584 : ℚ), (18011626809484796268850587793224380196091/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_7 (hDS : dsval.Contains deltaS) : (⟨(-9460706106972849947446169633404244081489/19807040628566084398385987584 : ℚ), (-8836426805124899047879343185713603417845/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_8 (hDS : dsval.Contains deltaS) : (⟨(8023586556734161067136150167246785413781/39614081257132168796771975168 : ℚ), (31709198072754215481666900520520365462749/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (2/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_8 (hDS : dsval.Contains deltaS) : (⟨(-4249667237622676649963128726000131602199/9903520314283042199192993792 : ℚ), (-3582791071246496092912287049468329698369/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (2/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_9 (hDS : dsval.Contains deltaS) : (⟨(13012630429969071360230750466554754934401/79228162514264337593543950336 : ℚ), (1741326069088428843360538755028486256381/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (target .m200 (9/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_9 (hDS : dsval.Contains deltaS) : (⟨(-15209185081176500050566117164267910531413/39614081257132168796771975168 : ℚ), (-5802742387442824833185514233448157136017/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (9/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_10 (hDS : dsval.Contains deltaS) : (⟨(10540984107023050612494268822440031096785/79228162514264337593543950336 : ℚ), (763825713545613993103856439263813295327/2475880078570760549798248448 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_10 (hDS : dsval.Contains deltaS) : (⟨(-27116353780634660753179757656579458149973/79228162514264337593543950336 : ℚ), (-4694702543895559449032384265015489801339/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_11 (hDS : dsval.Contains deltaS) : (⟨(8531176821303385169453818845108172011065/79228162514264337593543950336 : ℚ), (10707310359789094896312607840201848091403/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (11/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_11 (hDS : dsval.Contains deltaS) : (⟨(-24092174772359707755745181330000426107057/79228162514264337593543950336 : ℚ), (-3795898028100936850495527967541537570443/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (11/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_12 (hDS : dsval.Contains deltaS) : (⟨(6899197964758223923899446126307875063419/79228162514264337593543950336 : ℚ), (18740153990823303078604441326306188379121/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_12 (hDS : dsval.Contains deltaS) : (⟨(-5334964211840313984219159338547267315873/19807040628566084398385987584 : ℚ), (-3068021010106181252995232811840548485853/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_13 (hDS : dsval.Contains deltaS) : (⟨(2787797198593959898242545041930445818921/39614081257132168796771975168 : ℚ), (16383007561301428142171711917278778477541/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (13/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_13 (hDS : dsval.Contains deltaS) : (⟨(-9424331644698701779675421363527802761347/39614081257132168796771975168 : ℚ), (-2479253700943716563701759859852384923761/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (13/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_14 (hDS : dsval.Contains deltaS) : (⟨(4503193675815690287294208187275851355575/79228162514264337593543950336 : ℚ), (14309432417051197677495621609597531347823/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_14 (hDS : dsval.Contains deltaS) : (⟨(-16604714543732914443270520585736786812475/79228162514264337593543950336 : ℚ), (-2003390240634515913868614084257909091791/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_15 (hDS : dsval.Contains deltaS) : (⟨(3635092207861966220477624966660463748391/79228162514264337593543950336 : ℚ), (12488258000087734111699027707203631076843/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (15/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_15 (hDS : dsval.Contains deltaS) : (⟨(-1824026666702713553890279981035718459729/9903520314283042199192993792 : ℚ), (-809488346765711678726283227490510289047/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (15/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_16 (hDS : dsval.Contains deltaS) : (⟨(2932916331321175783656703832535080914065/79228162514264337593543950336 : ℚ), (10891002260834352153705170267920759413407/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (4/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_16 (hDS : dsval.Contains deltaS) : (⟨(-12794358589839818688829204319339828882225/79228162514264337593543950336 : ℚ), (-81782938542834746945384484260007835891/4951760157141521099596496896 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (4/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_17 (hDS : dsval.Contains deltaS) : (⟨(591335451757120340685133404991370227789/19807040628566084398385987584 : ℚ), (1186481523841952340556145795968740646095/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .m200 (17/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_17 (hDS : dsval.Contains deltaS) : (⟨(-2673589821009975979127612317762546933117/19807040628566084398385987584 : ℚ), (-1106235397085051406958980739994450065423/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (17/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_18 (hDS : dsval.Contains deltaS) : (⟨(476711718508086218089912708606487009147/19807040628566084398385987584 : ℚ), (8267563628923968784442075334649016703171/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (9/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_18 (hDS : dsval.Contains deltaS) : (⟨(-9380407435739691848538204071460806530925/79228162514264337593543950336 : ℚ), (-111349038142505912636191676920372278911/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (9/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_19 (hDS : dsval.Contains deltaS) : (⟨(192083803242053917844884123679453849849/9903520314283042199192993792 : ℚ), (3598656596558884506905184840106823266985/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (19/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_19 (hDS : dsval.Contains deltaS) : (⟨(-8220439651202759875104627036418252748353/79228162514264337593543950336 : ℚ), (-22406840136248523095104102588477456225/2475880078570760549798248448 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (19/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_20 (hDS : dsval.Contains deltaS) : (⟨(618973451256944373259382960136464693659/39614081257132168796771975168 : ℚ), (782815555549141887643993716226107007645/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_20 (hDS : dsval.Contains deltaS) : (⟨(-3598936877293835095223997282442681574061/39614081257132168796771975168 : ℚ), (-288470470414207610831282259668995562635/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_21 (hDS : dsval.Contains deltaS) : (⟨(498495688763933901399541547590545923143/39614081257132168796771975168 : ℚ), (5446683078322473318906486530232465051791/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (21/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_21 (hDS : dsval.Contains deltaS) : (⟨(-1574399246905139230007220350391173694161/19807040628566084398385987584 : ℚ), (-464084957271551334796676874369415063197/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (21/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_22 (hDS : dsval.Contains deltaS) : (⟨(100338880175538985925432307483407355041/9903520314283042199192993792 : ℚ), (4735150910072053808093247288595984442005/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (11/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_22 (hDS : dsval.Contains deltaS) : (⟨(-5505913042800110548842297932868177575613/79228162514264337593543950336 : ℚ), (-93300832377360119753766441321609016385/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (11/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_23 (hDS : dsval.Contains deltaS) : (⟨(646122251730729027421784642433851863107/79228162514264337593543950336 : ℚ), (4114984953640999706568184915682554468385/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (23/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_23 (hDS : dsval.Contains deltaS) : (⟨(-4810467587717811216614414670848817010505/79228162514264337593543950336 : ℚ), (-75011881929154090617744856851881192755/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (23/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_24 (hDS : dsval.Contains deltaS) : (⟨(129988810540091572116036751163182887039/19807040628566084398385987584 : ℚ), (893691374761313682350263126030851745049/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (6/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_24 (hDS : dsval.Contains deltaS) : (⟨(-2100080534538734042877373295298058641019/39614081257132168796771975168 : ℚ), (-241181597934782995781679340037971192919/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (6/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_25 (hDS : dsval.Contains deltaS) : (⟨(418331281348663097506733581974714967377/79228162514264337593543950336 : ℚ), (1552217719228139832538035500260225787773/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (25/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_25 (hDS : dsval.Contains deltaS) : (⟨(-3665055156892345158902167023249797086711/79228162514264337593543950336 : ℚ), (-193829351418596492310368663592508409823/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (25/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_26 (hDS : dsval.Contains deltaS) : (⟨(10515608991177738092256756128407780357/2475880078570760549798248448 : ℚ), (2695151959511931922203340712279619029053/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (13/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_26 (hDS : dsval.Contains deltaS) : (⟨(-3196277348715382789013013625719155542465/79228162514264337593543950336 : ℚ), (-77874755170396298112874189228205035253/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (13/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_27 (hDS : dsval.Contains deltaS) : (⟨(270622613994994828695448044299492213601/79228162514264337593543950336 : ℚ), (2339150956684824927599457939580143198521/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (27/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_27 (hDS : dsval.Contains deltaS) : (⟨(-2785926872306638299230102518721564358275/79228162514264337593543950336 : ℚ), (-125133936780379766268531197460180029081/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (27/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_28 (hDS : dsval.Contains deltaS) : (⟨(217602926622368739053208215751108333369/79228162514264337593543950336 : ℚ), (2029624001797255822796005413823372035889/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_28 (hDS : dsval.Contains deltaS) : (⟨(-2426983997332101118431315472166585108005/79228162514264337593543950336 : ℚ), (-50262421637216062239916643893911868611/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
#print axioms host_M200_0
#print axioms host_dM_dz_0
#print axioms host_M200_1
#print axioms host_dM_dz_1
#print axioms host_M200_2
#print axioms host_dM_dz_2
#print axioms host_M200_3
#print axioms host_dM_dz_3
#print axioms host_M200_4
#print axioms host_dM_dz_4
#print axioms host_M200_5
#print axioms host_dM_dz_5
#print axioms host_M200_6
#print axioms host_dM_dz_6
#print axioms host_M200_7
#print axioms host_dM_dz_7
#print axioms host_M200_8
#print axioms host_dM_dz_8
#print axioms host_M200_9
#print axioms host_dM_dz_9
#print axioms host_M200_10
#print axioms host_dM_dz_10
#print axioms host_M200_11
#print axioms host_dM_dz_11
#print axioms host_M200_12
#print axioms host_dM_dz_12
#print axioms host_M200_13
#print axioms host_dM_dz_13
#print axioms host_M200_14
#print axioms host_dM_dz_14
#print axioms host_M200_15
#print axioms host_dM_dz_15
#print axioms host_M200_16
#print axioms host_dM_dz_16
#print axioms host_M200_17
#print axioms host_dM_dz_17
#print axioms host_M200_18
#print axioms host_dM_dz_18
#print axioms host_M200_19
#print axioms host_dM_dz_19
#print axioms host_M200_20
#print axioms host_dM_dz_20
#print axioms host_M200_21
#print axioms host_dM_dz_21
#print axioms host_M200_22
#print axioms host_dM_dz_22
#print axioms host_M200_23
#print axioms host_dM_dz_23
#print axioms host_M200_24
#print axioms host_dM_dz_24
#print axioms host_M200_25
#print axioms host_dM_dz_25
#print axioms host_M200_26
#print axioms host_dM_dz_26
#print axioms host_M200_27
#print axioms host_dM_dz_27
#print axioms host_M200_28
#print axioms host_dM_dz_28
