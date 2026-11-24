import SashimiEhHost
open SashimiEhHost SashimiEisensteinIntegral FloatLib.Numerics
set_option maxHeartbeats 0
set_option maxRecDepth 100000
private def cr := (⟨(4249493329856275766452151469309587/5192296858534827628530496329220096 : ℚ), (4249493336317676395081585996139171/5192296858534827628530496329220096 : ℚ)⟩ : RationalInterval)
private def pr := (⟨(3181581070066054162317212445488541/2596148429267413814265248164610048 : ℚ), (3181581076224908531977775316882823/2596148429267413814265248164610048 : ℚ)⟩ : RationalInterval)
private def sval := (⟨(-622293623953105633001/1180591620717411303424 : ℚ), (-2489174495812286061699/4722366482869645213696 : ℚ)⟩ : RationalInterval)
private def dsval := (⟨(46520819626288821688751531907475769318132852461283252819531286519425960906379635427021775848724518973733028165768698495514217608922591970205142058106529502275059774003754300935375556565052417463/24255965494136550832163448851745795054858258602564536480525197965139528997249891000428291668414425545409922793299748545130336044921982916090020259799982987293324558992174332959162061100628250000 : ℚ), (424311289150168421200937871341075823114764614088914008656618432127735431151699527222772862219009123978074225518715443719303294098244496915994299323570921681817827279189652713495115388346474828899/211763252438607759058745282782069224482923605673032420803432576298840231124788076840121943296813178520227435726699826731843852969518870821621442413058088119046378051280044019222259521484375000000 : ℚ)⟩ : RationalInterval)
private theorem hc : cr.Contains childRadius := radius_check_sound (bits:=128) (degree:=40) (by decide +kernel)
private def heathRecord : SashimiHeath.Row := ⟨(⟨(1/1 : ℚ), (1/1 : ℚ)⟩ : RationalInterval), sval, (⟨(7961909890118221830291/4722366482869645213696 : ℚ), (1990477472529555457573/1180591620717411303424 : ℚ)⟩ : RationalInterval)⟩
private theorem hSlope : sval.Contains (deriv (SashimiHeath.D matched.om) 0) := by
  have h := (SashimiHeath.check_sound (m:=matched.om) (z:=0) (n:=90) (out:=heathRecord) (by decide +kernel)).2.1
  simpa only [heathRecord, Rat.cast_zero] using h
theorem host_fit_q : (⟨(1142475666052732651622891732858873141851/340282366920938463463374607431768211456 : ℚ), (1142475666052732652448042417655414512359/340282366920938463463374607431768211456 : ℚ)⟩ : RationalInterval).Contains (qExpr.eval (fun _ => 0)) :=
  FloatLib.Numerics.Interval.Expr.contains_gridEnclosure? (show qExpr.gridEnclosure? 128 40 = some (⟨(1142475666052732651622891732858873141851/340282366920938463463374607431768211456 : ℚ), (1142475666052732652448042417655414512359/340282366920938463463374607431768211456 : ℚ)⟩ : RationalInterval) from by decide +kernel)
theorem host_fit_child : (⟨(6334494284973442513925343470422989737487045941513/21267647932558653966460912964485513216 : ℚ), (101351908559575080296006701677892744433518950803535/340282366920938463463374607431768211456 : ℚ)⟩ : RationalInterval).Contains (childExpr.eval (fun _ => 0)) :=
  FloatLib.Numerics.Interval.Expr.contains_gridEnclosure? (show childExpr.gridEnclosure? 128 40 = some (⟨(6334494284973442513925343470422989737487045941513/21267647932558653966460912964485513216 : ℚ), (101351908559575080296006701677892744433518950803535/340282366920938463463374607431768211456 : ℚ)⟩ : RationalInterval) from by decide +kernel)
theorem host_fit_radius : (⟨(278494795067639886701275761539527259951/340282366920938463463374607431768211456 : ℚ), (139247397533819943384192418056976859383/170141183460469231731687303715884105728 : ℚ)⟩ : RationalInterval).Contains (childRadiusExpr.eval (fun _ => 0)) :=
  FloatLib.Numerics.Interval.Expr.contains_gridEnclosure? (show childRadiusExpr.gridEnclosure? 128 40 = some (⟨(278494795067639886701275761539527259951/340282366920938463463374607431768211456 : ℚ), (139247397533819943384192418056976859383/170141183460469231731687303715884105728 : ℚ)⟩ : RationalInterval) from by decide +kernel)
theorem host_M200_0 (hDS : dsval.Contains deltaS) : (⟨(1000000000000/1 : ℚ), (1000000000000/1 : ℚ)⟩ : RationalInterval).Contains (target .m200 (0/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_0 (hDS : dsval.Contains deltaS) : (⟨(-47414779218269823130959427636016956442925/79228162514264337593543950336 : ℚ), (-22781483671690186850176427260128149264369/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (0/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_1 (hDS : dsval.Contains deltaS) : (⟨(34623261379869568118589560236139901791883/39614081257132168796771975168 : ℚ), (69588159919337608242053673903308407556787/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_1 (hDS : dsval.Contains deltaS) : (⟨(-5450377086480977554394124476900577336397/9903520314283042199192993792 : ℚ), (-20864566819421487987787249646311665195061/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_2 (hDS : dsval.Contains deltaS) : (⟨(7502097296822269738035757252883348199171/9903520314283042199192993792 : ℚ), (30299596501918870472087184883947835286043/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_2 (hDS : dsval.Contains deltaS) : (⟨(-9747154025471896144731278831128392649005/19807040628566084398385987584 : ℚ), (-18556660960288396029908234141608364091243/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_3 (hDS : dsval.Contains deltaS) : (⟨(12929237997235228774551753506895617691075/19807040628566084398385987584 : ℚ), (13114998511152159517671595510085889166485/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_3 (hDS : dsval.Contains deltaS) : (⟨(-34404852302229308752422337131293613171697/79228162514264337593543950336 : ℚ), (-32561162070948384615700032235773780078427/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_4 (hDS : dsval.Contains deltaS) : (⟨(44376468010162103889683136372893299841333/79228162514264337593543950336 : ℚ), (88315824714988741252843833880540472919/154742504910672534362390528 : ℚ)⟩ : RationalInterval).Contains (target .m200 (1/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_4 (hDS : dsval.Contains deltaS) : (⟨(-30120714731162894444262364203385035628071/79228162514264337593543950336 : ℚ), (-14169747947221293482378462544790313473423/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (1/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_5 (hDS : dsval.Contains deltaS) : (⟨(37954815361039499893762301988746175159807/79228162514264337593543950336 : ℚ), (38846507992586879590880272282708734247313/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_5 (hDS : dsval.Contains deltaS) : (⟨(-26222055161613592770321967305822596738599/79228162514264337593543950336 : ℚ), (-12263873157307166102168561402283085830943/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_6 (hDS : dsval.Contains deltaS) : (⟨(32379655777562496841632782731386089062887/79228162514264337593543950336 : ℚ), (33286041192616425171141095801337224472327/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_6 (hDS : dsval.Contains deltaS) : (⟨(-11364075427684925643208045999516800010121/39614081257132168796771975168 : ℚ), (-21137961159521192568089494911082896758867/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_7 (hDS : dsval.Contains deltaS) : (⟨(13783214236768367736669916377361909806017/39614081257132168796771975168 : ℚ), (28461343251513554148180236167867272461195/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_7 (hDS : dsval.Contains deltaS) : (⟨(-19629589662602336506485133192511233852595/79228162514264337593543950336 : ℚ), (-4538522544351533643141030516467809895345/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_8 (hDS : dsval.Contains deltaS) : (⟨(23428692103518899848786072918217181084575/79228162514264337593543950336 : ℚ), (24293580148995194976312046465208678653457/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (2/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_8 (hDS : dsval.Contains deltaS) : (⟨(-16903080564357771288124011169418093582299/79228162514264337593543950336 : ℚ), (-15547301196256515477067925490800885103773/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (2/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_9 (hDS : dsval.Contains deltaS) : (⟨(1242720886790499971015664317991479792053/4951760157141521099596496896 : ℚ), (20705793713245531344976635320704800864993/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (9/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_9 (hDS : dsval.Contains deltaS) : (⟨(-7259323615550727831482372545849689503695/39614081257132168796771975168 : ℚ), (-6641633036692620557670607064783156305763/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (9/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_10 (hDS : dsval.Contains deltaS) : (⟨(4213562145268553952351802816893168543127/19807040628566084398385987584 : ℚ), (17625879843175476883911460902676503955617/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_10 (hDS : dsval.Contains deltaS) : (⟨(-1555459646748725555728259551298238723437/9903520314283042199192993792 : ℚ), (-11326118843701861076917442801725050537705/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_11 (hDS : dsval.Contains deltaS) : (⟨(7135733991079740695367176633905744253353/39614081257132168796771975168 : ℚ), (14987978638327404923148555262104416679627/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (11/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_11 (hDS : dsval.Contains deltaS) : (⟨(-5322670944060508125484138676972424386949/39614081257132168796771975168 : ℚ), (-2410179050771633189693883164982077753501/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (11/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_12 (hDS : dsval.Contains deltaS) : (⟨(12073420174070701370233738478751405669771/79228162514264337593543950336 : ℚ), (12732947412469599495397568198282637485667/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (3/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_12 (hDS : dsval.Contains deltaS) : (⟨(-4546020962100137493955673740566621032429/39614081257132168796771975168 : ℚ), (-4096954802053780386441174196198213982837/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (3/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_13 (hDS : dsval.Contains deltaS) : (⟨(5102848107940787963278546172850767914771/39614081257132168796771975168 : ℚ), (5404154492578150934112023268667246975999/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (13/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_13 (hDS : dsval.Contains deltaS) : (⟨(-7754229545522359491329655363878491278941/79228162514264337593543950336 : ℚ), (-6955206485878072637363008046401351762079/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (13/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_14 (hDS : dsval.Contains deltaS) : (⟨(8620760367056006177686570288881789982539/79228162514264337593543950336 : ℚ), (2291977430652685395352117844742655731925/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_14 (hDS : dsval.Contains deltaS) : (⟨(-25800103941856164674253567507541872237/309485009821345068724781056 : ℚ), (-2948523930738948643904459828362896432169/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_15 (hDS : dsval.Contains deltaS) : (⟨(7277338839777788904764185330828074790661/79228162514264337593543950336 : ℚ), (1942857270190131521517958390073759921535/19807040628566084398385987584 : ℚ)⟩ : RationalInterval).Contains (target .m200 (15/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_15 (hDS : dsval.Contains deltaS) : (⟨(-1404845116117607336510155469606845483905/19807040628566084398385987584 : ℚ), (-4994847227290423562728820492713466815819/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (15/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_16 (hDS : dsval.Contains deltaS) : (⟨(3069884457338066023959098837909585818735/39614081257132168796771975168 : ℚ), (6583826485690768835161275775908809367883/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (4/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_16 (hDS : dsval.Contains deltaS) : (⟨(-4776053708367197308350770226487358008249/79228162514264337593543950336 : ℚ), (-4226885408303423674346968934904669305825/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (4/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_17 (hDS : dsval.Contains deltaS) : (⟨(5177356032409043903086149498618133758745/79228162514264337593543950336 : ℚ), (696847124198492143792020965838007856155/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .m200 (17/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_17 (hDS : dsval.Contains deltaS) : (⟨(-506016230671522780311009205819963133285/9903520314283042199192993792 : ℚ), (-3585092795006303077034719592558178126485/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (17/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_18 (hDS : dsval.Contains deltaS) : (⟨(4363765541380425267312736232144635664745/79228162514264337593543950336 : ℚ), (4718125739782011165053224477747845751963/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (9/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_18 (hDS : dsval.Contains deltaS) : (⟨(-3435563690952508088266230981446716046001/79228162514264337593543950336 : ℚ), (-1514387316137597913519393514244210786323/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (9/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_19 (hDS : dsval.Contains deltaS) : (⟨(3676462820419112211509218011650409703559/79228162514264337593543950336 : ℚ), (3991377742605993070835193718598650143075/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (19/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_19 (hDS : dsval.Contains deltaS) : (⟨(-2913880602814919381053051451154574386891/79228162514264337593543950336 : ℚ), (-2557262977108739475941523868665410161971/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (19/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_20 (hDS : dsval.Contains deltaS) : (⟨(3096207507895974146452643506403911162583/79228162514264337593543950336 : ℚ), (3375231198439107681832698820059623862887/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (5/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_20 (hDS : dsval.Contains deltaS) : (⟨(-2470012677955750383078392177197478621065/79228162514264337593543950336 : ℚ), (-1078993352583342066131994795223804749093/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (5/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_21 (hDS : dsval.Contains deltaS) : (⟨(1303301361218816457052914143911483546895/39614081257132168796771975168 : ℚ), (2853157165489975142891681440797858721223/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (21/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_21 (hDS : dsval.Contains deltaS) : (⟨(-1046335355391254008517333018000479337973/39614081257132168796771975168 : ℚ), (-227518817526320466347621053539517793087/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (21/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_22 (hDS : dsval.Contains deltaS) : (⟨(274212163905111532942475748111443771251/9903520314283042199192993792 : ℚ), (1205512671690061795957596851001249099439/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (11/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_22 (hDS : dsval.Contains deltaS) : (⟨(-886063870296625247631682172363442807741/39614081257132168796771975168 : ℚ), (-1534506869254941958836737830464686280625/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (11/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_23 (hDS : dsval.Contains deltaS) : (⟨(922818792352433340716286919103892941313/39614081257132168796771975168 : ℚ), (2036773785820355357621661926439535515573/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (23/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_23 (hDS : dsval.Contains deltaS) : (⟨(-1500021854551648486310162150117490088523/79228162514264337593543950336 : ℚ), (-1293150581672988962373991791757884976679/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (23/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_24 (hDS : dsval.Contains deltaS) : (⟨(1552364225122297063576281369513605191645/79228162514264337593543950336 : ℚ), (1720119364863463755595512485055572894711/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (6/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_24 (hDS : dsval.Contains deltaS) : (⟨(-1269178753537622134665143290276338562541/79228162514264337593543950336 : ℚ), (-1089336559683511300811305937241880639923/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (6/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_25 (hDS : dsval.Contains deltaS) : (⟨(1305349791088320734278643856972099911203/79228162514264337593543950336 : ℚ), (726152665482645776951876974161861462589/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (25/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_25 (hDS : dsval.Contains deltaS) : (⟨(-536726823036112024810003001129369249345/39614081257132168796771975168 : ℚ), (-917318680240788597515658940297166202001/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (25/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_26 (hDS : dsval.Contains deltaS) : (⟨(548686140990212413823173148140562868099/39614081257132168796771975168 : ℚ), (1225882182179737221249182164374381533595/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (13/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_26 (hDS : dsval.Contains deltaS) : (⟨(-907591453263057458827090464259883093937/79228162514264337593543950336 : ℚ), (-96526093960423033628605699358594585665/9903520314283042199192993792 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (13/2 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_27 (hDS : dsval.Contains deltaS) : (⟨(922320450620233248829761346071672868725/79228162514264337593543950336 : ℚ), (1034518135673052090808425134765160150151/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .m200 (27/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_27 (hDS : dsval.Contains deltaS) : (⟨(-767103978098030596941676635876247660219/79228162514264337593543950336 : ℚ), (-649853598382759464826588352189264214069/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (27/4 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_M200_28 (hDS : dsval.Contains deltaS) : (⟨(387513437088597074146700954104904369565/39614081257132168796771975168 : ℚ), (436417849025845342270485929827241279985/39614081257132168796771975168 : ℚ)⟩ : RationalInterval).Contains (target .m200 (7/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
theorem host_dM_dz_28 (hDS : dsval.Contains deltaS) : (⟨(-162040633973640721899127351481272768371/19807040628566084398385987584 : ℚ), (-546728447144783930164933362041108118427/79228162514264337593543950336 : ℚ)⟩ : RationalInterval).Contains (target .virialRate (7/1 : ℚ)) :=
  eval_sound (bits:=96) (degree:=32) hSlope hDS (by decide +kernel)
#print axioms host_fit_q
#print axioms host_fit_child
#print axioms host_fit_radius
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
