# SASHIMI のスカラー評価

`prototypes/physicslean` を作業ディレクトリとして、リポジトリ指定の Lean/Lake と Python 3.10 以上を使う。各呼び出しは必要な Lean モジュールをビルドし、精度と区間幅を指定して評価する。Python 側に物理公式の別実装はない。

```sh
python3 scripts/sashimi.py --list
python3 scripts/sashimi.py c200-correa --describe
python3 scripts/sashimi.py hubble z=2 --replay
python3 scripts/sashimi.py growth-derivative z=2 --replay
python3 scripts/sashimi.py age z=2
python3 scripts/sashimi.py wdm-q5 alpha=1/10 k=1
python3 scripts/sashimi.py nfw-mass r=1/100 rs=1/50 rho_s=1e15 --tol 0.001
python3 scripts/sashimi.py cvir-root c=7 delta=100 lo=1 hi=20 --tol 1e-5
python3 scripts/sashimi.py g-rho-fit T=1 --tol 1e-4
```

Python からは、正確な数の文字列、整数、`Decimal`、`Fraction` を入力できる。バイナリ浮動小数点数を渡した場合はエラーになる。

```python
from scripts.sashimi import evaluate

result = evaluate("growth", {"z": "2"}, tol="1e-5", bits=48, degree=16, replay=True)
# interval = { "lower": "29387730859161/70368744177664",
#              "upper": "14693865631317/35184372088832" }
```

`enclosed` は Lean の実行可能な評価器が返した包含区間で、幅は出力単位で指定した `tol` 以下。一般の健全性定理は各評価器のビルドでカーネルが検査する。`--replay` は個別入力での評価結果と両端点の完全一致も `decide +kernel` で検査して `kernel_replayed: true` を返す。`unresolved` は評価器の定義域、区間幅、根の初期区間・反復数等の条件を満たさず区間が返らなかった状態。CLI の終了値は順に 0、1。引数や Lean 実行にエラーがあった場合は 2。

## 背景モデルの選択

| `--profile` | デフォルトの係数 | 対象 |
|---|---|---|
| `native-flat` | Ωm=0.315, h=0.674, ΩΛ=1−Ωm | C/F/SI と ITAMAE native の平坦背景 |
| `w-consistent` | Ωm=0.27, h=0.70, ΩΛ=1−Ωm | W の整合的な平坦背景 |
| `w-legacy` | Ωm=0.27, h=0.70, Ωb=0.0469, Ωr=0 | W の公開式の係数を使う背景・成長 |
| `flrw` | Ωm=0.315, h=0.674, Ωr=Ωk=0, ΩΛを指定 | 放射・曲率を含む背景量 |

`omega_m`、`h` を上書きできる。背景の `G` は標準値 4.30091×10⁻⁹ Mpc (km/s)²/Msun。`w-legacy` では `omega_b`、`omega_r`、`flrw` では `omega_r`、`omega_k`、`omega_lambda` を明示できる。各組の H スケールは 100h km/s/Mpc。自由係数の和が 1 のときにこれが H(0) と一致する。W の歴史的 `g` は `g(0)=1+Ωb` で、`Hz` は放射項を省き、`Omegaz` は第2係数に二乗を掛ける。SI等への単位変換・SASHIMI固有の重力定数による再規格化は呼出側が指定する。

成長は Carroll–Press–Turner フィットの値を包む。`growth-derivative` はその解析微分、`growth-derivative-f-legacy` は歴史的 F の h⁻² を含む計算式。`growth-secant z1=... z2=...` は F の consistent モードで使うような指定した2点の差分商を評価する。差分商から導関数までの打切り誤差は別に評価する。`growth-ratio` は D(z₂)/D(z₁)。`age` は平坦モデルの z から無限大まで、`lookback` は 0 から z、`elapsed` は z₁ から z₂ の時間積分を Gyr で返す。SI の有限な上限を再現するには `elapsed z1=... z2=1000` を使う。

## スペクトル・ハロー・熱的自由度

- `wdm-alpha` は原式の Mpc/h を返す。`wdm-q5`、`wdm-q10`、`wdm-amplitude`、半モードは α を Mpc、k を Mpc⁻¹ で受ける。Mpc/h の α を h で割るとこの規約の Mpc になる。Ωc は W の全 dark matter が warm であるという仮定での温かい成分の割合。q5 は歴史的なパワー比、q10 は振幅伝達関数を二乗したパワー比。`wdm-half-mode` と `wdm-standard-power-half` はそれぞれ振幅 1/2 と標準パワー 1/2 の点。
- `fdm-ratio` の `m22` は 10⁻²² eV を単位とする。Ludlow の質量では `pivot` と `mass` を同じ Msun 単位で渡す。テーブルからの積分で必要なスペクトル、観測誤差、フィットの物理的な較正誤差は個々のスカラーの演算誤差とは別の入力である。
- NFW の半径は物理 Mpc、密度は Msun/Mpc³、質量は Msun。`nfw-invert-mass` の標的は f(c)、`nfw-invert-mean` は f(c)/c³。`lo` と `hi` は無次元の濃度の探索区間。`mvir-root` は M200、c200、Δを与え、cvir に対する `concentration_tol` と質量に対する `--tol` をともに満たす Mvir 区間を返す。`cvir-root` で濃度を直接得られる。`mvir-fit` は Hu–Kravtsov フィットを包む。
- `c200-w-ludlow` は W の100点グリッドによる濃度。入力 `dS_lo` と `dS_hi` は選択した質量ビンについて規格化したトップハット分散差 S(0.02 M)−S(M) を含む区間で、固定値なら同じ有理数を両端に指定する。
- `g-rho-fit` と `g-s-fit` は T [GeV] における Saikawa–Shirai の解析フィットを包む。原表の誤差列を使った物理的な帯は `scripts/thermal_table.py` が扱う。

## フィットの組合せと単位

| 演算 | モデルの入力と意味 |
|---|---|
| `cold-m200`, `cold-virial-rate` | C の通常モードおよび SI の CDM 背景。z=0 の M₀ [Msun]、0≤zi≤z、Ωm と h を指定する。Ludlow の CDM 分散、規格化 CPT とその解析微分、Correa の c200 と近似 virial 質量変換を組み合わせる。出力はそれぞれ M200(z) [Msun]、dM/dz [Msun/赤方偏移]。 |
| `mass-history`, `mass-rate` | 選択したモデルの D(zi)、dD/dz(zi)、正の ΔS=S(Mᵢ/q,zi)−S(Mᵢ,zi) と dz=z−zi を受ける。FDM や WDM の分散を CDM フィットで置き換えないため、入力を明示する。 |
| `c200-correa` | M200 [Msun] と z を入力し、z≤4 と z>4 の較正式を切り替える。 |
| `wdm-jeans`, `wdm-collapse` | W の熱的粒子質量 [keV]、密度割合、h を用いた Jeans 質量 [Msun] と、Benson の質量依存障壁。通常の球対称崩壊障壁は `collapse-threshold`。 |

`wdm-alpha` の `omega_c`、`mass-history` の ΔS、`c200-w-ludlow` の ΔS は物理モデルの入力である。それらを区間の代表値で置き換えると、証明は選んだ代表値について成立する。W のグリッド評価器は両端点に挟まれる任意の ΔS を扱える。

## ITAMAE からの利用

ITAMAE のソースチェックアウトからは同梱のプロトタイプ、wheel からは同じバージョンのビルド済み Lake プロジェクトを使う。wheel で使う場合は下記の各関数に `project="/absolute/path/to/prototypes/physicslean"` を渡す。どちらも PATH 上の Lake を利用する。

```python
from itamae.cosmology.certified import evaluate_mass, evaluate_scalar, evaluate_table

prediction = evaluate_scalar("hubble", {"z": "2"}, tolerance="1e-6", replay=True)
print(prediction.lower, prediction.upper, prediction.unit, prediction.kernel_replayed)

power = evaluate_table(
    "power.dat", layout="canonical", k="0.1", tolerance="1e-4",
    interpolation="linear", transfer="cdm", placement="continuous",
)
print(power.status, power.width, power.metadata["table"]["sha256"])

# power.dat の k 範囲が [0.01, 1] Mpc⁻¹ を含む場合
slope = evaluate_mass(
    "power.dat", layout="canonical", quantity="top-hat-derivative",
    mass="1e12", density="4e10", k_min=".01", k_max="1", growth=1,
    tolerance="1e-11", interpolation="linear", transfer="cdm", placement="continuous",
)
print(slope.status, slope.unit)  # dS/dM [1/Msun]
```

有限帯域の分散で規格化するには `sigma8="0.8"` を渡す。CAMB・WMAP7 では半径は既定で 8/h Mpc、既に Mpc 表記の canonical 表は `anchor_radius` [Mpc] を明示する。`anchor_k_min` と `anchor_k_max` [Mpc⁻¹]、`anchor_depth` を個別指定できる。両端を省くと表の全域を積分する。半径 0 は sharp-k の有限帯域積分を指定する。入力ターゲットと窓・積分域は `result.metadata["normalization"]` で読める。

`EnclosureResult` の `lower`、`upper`、`width`、`tolerance` は `Fraction`。精度未達時は `status=="unresolved"` で端点と幅は `None` になる。`unit` と `metadata` は物理単位やテーブルのハッシュ等を持つ。`evaluate_scalar` と `evaluate_table` の `replay=True` は個別端点を Lean カーネルで再検査する。テーブルでは元の数値行、列順、単位、規格化を含むリクエストを再実行する。成功時は `kernel_replayed` が真になる。質量依存の評価器は検証済み定理に基づき実行する。数値バックエンドと照合するときは、同じ Ωm、h、単位と規格化を指定する。

`evaluate_mass` の `quantity` は `sharp` / `top-hat` で S、末尾の `-sigma` で √S、`-derivative` で dS/dM を指定する。質量は Msun、平均密度は Msun/Mpc³。`growth` に表の基準時からの成長比を与える方法と、`redshift`、`omega_m`（必要なら `omega_lambda`）で z=0 の表からの CPT 成長を評価する方法がある。標準の dS/dM は D² に比例する。歴史的な F の一次成長係数による微分式を照合する場合は `derivative_growth="legacy-f-linear"` を選ぶ。物理的 WDM を指定するには `wdm_mass_kev` [keV]、`omega_wdm`、`thermal_h` を渡す。既知の α [Mpc] は粒子入力を省いて `parameter` に渡す。この指定は `evaluate_mass` と `evaluate_table` で共通する。質量クエリの規格化は入力表の値に従う。

`scale-factor`、`omega-lambda-flat`、`rho-cdm-flat`、`dt-dz-flat` は平坦背景の派生量で、`dt-dz-flat` は負の dt/dz [Gyr/赤方偏移]。剥離モデルの `stripping-*-csf` と `stripping-*-w` は各系列の係数と年の長さを用いる。入力の E は H(z)/(100h)、Δは臨界密度基準の過密度、`mvir` は Msun である。

## Colossus・CAMB との照合

Colossus 1.4.0、SciPy 1.17.1、NumPy 2.4.6 で `planck18` の Ωm=0.3111、h=0.6766 を ITAMAE の `NativeFlatLCDM` にも与えた。z=2 における通常の倍精度出力を示す。

| モデルと設定 | H [km/s/Mpc] | D(z) | 宇宙年齢 [Gyr] |
|---|---:|---:|---:|
| Colossus、放射あり・補間あり（既定） | 204.058689237 | 0.418238001 | 3.276322072 |
| Colossus、放射あり・補間なし | 204.058689237 | 0.418233672 | 3.276318457 |
| Colossus、放射なし・補間なし | 203.976665769 | 0.418233672 | 3.280463062 |
| NativeFlatLCDM、放射なし・CPT 成長 | 203.976665769 | 0.418753159 | 3.280392916 |

この比較は Colossus の `Cosmology(name="check", **dict(cosmologies["planck18"]), persistence="", relspecies=..., interpolation=...)` の `Hz(2)`、`growthFactor(2)`、`age(2)` と、ITAMAE の `H(2)`、`growth_factor(2)`、`cosmic_time(2)` で再現できる。既定 Colossus から読み出した Ωr=0.00009138647952184145、ΩΛ=0.6888086135204782 を `flrw` に与えると、Lean 区間の両端を有効12桁に丸めた値は 204.058689237 になった。数値比較で同じ精度の丸め値が得られることと、外部コードの演算に対する保証は区別している。

SASHIMI-F の Colossus モードは専用の Planck18 オブジェクトで表を読み、質量を Mpc/h の平滑化半径へ換算する。sharp-k の係数 c は半径を c で割って反映する。F が Colossus 用に書き出すキャッシュは log₁₀ k と log₁₀ P の2列である。Colossus はその表をスプライン補間し、8 Mpc/h のトップハットで cosmology の σ₈ に再規格化する。`sigma(derivative=True)` は d ln σ / d ln R を返す。PhysicsLean の `layout="camb"` には線形値 (k,P) の元表を渡す。Colossus 実行の再現にはバージョン、全宇宙論パラメータ、放射と補間の設定、フィルタ、表・伝達関数・キャッシュの識別子が必要となる。[F のキャッシュ形式](https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L506)・[分散と微分](https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L768)

PhysicsLean のテーブル検査は入力ファイル、単位変換、指定した補間関数、有限の積分域、窓と規格化を固定する。CAMB が計算した物質パワーの物理的精度、表の外側の寄与、Colossus のスプラインやソフトウェア実装の誤差は外部に残る。CAMB を再実行して表を再現するには、生成器のバージョン、設定ファイル（成分密度・ニュートリノ・原始スペクトル・非線形補正を含む）、赤方偏移、出力列の単位、k 格子、外挿方法とファイルのハッシュを記録する。[C の入力表の読込み](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L41)・[W の WMAP7 表](https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L38)
