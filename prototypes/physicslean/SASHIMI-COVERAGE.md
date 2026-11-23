# SASHIMI の宇宙論量と利用箇所

2026年11月24日（日本時間）に、SASHIMI family が組合せとして固定したコード、各 gomeshun リポジトリの公開 `main`、PhysicsLean の分岐点を照合した。対象は Python の計算経路に現れる背景量、線形成長、物質パワー、ハローへの受け渡しである。「直接」はバックエンドあるいはモデルのメソッドで評価する量、「派生」はその戻り値を使って後段が計算する量とする。優先度 P1 は共有の背景と次元換算、P2 は複数モデルの構造形成、P3 は個別モデルの処方を示す。優先度は計算依存順を表す。

## 固定した版

| リポジトリ | family が固定する版 | 当日確認した main | family 側のみ／main 側のみのコミット数 |
|---|---|---|---:|
| SASHIMI family | [1112d1c][family] | 同一 | 0／0 |
| Itamae | [8bd6d13][ipin] | [863cb06][imain] | 57／0 |
| SASHIMI-C | [9334d14][cpin] | [2ffefc8][cmain] | 75／17 |
| SASHIMI-W | [d234ae2][wpin] | [5b845b5][wmain] | 38／18 |
| SASHIMI-SI | [cf1b190][spin] | [8197ca5][smain] | 40／11 |
| SASHIMI-F | [fc4cd83][fpin] | [f4893a3][fmain] | 54／18 |

件数は各履歴で `git rev-list --left-right --count <pin>...<main>` を実行した値である。作業開始時の PhysicsLean ブランチ `12026fa` は Itamae の `main` から4コミット進んでいた。family の Itamae とは `863cb06` を共通祖先として分岐している。family 版にはスペクトル積分や移行用パイプラインが加わる。今回の独立した Lean パッケージを組み込む際は、各モデルの互換 SHA、パラメータ、単位変換に対して検査する必要がある。family の固定版と、さらに開発が進んだ各 `itamae-migration` の先頭も異なる。[対応表][manifest]・[PhysicsLean の開始点][physicsbase]・[family 版の分散計算][ivar]

## バックエンドからモデルまで

式中で $m=\Omega_{m0},\ b=\Omega_{b0},\ \lambda=\Omega_{\Lambda0},\ s=1+z$ と置く。平坦な物質＋Λでは $\lambda=1-m$。C・SI・F の基準値は $(m,b,h)=(0.315,0.049,0.674)$ で、W は WMAP7 の入力表から $(0.27,0.0469,0.7)$、$\sigma_8=0.82$ を読む。[C][cbackground]・[SI][sbackground]・[F][fbackground]・[W の入力表][wdata]

| 優先 | 利用量と式 | 種別、利用箇所、単位 |
|---|---|---|
| P1 | $a/a_0=1/s,\ E^2=m s^3+\lambda,\ E=\sqrt{E^2},\ H=100h E$ | 直接。Itamae の `e2/H`、C・SI・F の `g/Hubble`、W の `g/Hz`。$h$ と $E$ は無次元。Itamae の $H$ は km/s/Mpc、各レガシー計算の $H$ と $H_0$ は内部単位で s⁻¹。$c$ を掛ける距離積分をこれらのハロー計算は要求しない。[Itamae][inative]・[C][cbackground]・[W][wbg] |
| P1 | $\rho_c(z)=3H(z)^2/(8\pi G),\ \rho_m(z)=m\rho_c(0)s^3,\ \Omega_m(z)=ms^3/E^2,\ \Omega_\Lambda(z)=\lambda/E^2$ | 直接および派生。Itamae の `rho_crit/rho_m/omega_m`、各モデルの `rhocrit`。$\rho_m(0)$ は質量と平滑化半径の変換、$\Omega_m(z)$ は成長・ビリアル過密度に使う。C の prompt cusp は冷たい成分 $\rho_{c\mathrm{dm},0}=(m-b)\rho_c(0)$ を使い、半径用の全物質密度とは異なる。物理密度の契約単位は $M_\odot/\mathrm{Mpc}^3$。[バックエンド契約][iproto]・[C の平滑化][cfilter]・[C の cusp][ccusp] |
| P1 | $K(z)=1/[sH(z)]$, $dt/dz=-K$、区間の経過時間 $\int_{z_1}^{z_2}K\,dz$、宇宙年齢 $t(z)=\int_z^\infty K\,dz$、lookback $t(0)-t(z)$ | 直接および派生。Itamae は Gyr への換算定数 $1.0227121650537077\times10^{-3}$ を $H$ に掛け、SciPy `quad` を用いる。SI は観測の基準赤方偏移 $z_r$ を起点とし、$t_L=\int_{z_r}^{z}K$ を補間、$t_U=\int_{z_r}^{1000}K$ を求める。以後 $t_U-t_L(z_f)$ を形成時刻等に使う。SI の積分上端と補間範囲（$z_r$～10）、単位系を保存する。$m=0$ の純Λ模型では有限区間の差は定義できるが、過去無限までの年齢は発散する。[Itamae][inative]・[SI の実装][stime]・[SI main][smainage] |
| P2 | $u=\Omega_m(z), v=\Omega_\Lambda(z),\ \phi=u^{4/7}-v+(1+u/2)(1+v/70),\ g_{\rm CPT}=5u/(2\phi),\ D=g_{\rm CPT}(z)/[s\,g_{\rm CPT}(0)]$ | 直接。Itamae の `growth_factor`、C・SI・F・W の `growthD`。W の `linear_growth_factor` は2要素入力では $D(z_2)/D(z_1)$ を返す。成長抑制関数 $g_{\rm CPT}$ と上のモデル名 `g(z)=E(z)^2` は別の量である。ここでの数値保証は Carroll–Press–Turner の近似式に対するもの。$m=0$ は $D$ のこの正規化の定義域から外れる。[Itamae][inative]・[C][cbackground]・[W の2時刻の比][wratio] |
| P2 | $v'=-3\lambda m s^2/(ms^3+\lambda)^2$, $u'=-v'$, $dD/dz$、$\delta_c(z)=1.686/D(z)$ | 直接および派生。$D'$ は C・W・SI・F の `Mzi/Mzzi/dMdz` が使う。`dOdz` は暗黒エネルギー分率の導関数。$\delta_c$ は Itamae の `collapse_threshold` と CDM の EPS の崩壊しきい値。W は質量依存の Jeans 補正も使う。[C の履歴][chistory]・[F の例外][fderiv]・[W の障壁][wbarrier] |
| P2 | $x=\Omega_m(z)-1,\ \Delta_{\rm vir}=18\pi^2+82x-39x^2$, $r_{200}^3=3M_{200}/(4\pi\,200\rho_c)$、$\rho_{\rm vir}=\Delta_{\rm vir}\rho_c$ | 派生。SASHIMI の `Delc(x)` は **臨界密度** に対するビリアル過密度。$1.686/D$ と区別する。NFW 形状 $f(y)=\log(1+y)-y/(1+y)$、濃度 $c_{200}$、$r_s=r_{200}/c_{200}$、$\rho_s=M_{200}/(4\pi r_s^3 f(c_{200}))$ が $M_{200}\to M_{\rm vir}$ の求根版・近似版に入る。C・SI・F の CDM 濃度は Correa 型の赤方偏移で切替えるフィット。W では抑制された分散と形成時刻から濃度を求める経路もある。[C の変換][cconvert]・[W 固定版の濃度][wconcpin]・[W main の濃度][wconc]・[Itamae の NFW][infw] |
| P2 | $t_{\rm dyn}=(1.628/h)\,\mathrm{Gyr}\,(\Delta_{\rm vir}/178)^{-1/2}/E(z)$、$\Phi=A/(t_{\rm dyn} Hs)$、$dm/dz=\Phi\,m(m/M_{\rm vir})^\zeta$ | 派生。背景を潮汐剥離の赤方偏移 ODE に渡す。$dM/dz$（ホストの質量集積）と $dm/dz$（サブハローの剥離）は役割が異なる。質量・赤方偏移依存の $A,\zeta$、濃度処方、初期値、ODE/Picard/摂動法の選択が最終質量に影響する。[C の係数と ODE][ctide]・[W のソルバ][wstrip]・[F][ftide] |
| P2 | $\Delta^2(k)=k^3P(k)/(2\pi^2)$、$R^3=3M/(4\pi\rho_{m0})$、$S_0=\int\Delta^2 W^2\,d\ln k$、$\sigma=\sqrt{S_0}D$、$S=S_0D^2$、$dS/dM$ | 直接および派生。C・SI は Ludlow 型の質量フィットを使う経路がある。C の prompt cusp、W、F はスペクトルから評価する経路も持つ。top-hat の $W=3(\sin x-x\cos x)/x^3, W(0)=1$ と sharp-$k$ の窓、較正半径、積分域、表の補間、$\sigma_8$ の正規化は入力の一部。$R$ が Mpc なら $k$ は Mpc⁻¹。W の表は $k$ が $h/\mathrm{Mpc}$、$P$ が $(\mathrm{Mpc}/h)^3$、フィルタ質量が $M_\odot/h$。Itamae family の `IntegratedVarianceModel` には有限積分域と窓の指定がある。[family 版の積分][ivar]・[C フィット][cfit]・[W の sharp-k][wfilter]・[F の分散][fvar] |
| P3 | C のカットオフとスペクトルモーメント $\sigma_j^2=\int k^{2j}\Delta^2(k)\,d\ln k$、$R_*$、$\gamma$、$a_{\rm coll}$ | prompt cusp オプションが消費。入力表は $z=31$ の $\Delta^2$ で、$z$ での振幅二乗には $(32a)^{2\times0.901}$ とガウス型遮断を掛ける。$R_*=\sqrt{3}\sigma_1/\sigma_2$、$\gamma=\sigma_1^2/(\sigma_0\sigma_2)$。小スケールの $\delta(a)=\delta a^{0.901}$ と障壁 $1.686 f_{ec}$ は、このオプション固有の近似。ピーク乱数、選別、崩壊確率・残存率が後段の個数や $J$ 因子に入る。[C のモーメント][cmoments]・[選別と規格化][ccusp] |
| P3 | W の $P_{\rm WDM}=P_{\rm CDM}[1+(\alpha k)^{2\nu}]^{-q/\nu}$、$\nu=1.12$、$\alpha\propto m_{\rm WDM}^{-1.11}\Omega_c^{0.11}h^{1.22}$、half-mode、$1+z_{eq}^{\rm fit}=3600\,m h^2/0.15$、$M_J$、$\delta_{\rm WDM}(M,z)$ | 直接。$q=5,10$ は下記の版差に従う。粒子質量の単位は keV、$\alpha$ は Mpc/$h$。$M_J\propto(1+z_{eq}^{\rm fit})^{3/2}(m h^2)^{1/2}g_x^{-1}m_{\rm WDM}^{-4}$。この等価赤方偏移は Jeans 処方のフィットで、放射を省いた Hubble 模型から導出していない。[W の背景と伝達][wfilter]・[half-mode 規約][wmigratepower]・[W の障壁][wbarrier] |
| P3 | F の $k_J=9\sqrt{m_{22}}\,\mathrm{Mpc}^{-1}$、$A=2.22m_{22}^{0.04-\ln(m_{22})/1000}$、$B=0.16m_{22}^{-0.05}$、$x=Ak/k_J$, $T=\sin(x^{2.5})/[x^{2.5}(1+Bx^{3.5})]$、$P_{\rm FDM}=T^2P_{\rm CDM}$ | 直接。$m_{22}$ は粒子質量を $10^{-22}$ eV 単位で表したもの。F はネイティブの top-hat、sharp-$k$ と Colossus の分散・導関数を選べる。$k=0$ の連続延長と、Python が式を割り算で実装する入力範囲を区別する。[F の伝達][ftransfer]・[F の Colossus][fcol] |
| P3 | $M(z),dM/dz$, $S(M_1)-S(M_2)$、$\delta_c(z_1)-\delta_c(z_2)$、$N_a$、濃度散乱、$M_{\rm vir}$、カタログの重み | 派生。これらは背景や線形成長の**利用先**で、最終のモデル予測には EPS カーネル、積分、フィットやサンプリングの仮定が重なる。SI の $(t-t_f)/t_c$ には断面積モデルと内部ハロー密度も必要。背景の包含定理を適用した時点で全カタログや散乱確率の保証が得られるわけではない。[C の EPS 積分][ceps]・[SI の時間利用][stime]・[F の EPS][feps] |

## 実装間の差

- **W の dark-energy 規約**。family 固定版の標準クラスでは $\Omega_c=0.2231,\ \lambda=1-\Omega_c=0.7769$ で、$E(0)^2=1.0469$ になる。標準 `growthD` の零点正規化もこの影響を受ける。family の移行クラスの `consistent` と公開 main の標準クラスでは $\lambda=1-m=0.73$。固定版の `g` は放射の係数を $s^4$ に掛け、`Omegaz(p,z)` の第2係数は $s^2$ に掛かり、`Hz` は物質とΛだけを含む。配布入力で当該係数はゼロ。これを放射つきの一貫した宇宙論として拡張するには、3式を改めて指定する。[固定版の定数][wdataexpr]・[固定版の3式][wbg]・[移行版][wmig]・[公開 main][wmainbg]
- **W のパワーと分散**。固定版の標準は $q=5$ の式そのものをパワー比に用いる。移行クラスでは背景の `physics_mode` と $q$ が独立。公開 main は $q=10$ のパワー比を定義し、振幅はその平方根。半値の条件が伝達振幅かパワー比かも結果に添える。固定版標準の $dS/dM$ は $D$ に比例し、consistent と公開 main の $dS/dM$ は $D^2$ に比例する。[固定版の宣言][wdataexpr]・[移行版の規約][wmigratepower]・[公開 main のパワー][wmainpower]・[公開 main の微分][wmaindiff]
- **F の導関数**。固定版の `dDdz` は中間の $\Omega_\Lambda(z)$ に $m h^{-2}s^3$ を使うが、その $d\Omega_\Lambda/dz$ と `growthD` は $m s^3$ を使う。固定版の移行クラスでは `legacy` がこれを再現し、`consistent` は選んだ成長バックエンドを有限差分する。公開 main の解析式は余分な $h^{-2}$ を取り除いた。F の `xi(M)` は Ludlow 型 $\sigma(M,z)$ に入る無次元の逆質量 $10^{10}M_\odot/(hM)$ で、調査した実行経路で2点相関関数 $\xi(r)$ の計算は見つからなかった。F の固定版ネイティブ分散微分にも $D$ と $D^2$ のモード差がある。解析導関数の保証と有限差分の打切り誤差の保証は別に扱う。[F 固定版][fderiv]・[F の xi と分散][fxi]・[移行クラス][fmig]・[公開 main][fmaingrowth]
- **単位と外部バックエンド**。C・SI・F の既存移行アダプタには $\rho_{c,\rm legacy}(0)/\rho_{c,\rm backend}(0)$ の規格化がある。各実装が採用する $G$、pc と cm の換算、組成、放射・ニュートリノの有無、$h$ を明記して照合する。Colossus の `planck18` は実測時に $h=0.6766,m=0.3111$、相対論的成分ありだった。C と F が CAMB 表を読む際の単位を物理 Mpc に統一すると $k=h k_h$, $P=P_h/h^3$。F の $x=Ak/k_J$ は $k$ と $k_J$ ともに物理 Mpc⁻¹、$k_J=9\sqrt{m_{22}}$。F の固定版は $k=0$ も `np.sinc` で処理する一方、確認した main は分母に $x^{2.5}$ を持つ。F の Colossus 経路は入力密度に $h^2/\mathrm{kpc}^3$ を掛ける。[C の密度][cmig]・[C の CAMB][cfilter]・[F の CAMB][funits]・[F 固定版の伝達][ftransferpin]・[F main の伝達][ftransfer]・[F の密度][fcol]・[SI の密度][smig]
- **Colossus 密度の修正**。作業開始時の Itamae は `rho_c/rho_m` が返す $M_\odot h^2/\mathrm{kpc}^3$ に $10^9$ だけを掛けていた。契約の $M_\odot/\mathrm{Mpc}^3$ にそろえるため、このブランチの両メソッドで $10^9 h^2$ を掛ける。Colossus 1.4.0 の planck18、$z=0$ の $\rho_c$ は修正前 $2.77536627246\times10^{11}$、修正後 $1.27052815397\times10^{11}$。$G=4.30091\times10^{-9}$ を使う独立式は $1.27053030161\times10^{11}$。回帰試験は $z=0,1,3$ の $\rho_c=3H^2/(8\pi G)$ と $\rho_m=\Omega_{m0}\rho_c(0)s^3$ を調べる。旧実装では失敗し、修正後は Itamae の26テストが合格した。[変更前][icoldensity]・[上流移行版の単位説明][imigunits]・[修正実装](../../src/itamae/cosmology/colossus.py)・[試験](../../tests/test_cosmology.py)

**F の分散バックエンドと表の意味**。family 固定版と gomeshun main ともに `sigmaMz` と `dsdm` の既定引数は `filt="tophat", backend="sashimi"`。固定版のネイティブ top-hat は元の $k$ 標本で Simpson 積分し、質量グリッド上の $\log M,\log\sigma$ を一次補間する。範囲外の質量ではグリッドを張り直す。main の `sigma_TopHat` は指定質量の積分を直接呼ぶ。固定版の移行層には、基礎表を補間した後に $T^2(k)$ を連続的に掛ける別構成もある。Colossus を明示した場合は $R$ を Mpc/$h$ で渡し、Colossus 自身の宇宙論、パワーの正規化、数値積分、$d\ln\sigma/d\ln R$ を使う。F の移行アダプタは遮断質量より下で導関数の符号差を報告し、Colossus 利用に明示指定を要求する。PhysicsLean が囲むのは、読込んだ表と指定した補間・転送位置・窓で定義した連続積分である。元の離散 Simpson 値、質量グリッド補間、Colossus 内部の浮動小数点実行を追跡する主張とは区別する。[既定値と分岐][fdispatch]・[固定版の再標本化][fgrid]・[main の直接計算][fmainth]・[移行版の連続転送][fintegrated]・[Colossus の呼び出し][fcol]・[符号差に関する記録][fcolwarning]

## 分散の規格化とハローに渡す量

| 対象 | $\sigma_8$ を適用する場所とフィルタ |
|---|---|
| W family 固定版、標準 | EPS では $k_{\rm cut}=(9\pi/2)^{1/3}\,2.5/R_{\rm TH}$ の sharp-$k$ 積分を100個の質量で計算し、$M_8=4\pi\bar\rho_m(8/h)^3/3$ で補間した $\sigma(M_8)$ を 0.82 に合わせる。濃度では top-hat 分散を別に計算し直し、対応する質量補間値で規格化する。[sharp-$k$][wfilter]・[濃度の top-hat][wconcpin] |
| W family 固定版、移行クラスの consistent | 物理 Mpc で表と $\alpha$ を表し、選択した $q$ の生のパワーを sharp-$k$ で積分した $\sigma(M_8)$ で規格化する。[単位と規格化][wpinanchor] |
| W 公開 main | 同じ補間・伝達による生の $P_{\rm raw}$ に対して、表の**全域**で $V_8=\int k^2P_{\rm raw}W_{\rm TH}(k\,8/h)^2\,dk/(2\pi^2)$ を定義する。$P=P_{\rm raw}\sigma_8^2/V_8$ を EPS と濃度の双方へ渡す。EPS 自体の sharp-$k$ は $k_{\rm cut}=2.5/R_{\rm TH}$ で、旧カットオフの追加係数を使わない。[規格化関数][wmainanchor]・[設計と積分域][wpresc] |
| PhysicsLean | `anchor?` は指定した表、転送、補間、積分域から規格化分散を囲む。`apply?` は分母区間の正値性を検査し、符号を持つ $dS/dM$ を含む任意の分子に $\sigma_8^2/V_8$ を適用する。`spectrum_table.py --sigma8` は対象と同じ転送・補間を強制し、対象の波数・積分域と規格化の積分域を別に受ける。$q$ を変えて毎回規格化する場合、両モデルの $V_8$ も変わるため、生のパワーの点ごとの順序だけから正規化後のすべての波数の順序は決まらない。[規格化の証明](src/SashimiNormalization.lean)・[同一規格化における順序](src/SashimiWdm.lean) |

WMAP7 の全3001行（$h=0.7$, $\sigma_8=0.82$）を使い、$\alpha=0.045$ 物理 Mpc、$k=0.1\,\mathrm{Mpc}^{-1}$ の正規化済みパワーを計算した。幅予算 5000 Mpc³、各積分セルを2分割したときの外向き表示は q5 が $[8925.49,10595.42]$、q10 が $[8925.91,10595.88]$ Mpc³。公開 main の16点 Gauss–Legendre 規格化による q10 参照値 $9953.90503$ Mpc³ は区間に含まれた。q5 も同一 top-hat 規格化を用いた比較実験であり、固定版の質量ビン規格化は別の処方である。両全表結果は一般の包含定理をビルドした Lean の実行結果。個別の全3001行はカーネルで再演していない。独立した正規化帯域を使う3行テストは両端点までカーネルで再演する。[全入力条件・有理数端点](results/sashimi-w-normalization.json)・[回帰試験](tests/test_spectrum_table.py)・[表][wdata]

同じ表ファイルへのパスを `WMAP7_TABLE` に指定して、両方の実行を再現できる。

```sh
lake build +SashimiTableRequest
for q in wdm-q5 wdm-q10; do
  python3 scripts/spectrum_table.py "$WMAP7_TABLE" --format wmap7 --k 0.1 \
    --interpolation linear --placement continuous --transfer "$q" --parameter .045 \
    --sigma8 .82 --depth 0 --anchor-depth 1 --bits 112 --degree 18 --tolerance 5000
done
```

C と SI の `dMdz` は $(\beta+\alpha/(1+z-z_i))M_{\rm vir,fit}$ を返す。非線形変換 $M_{\rm vir,fit}(M_{200}(z),z)$ を合成関数として微分する操作とは区別する。W の濃度は線形質量で最近傍のビンを選び、$f=0.02$ の分散差、100個の試行濃度、有限な残差の選別とゼロの補間を使う。最低質量の外側に及ぶ場合、レガシーは一次スプラインの外挿を含む。また W の濃度では物理質量から表上の数値への換算が固定版の $(M/M_\odot)/h$ から main の $(M/M_\odot) h$ に変わったため、同じ物理質量で選ぶビンも異なり得る。Correa et al. は濃度予測を $10^{-2}$～$10^{16}M_\odot$、$z=0$～20 で提示し、濃度計算に使うシミュレーションによる比較は別に記述している。数値包含区間はこれらの経験的な処方に対する計算誤差を表す。[C の履歴][chistory]・[W 固定版][wconcpin]・[W main の質量換算][wconc]・[Correa et al.](https://doi.org/10.1093/mnras/stv1363)・[ハロー処方](src/SashimiHalo.lean)

## WDM 伝達関数の一次文献の照合

同じ $\alpha>0$、$\nu=1.12$ を固定し、$F(k)=[1+(\alpha k)^{2\nu}]^{-5/\nu}$ と置く。問題は、同じ $F$ を密度揺らぎの振幅比 $T$ に代入するか、パワー比 $R=P_{\rm WDM}/P_{\rm CDM}$ に代入するかである。$P$ は振幅の二乗を平均した量なので、$T$ が振幅比なら $R=T^2$ となる。

| 照合した一次資料 | 紙面に書かれた定義 | 同じスケールでの帰結 |
|---|---|---|
| Viel ほか (2005)、arXiv PDF p.2 式(4)、p.3 式(6),(7) | $T=\sqrt R$、$T=F$。したがって $R=F^2$。係数 0.049 の $\alpha$ は熱的残存粒子向けで、論文は $\nu=1.12$ のフィットの範囲として $k<5\,h/\mathrm{Mpc}$ を記す。 | $F=1/2$ なら $R=1/4$。[定義](https://arxiv.org/pdf/astro-ph/0501562#page=2)・[フィット](https://arxiv.org/pdf/astro-ph/0501562#page=3) |
| Viel・Markovič・Baldi・Weller (2012)、arXiv v2 PDF p.3 式(1) | 2005年の手続きに従うと説明した上で $T_{\rm lin}^2=R=F$ と印刷。同じ $\alpha$ の係数と指数を使う。 | $F=1/2$ なら $R=1/2$。2005年の2式をそのまま合成した結果と食い違う。[式(1)](https://arxiv.org/pdf/1107.4094v2#page=3) |
| Dekker ほか (SASHIMI-W)、arXiv v2 PDF p.3 式(2)、p.4 半モードの定義 | 引用[60]と[61]（上の2論文）を併記し、$T^2=R=F$ と記す。半モードはパワー比 $R=0.5$ で定義する。 | 自論文の2箇所は整合し、2012年の**印刷式**と一致する。2005年の振幅フィットをパワー比に使う場合には指数が2倍必要。[式(2)](https://arxiv.org/pdf/2111.13137v2#page=3)・[半モード](https://arxiv.org/pdf/2111.13137v2#page=4) |
| Schneider ほか (2012)、arXiv v2 PDF p.3 式(8),(9) | 振幅 $T=1/2$ となる波長を定義し、$M_{\rm hm}=(4\pi/3)\bar\rho(\lambda_{\rm hm}/2)^3$ とする。 | $T=F$ の規約では $R=1/4$。Dekker の $R=1/2$ と同じ記号の使用に注意する。[式(8),(9)](https://arxiv.org/pdf/1112.0330v2#page=3) |

Ono ほか (2025) のシミュレーション説明でも arXiv PDF p.3 式(2),(3) は $P_{\rm WDM}=T^2P_{\rm CDM}$ と $T=F$ を併記し、同頁 Table 3 は振幅半値を指定する。$\alpha$ については zoom-in と一様ボックスで異なる処方を使ったと明記されている。同頁の式(4)に印刷された質量指数 $-1.15$ や $g_{\rm WDM}$ 依存は Viel (2005) 式(7)の熱的残存粒子向け第2行と異なる。論文が比較する $\alpha$ の違いと、今回の q5/q10 の二乗の違いは個別に照合する必要がある。[Ono ほかの式(2)～(4)と Table 3](https://arxiv.org/pdf/2505.15226#page=3)

半モードの違いは解析的に検算できる。$k_*=(2^{\nu/5}-1)^{1/(2\nu)}/\alpha$ とすれば $F(k_*)=1/2$ なので、同じ $k_*$ で q5 のパワー比は $1/2$、q10 では $1/4$ である。q10 において**同じパワー半値**を求めるなら $k_{R=1/2}=(2^{\nu/10}-1)^{1/(2\nu)}/\alpha$ とする。$k>0$ では $0<F<1$ であり、同一入力・同一規格化で比そのものを比較すれば q10 のほうが強く抑制する。[Lean の規約と半値の等式](src/SashimiSpectrum.lean)

実装の採用規約を論文から切り離して検査した。family 固定版の標準クラスは $R=F$、同じ固定版の移行クラスは q5 と q10 を選択できる。著者の公開リポジトリ shinichiroando/sashimi-w でも積分の各項に $F$ を1回掛ける。gomeshun の公開 main は $R=F^2$。また Dekker 式(2)直後の記述は $\Omega_{\rm WDM}=0.27$ だが、リポジトリの $\alpha$ に入る $\Omega_c=0.27-0.0469=0.2231$ とは異なる。この差は指数の2倍と独立である。[著者による公開実装][wauthor]・[固定版][wfilter]・[移行版][wmigratepower]・[gomeshun main][wmainpower]・[main の検証記録](https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/docs/w-power-prescription.md)

一次資料で直接確定できるのは掲載式の不整合と上記の代数的帰結である。今回確認した APS の掲載ページ、arXiv の版履歴、書誌検索では、この指数に関する正式な正誤表は見つからなかった。2012年の実シミュレーションの初期条件や Dekker 論文の全図表の計算設定、衛星数に由来する質量下限の変化は、掲載式だけでは確定できない。定量的な制限の更新には、正規化、鋭い $k$ フィルタ、濃度、EPS、選択関数をそろえた再計算が必要になる。[Dekker 掲載ページ](https://doi.org/10.1103/PhysRevD.106.123026)・[Dekker の版履歴](https://arxiv.org/abs/2111.13137)・[Viel ほかの版履歴](https://arxiv.org/abs/1107.4094)

## 熱史との接点

C の prompt cusp は $T_{\rm kd}=0.03$ GeV、$a_{\rm kd}=5.3\times10^{-12}$ を外部定数として使う。W の Jeans 式の $g_x=1.5$ は粒子の有効自由度である。4モデルの調査範囲では、温度に依存する標準模型プラズマの $g_{*\rho}(T),g_{*s}(T)$ を呼び出す計算経路は確認しなかった。PhysicsLean の熱史モジュールでは $\rho_{\rm rad}=\pi^2g_{*\rho}T^4/30$、エントロピー密度 $s_{\rm ent}=2\pi^2g_{*s}T^3/45$ を自然単位で扱う。これらを $a_{\rm kd}$ へ結び付けるには、対象成分の断熱膨張などの仮定が必要になる。[C の定数](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L88-L99)・[W の $g_x$][wbarrier]・[追加した熱史式](src/SashimiThermal.lean)

## 保証の読み方

- 背景・成長は平坦な物質＋Λの物理量と CPT 近似に対して区間を求める。CPT の赤方偏移微分は実際の導関数に接続する。時刻は有限区間と無限上端の年齢それぞれに専用の等式があり、SI の上限 1000 は有限区間として指定できる。W の歴史的な背景・成長規約も選択できる。[プロファイルと実行法](SCALARS.md)・[成長と導関数](src/SashimiGrowth.lean)・[年齢と経過時間](src/SashimiTime.lean)
- ハローについては NFW 変換の一意な根とビリアル質量・誤差幅を検査し、経験的な近似変換と識別する。W 濃度の検証器は入力した分散差で100個の試行から濃度を求める。その分散差をどの表・カットオフ・質量ビンから得たかも結果の条件である。C の任意オプションの数値例は入力したスペクトルからピーク統計を評価する。固定版の元の `powerspectrum31.txt` は配布物から取得できず、同梱した CAMB 検査表の記載帯域で検証した。表の観測誤差、宇宙論フィットの較正誤差、最終カタログの確率分布は、この数値区間の幅には含まれない。[NFW 反転と質量](src/SashimiHaloRoots.lean)・[W 濃度](src/SashimiHaloWarm.lean)・[C の対象量と入力範囲](CUSP.md)
- Physlib に接続された $H(z)$ の区間は、このプロトタイプでカーネル検査された定理である。上表にある外部 Python の出力に同じ定理を使うには、パラメータ、式、単位、丸めを接続する必要がある。とくに、CPT 成長やハロー処方の数値誤差を証明しても、近似式と観測の一致は別の研究課題として残る。[Hubble 基盤](src/HzCertificate.lean)・[カーネル検査の記録](results/verification.txt)
- 範囲評価で各小区間の**全体**を押さえた有限積分に対しては、FloatLib 式の可測性とセルごとの有界性から積分可能性も従う。`integrate_sound_auto` は評価が成功したとき積分の包含と幅を同時に証明する。特異点、区間分割不足、精度ゲートで値を返せない場合は証明結果を採用しない。任意のスペクトル表の観測誤差や無限区間への外挿は、表と積分領域の追加仮定を要する。[有限積分の検証器](src/SashimiNumerics.lean)・[積分可能性と包含](src/SashimiIntegrable.lean)
- 宇宙論的光度距離、角径距離、BAO の音響地平線、観測尤度は、照合した4モデルの上記実行経路の入力には現れなかった。サブハローの空間的な半径や annihilation boost の `luminosity_ratio` はそれぞれハロー構造と光度の比である。[SASHIMI-F の boost 算出][fboost]

[ipin]: https://github.com/gomeshun/itamae/tree/8bd6d135003e91b6ca2717655df250c8f808cd24
[imain]: https://github.com/gomeshun/itamae/tree/863cb0643c61ca4ea7458d7a221f48e5c5bfe515
[cpin]: https://github.com/gomeshun/sashimi-c/tree/9334d14622a0d66082375ea9b66ce21a6d8f1ad2
[cmain]: https://github.com/gomeshun/sashimi-c/tree/2ffefc8f02ff1db2f634aae73c1b2a394914a98e
[spin]: https://github.com/gomeshun/sashimi-si/tree/cf1b19021315d16b29351a77b2306ad646f88265
[smain]: https://github.com/gomeshun/sashimi-si/tree/8197ca598c019a603663f7dc73d3f83fc9d4bfd0
[fpin]: https://github.com/gomeshun/sashimi-f/tree/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d
[fmain]: https://github.com/gomeshun/sashimi-f/tree/f4893a3c24a7b6d3c34983cbc5099148ce770d91
[wpin]: https://github.com/gomeshun/sashimi-w/tree/d234ae250f4343e465be483c5a960d12ea5d3dbf
[wmain]: https://github.com/gomeshun/sashimi-w/tree/5b845b52080b30ab16e29a2de93da7ec347682e0
[family]: https://github.com/gomeshun/sashimi-family/tree/1112d1c422368bde7cfacd7bb94398a3dc8d62de
[manifest]: https://github.com/gomeshun/sashimi-family/blob/1112d1c422368bde7cfacd7bb94398a3dc8d62de/compatibility.toml
[physicsbase]: https://github.com/gomeshun/itamae/tree/12026fa44883f4f8cf13d5672edcfcd7375995f8
[ivar]: https://github.com/gomeshun/itamae/blob/8bd6d135003e91b6ca2717655df250c8f808cd24/src/itamae/variance/integrated.py#L117-L235
[inative]: https://github.com/gomeshun/itamae/blob/863cb0643c61ca4ea7458d7a221f48e5c5bfe515/src/itamae/cosmology/native.py#L36-L94
[iproto]: https://github.com/gomeshun/itamae/blob/863cb0643c61ca4ea7458d7a221f48e5c5bfe515/src/itamae/protocols/backends.py#L7-L46
[infw]: https://github.com/gomeshun/itamae/blob/863cb0643c61ca4ea7458d7a221f48e5c5bfe515/src/itamae/halo/nfw.py#L11-L62
[icoldensity]: https://github.com/gomeshun/itamae/blob/863cb0643c61ca4ea7458d7a221f48e5c5bfe515/src/itamae/cosmology/colossus.py#L44-L53
[imigunits]: https://github.com/gomeshun/itamae/blob/05c3473381de637937c5355bd8f96c01eb8a7c91/src/itamae/cosmology/colossus.py#L135-L185
[cbackground]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L58-L116
[chistory]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L228-L256
[cconvert]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L182-L225
[cfit]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L140-L175
[ctide]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L348-L397
[ceps]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L651-L709
[cfilter]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L30-L86
[ccusp]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L386-L437
[cmoments]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L101-L132
[cmig]: https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c_itamae_migration.py#L197-L245
[sbackground]: https://github.com/gomeshun/sashimi-si/blob/cf1b19021315d16b29351a77b2306ad646f88265/sashimi_si.py#L19-L82
[stime]: https://github.com/gomeshun/sashimi-si/blob/cf1b19021315d16b29351a77b2306ad646f88265/sashimi_si.py#L1353-L1383
[smainage]: https://github.com/gomeshun/sashimi-si/blob/8197ca598c019a603663f7dc73d3f83fc9d4bfd0/sashimi_si.py#L1448-L1479
[smig]: https://github.com/gomeshun/sashimi-si/blob/cf1b19021315d16b29351a77b2306ad646f88265/sashimi_si_itamae_migration.py#L118-L145
[fbackground]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L110-L170
[fderiv]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L152-L169
[fcol]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L755-L822
[fvar]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L650-L748
[fmaingrowth]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L135-L151
[ftransfer]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L495-L529
[ftide]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L927-L962
[feps]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L1274-L1305
[fboost]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L1903-L1929
[fmig]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f_itamae_migration.py#L131-L168
[wdata]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/WMAP7_camb_matterpower_z0_extrapolated.dat#L1-L5
[wdataexpr]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L46-L68
[wbg]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L335-L420
[wratio]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L185-L203
[wfilter]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L75-L130
[wbarrier]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L434-L450
[wmainbg]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w.py#L51-L63
[wmainpower]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w.py#L239-L291
[wmaindiff]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w.py#L1146-L1151
[wstrip]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w.py#L2127-L2165
[wconc]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w.py#L459-L550
[wmig]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w_itamae_migration.py#L205-L258
[wmigratepower]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w_itamae_migration.py#L159-L182
[wconcpin]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w.py#L204-L330
[wpinanchor]: https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/sashimi_w_itamae_variance.py#L170-L230
[wmainanchor]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/sashimi_w_power.py#L42-L113
[wpresc]: https://github.com/gomeshun/sashimi-w/blob/5b845b52080b30ab16e29a2de93da7ec347682e0/docs/w-power-prescription.md#L12-L39
[fxi]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L180-L188
[funits]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L445-L453
[ftransferpin]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L628-L641
[wauthor]: https://github.com/shinichiroando/sashimi-w/blob/99dfc3632eec0126080c0273ddf78f84fe09216c/sashimi_w.py#L85-L108
[fdispatch]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L831-L925
[fgrid]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f.py#L597-L600
[fmainth]: https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L530-L563
[fintegrated]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/sashimi_f_itamae_variance.py#L149-L189
[fcolwarning]: https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/README.md#L128-L137
