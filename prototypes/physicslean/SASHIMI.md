# SASHIMI で使う宇宙論量を Lean で計算する

ITAMAE の物理単位で赤方偏移、粒子質量、ハロー質量やパワースペクトルの表を入力し、指定したモデルの実数値を挟む有理数の区間を返す。絶対幅の要求に達しない計算は `unresolved` を返す。Lean からは `import Sashimi`、Python からは `itamae.cosmology.certified` が入口になる。[演算ごとの引数・単位](SCALARS.md)・[4系列の依存箇所と参照版](SASHIMI-COVERAGE.md)

## どこまで物理理論とつながるか

| 対象 | 証明・入力モデル |
|---|---|
| FLRW 背景 | Physlib の `reducedHubble`、膨張解に対する Friedmann 方程式、赤方偏移とスケール因子を接続した。放射・曲率を含む背景の評価、臨界密度、物質密度が使える。`friedmann_hubble_enclosure` は、2時刻の物質と放射のスケーリング則と基準時の係数・単位の一致から、数値区間を膨張解の $H(t)$ に移す。平坦な物質＋Λの年齢は赤方偏移から無限までの積分に対する証明を持つ。純Λでは有限の時間差を扱う。[背景](src/SashimiBackground.lean)・[時間](src/SashimiTime.lean) |
| 成長と構造形成 | 規格化 CPT 成長率、解析微分、崩壊しきい値、Ludlow 分散フィット、NFW の質量・濃度変換と逆変換、Correa の濃度と質量集積、W の Jeans 障壁を計算する。W のグリッド濃度は選択したビンの分散差の区間を受け取る。較正した経験式のモデル誤差は呼び出し側が与える。[成長](src/SashimiGrowth.lean)・[ハロー](src/SashimiHalo.lean)・[W の濃度](src/SashimiHaloWarm.lean) |
| 線形パワーと質量分散 | 指定した表と補間規則、転送規約、窓、積分帯域に対する積分を評価する。物理単位は k が Mpc⁻¹、P が Mpc³。全セルで関数値を囲い、積分誤差を合計する。赤方偏移ゼロのパワーを前提に成長を掛ける場合、σには D、S には D² が掛かる。[表と積分](src/SashimiTable.lean)・[dk と d ln k の同値性・帯域外の仮定](src/SashimiSpectrumVariance.lean) |
| Prompt cusp と熱史 | C の z=31 入力は無次元の Δ²。指定された表から σ₀、σ₁、σ₂、R*、γ を計算する。熱史には Saikawa–Shirai の中心フィットと表の誤差を採用でき、物理的な真値との関係は定理の仮定 `SaikawaShiraiDoF` に現れる。[C の表と検証範囲](CUSP.md)・[熱史・誤差の伝播](THERMAL.md) |

歴史的な SASHIMI-W の背景は式ごとに異なる密度パラメータを使い、SASHIMI-F の旧成長微分にも現行版との差がある。これらは再現用の名前付きプロファイルで選ぶ。光度距離、観測尤度、サブハローカタログの乱数生成や剥離 ODE 全体の誤差保証は本ライブラリの対象外である。[ソース対応表](SASHIMI-COVERAGE.md)

## WDM の転送規約を指定する

Viel ほか (2005) の式 (4),(6) では $F=[1+(\alpha k)^{2\nu}]^{-5/\nu}$ は振幅であり、パワー比は $F^2$。Viel ほか (2012) の式 (1)、Dekker ほかの式 (2) には同じ指数がパワー比として現れる。ライブラリは `wdm-q10` と `wdm-q5` を明示して呼ぶ。同じ $\alpha$、$\nu=1.12$ なら振幅半値点で前者のパワー比は 1/4、後者は 1/2。この等式は Lean が検査する。物理的なパワーへの適用範囲と誤差は `WdmFitAgreement` の引数に現れる。論文が報告した質量制限への影響について、指数だけから数値を導くことはできない。[一次文献・該当式・後続研究](SASHIMI-COVERAGE.md)・[半値の定理](src/SashimiSpectrum.lean)・[同一規格化での積分の順序と物理的仮定](src/SashimiWdm.lean)

## SASHIMI の入力表での実測

以下は検査対象の表から読み込んだ**同じ参照パワー**と、指定した十進有理数に対する結果。区間の十進端は外側へ丸めた。W の2行は規約差を調べる共通の未規格化の計算条件。ここでの上端 7 Mpc⁻¹ は $h=0.7$ で 10 $h$/Mpc に当たり、Viel (2005) が $\nu=1.12$ について述べた $k<5\,h/\mathrm{Mpc}$ を越える。数値はこの帯域へ指定通り延長した式の積分である。モデル固有の σ₈ 規格化を採用する場合は、同じ転送規則を適用した表からアンカー分散も計算する。[WMAP7 入力](https://github.com/gomeshun/sashimi-w/blob/d234ae250f4343e465be483c5a960d12ea5d3dbf/WMAP7_camb_matterpower_z0_extrapolated.dat)・[Planck 入力](https://github.com/gomeshun/sashimi-f/blob/fc4cd832acf4b4ae4b7ff02691d326e63945dc6d/data/Planck2018_CAMB_extrap.dat)

| データと条件 | 量 | 包含区間 |
|---|---|---|
| WMAP7、α=0.07 Mpc、q5 | k=7×10⁻¹³〜7 Mpc⁻¹、窓なしの有限帯域分散 | [13.4479708, 14.1982982] |
| 同じ表と帯域、q10 | 同じ有限帯域分散 | [11.9066902, 12.5917682] |
| Planck、k=0.674 Mpc⁻¹、CDM | P(k) [Mpc³] | [220.2605593, 220.2605602] |
| 同じ表、FDM、m=10⁻²² eV、節点で補正して対数補間 | P(k) [Mpc³] | [220.1191155, 220.1191221] |

`scripts/verify_sashimi_tables.py WMAP7.dat Planck.dat` で表の SHA-256 を確認し、4件の有理数端点と測定時間を出力できる。時間には入力処理と Lean の起動を含む。長波長・短波長端も含めてトップハットで積分する場合は、安定な小引数評価と振動尾部の上界を備えた `scripts/table_window.py` が使える。[回帰計算](scripts/verify_sashimi_tables.py)・[実測ログ](results/sashimi-tables.json)

質量を入力する `scripts/mass_variance.py` は窓の半径や sharp-k の積分端を有理数に丸めずに評価する。WMAP7 の q5、粒子質量 1.5 keV、ハロー質量 $10^{14} M_\odot$ の例では、0.692〜0.708 Mpc⁻¹ の寄与が $S\in[0.02300039,0.02388015]$、その質量微分が $[-9.609772,-9.608958]\times10^{-15} M_\odot^{-1}$ になった。Planck 表全域、FDM の $m_{22}=1$、$R\in[1,1.001]$ Mpc ではトップハット分散が [7.6474,8.4980] に入る。[質量入力と厳密な端点](results/mass-variance.json)・[smooth/sharp-k の微分定理](src/SashimiMassCalculus.lean)・[W の規格化試験](results/sashimi-w-normalization.json)・[F の全表および部分表の検査記録](results/sashimi-table-window.json)

F の全表について、粗い分割でも実行結果と同じ有理数端点にカーネルが到達することを確認した。3,000 行の単位換算と順序、2,999 セルの積分を188個の部分証明にまとめると、検査は 140.15 秒、区間は [7.24975,9.39316]。各セルを8分割した細かい計算では 23,992 部分区間に対し、証明書付きの計算に 19.95 秒、端点のカーネル再検査に 296.40 秒を要した。前段の [7.6474,8.4980] は、この再検査にも成功した区間を外側へ丸めた値である。元の F 表はアクセス権のあるチェックアウトから同じ SHA-256 を照合できる。公開 CI では、同梱の独立した CAMB 1.6.4 の3,000行でも粗い分割の全表検査を行い、厳密な端点・入力ハッシュ・所要時間を `fdm-kernel.json` に保存する。[全表再検査の実測](results/sashimi-table-window.json)

SASHIMI-W の sharp-k 質量較正も別々に指定する。$R$ を同じ物理的トップハット半径とすると、family 固定版は $k_c=(5/2)(9\pi/2)^{1/3}/R$、公開 main は $k_c=(5/2)/R$ になる。[名前付き較正と出典](SCALARS.md)

同一の Planck 表、$m_{22}=1$、$M=5\times10^{11} M_\odot$、$k=0.690$〜$0.706$ Mpc⁻¹ の例を $z=2$ へ進めると、規格化 CPT 成長と積分の両方の誤差を含めて $S\in[0.0109226,0.0118522]$、$\sigma\in[0.1045113,0.1088677]$ になる。`--derivative-growth legacy-f-linear` は旧 F が返す単一の $D$ 倍を選ぶ。同じ $S=D^2S_0$ の数学的な質量微分には $D^2$ を用いる。[赤方偏移の合成定理](src/SashimiGrowthVariance.lean)・[実表入力の個別証明](tests/SashimiGrowthVarianceCases.lean)

## ビルドと信頼境界

```sh
cd prototypes/physicslean
python3 scripts/fetch_mathlib_cache.py HzAudit
./run_demo.sh
python3 scripts/sashimi.py hubble z=2 --replay
python3 scripts/sashimi.py growth-derivative z=2 --replay
```

スカラーと `scripts/spectrum_table.py` に `--replay` を渡すと、指定した値から出力の有理数端点を Lean カーネルで再計算する。表ではファイル順の有理数列、$h$ による単位換算、補間、転送、有限帯域とアンカーを検査に含める。成功時は `kernel_replayed` が `true` になる。計算資源の上限を超えた場合や端点が一致しなかった場合にはコマンドがエラーを返し、その入力に対する個別証明は得られない。

大きな計算では、次の環境変数に正の整数を設定して再検査の上限を調整できる。

| 環境変数 | 既定値 | 意味 |
|---|---:|---|
| `PHYSICSLEAN_REPLAY_TIMEOUT` | 120 | 個別の Lean プロセスを待つ秒数 |
| `PHYSICSLEAN_REPLAY_MEMORY_MB` | 4096 | Lean のメモリ上限（MB） |
| `PHYSICSLEAN_REPLAY_HEARTBEATS` | 8000000 | 個別定理の Lean heartbeat 上限 |

- `lake build` は `Sashimi` の全モジュールと公理監査を実行する。許可する公理は `propext`、`Classical.choice`、`Quot.sound`。有効自由度や外部表の物理的な妥当性は、適用先の仮定・データとして読む。[公理監査](src/HzAudit.lean)・[依存の固定版](lake-manifest.json)
- 各計算の包含定理は Lean カーネルで検査する。通常の実行で返す有理数には Lean コンパイラとランタイムによる実行経路が入る。個別結果のカーネル再計算を行った記録は `kernel_replayed` または各ドライバの検査ログで確認できる。十進表示の再丸めには、有理数端点から外側への丸めを使う。[FloatLib の信頼境界](https://github.com/lean-dojo/FloatLib/blob/d24844f039320a024c7548f992920562d9a71a60/TRUST.md)
- CAMB などが作った表を真値の近似と解釈する場合、生成時の宇宙論・モデル・離散化誤差が加わる。有限表を越える分散は上下の尾部の仮定を要する。Colossus の補間や経験式の物理精度を含めるには、その処方と誤差の入力が別途必要になる。[外部バックエンドとの比較](SCALARS.md)・[尾部の定理](src/SashimiSpectrumVariance.lean)

[総合ビルド、テスト、主要定理の検査ログ](results/sashimi-verification.txt)

公開 CI が使う WMAP7 表、独立した CAMB 表の生成設定と出典は[回帰テスト用データ](tests/data/README.md)に記録した。
