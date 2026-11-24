# SASHIMI-C の EPS 計算に使う宇宙論量

[4枚の図をまとめた PDF](eps-comparison.pdf) · [背景・成長・障壁](background.png) · [質量分散・条件付き分布](mass.png) · [ホスト質量](host.png) · [保証精度](precision.png) · [全標本の CSV](samples.csv)

SASHIMI-C の標準 CDM 経路（`prompt_cusps=False`）を選び、赤方偏移を 0 から 7 まで 0.25 刻み、物理質量を $10^{-6}$ から $10^{12}\,M_\odot$ まで 1 decade あたり4点で評価した。中間の質量は有効数字12桁の十進数を厳密な有理数として入力した。下端の質量と上端の赤方偏移はカタログ生成の既定値に合わせた。このホストでのカタログの既定の子質量上限は $10^{11}\,M_\odot$。図はそこから親質量まで延長し、分散差が小さくなる領域も含む。[標本座標](grid.json) · [SASHIMI-C の設定と式](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py)

青い縦線は各標本の実数値を含む区間で、濃紺の線は区間中心を結ぶ。下段の相対差は「比較値 / 区間中心 − 1」。負の導関数で正の相対差を取る点では、比較側の導関数の絶対値が大きい。精度図は $(U-L)/(2|c|)$ と $|x-c|/|c|$ を重ねた（$c=(L+U)/2$、$x$ は SASHIMI-C が返した binary64 値）。これらの差や区間内判定は作図前に有理数で計算した。CSV の `lo_exact` と `hi_exact` は元の端点、`relative_half_width` は表示用の相対半幅。線で結んだ標本間の点は証明の対象外。フィットの物理的な誤差と宇宙論パラメータの不確かさは区間に含めていない。

## 比較した関数

| 量 | ITAMAE が囲む実数式 | 比較値の出典 |
|---|---|---|
| $H(z)$ [km/s/Mpc] | 平坦な物質＋Λ。$H_0=67.4$, $\Omega_{m0}=0.315$, $\Omega_{\Lambda0}=0.685$ | SASHIMI-C の `Hubble` を同一単位に換算。Colossus の `Hz`、Astropy の `FlatLambdaCDM.H` |
| $D,\ dD/dz,\ \delta_c$ | 正規化した Carroll–Press–Turner 成長フィット、その解析導関数、$1.686/D$ | SASHIMI-C の `growthD/dDdz/deltac_func`。Colossus の `growthFactor` とその微分。Astropy の膨張史 $E$ を使って Heath の成長積分 $D\propto E(a)\int_0^a da'/(a'^3 E(a')^3)$ を SciPy で評価し、その解析導関数を計算。各成長関数に同じ障壁係数 1.686 を適用 |
| $\sigma(M,0),\ S=\sigma^2,\ dS/dM$ | $\xi=10^{10}M_\odot/(hM)$ の Ludlow フィットと同じ $S$ の質量微分 | SASHIMI-C の `sigmaMz/dsdm`。Colossus の Eisenstein–Hu パワースペクトル、球形トップハット、$\Omega_{b0}=0.049,\ n_s=0.965,\ \sigma_8=0.811$。物理質量 $M$ を `lagrangianR` に $hM$ として渡し、$dS/dM=2S(d\ln\sigma/d\ln R)/(3M)$ を使用 |
| $M_{200}(z)$, `dMdz` | $M_{200}(0)=10^{12}\,M_\odot$, $z_i=0$ の質量履歴と、$M_{\rm vir,fit}\,d\ln M_{200}/dz$ | SASHIMI-C の `Mzzi/dMdz` を直接評価。ビリアル質量への変換は後者の係数に含まれる |
| 条件付き質量割合密度 $F$（単位 $S$ あたり） | 親 $(M_0,z_0)=(10^{12}\,M_\odot,0)$、子 $(M,1)$ に対する Yang ほかの式 (14) | SASHIMI-C の `Ffunc_Yang`。Colossus の曲線は同ライブラリの $S,D$ を同じ式に代入して算出 |

$F=\Delta\delta\,\exp[-\Delta\delta^2/(2\Delta S)]/(\sqrt{2\pi}\,\Delta S^{3/2})$ において、$\Delta\delta=\delta_c(1)-\delta_c(0)$、$\Delta S=S(M,0)-S(M_0,0)$。Lean は両差の正値性も検査する。$M=M_0$ は $\Delta S=0$ のため条件付き密度から除いた。SASHIMI-C の最終の降着個数 $N_a$ には、さらに正規化、質量制限、ホスト質量の散乱平均、ヤコビアンが入る。[Yang ほか (2011)](https://arxiv.org/abs/1104.1757)

$H$ には Physlib の Friedmann 解の仮定から真の膨張率に区間を移す定理がある。成長、分散、降着率、条件付き分布の証明は上表で選んだモデルの実数式に対する数値保証である。SASHIMI-C の Ludlow フィットの振幅は固定されており、Colossus に指定した $\sigma_8$ は追加のモデル入力である。最小質量の Colossus 半径は $1.2255\times10^{-6}$ Mpc/$h$。分散の補間表の範囲内だが、伝達関数をこの微小スケールまで評価する物理仮定と、粒子の自由流カットオフを含めない選択が効く。Astropy の質量分散は図に用いていない。`colossus.m` に該当する MATLAB 実装は取得したリポジトリに見つからず、比較には Python Colossus 1.4.0 を使用した。[Physlib との接続](../../src/SashimiBackground.lean) · [Ludlow ほか (2016)](https://arxiv.org/abs/1601.02624) · [Colossus](https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html) · [Astropy](https://docs.astropy.org/en/stable/api/astropy.cosmology.FlatLambdaCDM.html)

## 実測結果

- SASHIMI-C の binary64 値と区間中心との差は、全465点で相対 $6.03\times10^{-15}$ 未満だった。154点でその binary64 値自体も区間内にある。$H$ では最大相対全幅が $2.25\times10^{-29}$ に収まり、29点すべてで近傍の倍精度値と厳密な区間を分離できた。[精度図](precision.png) · [全量の集計](summary.json)
- 成長因子の相対差は、SASHIMI-C の CPT フィットに対し Astropy の膨張史から積分した解で最大 0.1663%。導関数の絶対値は $z=0$ で同じ積分が 2.482%、Colossus が 2.491% 大きい。[背景の比較](background.png)
- $10^{-6}\,M_\odot$ で Colossus の分散は 3.373%、質量微分の絶対値は 5.561% 大きく、条件付き密度は 4.522% 小さい。$3.1623\times10^{10}\,M_\odot$ では Colossus の条件付き密度が 4.486% 大きい。この点で分散差だけを Colossus に置換すると +4.176%、障壁差だけを置換すると +0.300% となる。[質量の比較](mass.png)
- $M=10^{10}\,M_\odot$ の条件付き密度 $F$ の証明端点を外向きに丸めると $[0.0118557112081697321,\ 0.0118557112081697562]$。SASHIMI-C は 0.011855711208169748、Colossus の入力を用いた同じ式は 0.012369368424262455 だった。[厳密な端点](kernel.json) · [SASHIMI-C 出力](sashimi.json) · [外部ライブラリ出力と設定](external.json)

## 再現

`prototypes/physicslean` から実行する。作図は Matplotlib 3.11.2、NumPy 2.4.6 と日本語フォント Noto Sans CJK JP を使用した。比較の再計算には SciPy 1.17.1、Astropy 8.0.2、Colossus 1.4.0 を用いた。Lean は 4.36.0、Physlib・Mathlib・FloatLib は Lake のロックファイルに固定している。[依存のコミット](../../lake-manifest.json)

保存済みの数値から図を生成するコマンドは次のとおり。

```sh
python scripts/plot_eps.py results/eps-plots /tmp/sashimi-eps-figures
```

数値と証明を作り直す場合は、`SASHIMI_C_SOURCE` に family 固定版のチェックアウト、または下記著者公開版の `sashimi_c.py` のパスを指定する。アダプタはソースの SHA-256 を照合する。両版のこのグリッドでの参照値は全465個が binary64 で一致した。[著者公開版のソース](https://github.com/shinichiroando/sashimi-c/blob/e07ef8fa13503bd80818995c0d366311104fb644/sashimi_c.py)

```sh
lake build
out=/tmp/sashimi-eps
mkdir -p "$out"
cp results/eps-plots/grid.json "$out/grid.json"
python scripts/eps_plot_background.py --grid "$out/grid.json" --output "$out/background.json"
python scripts/eps_plot_variance.py --grid "$out/grid.json" --output "$out/variance.json"
python scripts/eps_plot_kernel.py "$out/background.json" "$out/variance.json" "$out/kernel.json"
lake env lean "$out/kernel.lean"
python scripts/eps_plot_sashimi.py --source "$SASHIMI_C_SOURCE" --grid "$out/grid.json" --output "$out/sashimi.json"
python scripts/eps_plot_external.py "$out/grid.json" "$out/external.json"
python scripts/plot_eps.py "$out" "$out"
```

背景は96ビット・次数32、分散と EPS カーネルは96ビット・次数40で区間を求めた。個別証明は記録された有理数端点を Lean カーネルで再演する。条件付き密度の証明は入力の障壁と分散から再演できる。各段階の JSON に、依存ライブラリをビルド済みの状態での評価時間と再演時間を秒で記録した。[背景の174区間](background.json) · [分散と導関数の219区間](variance.json) · [条件付き密度の72区間](kernel.json) · [証明書](kernel.lean) · [公理依存](kernel-audit.txt)

| 依存ライブラリ構築後の処理 | 実数区間の計算 | 端点と証明の再検査 |
|---|---:|---:|
| 背景・成長・ホスト履歴（174区間） | 48.32 秒 | 32.67 秒 |
| 分散と質量微分（219区間） | 4.45 秒 | 27.13 秒 |
| 条件付き密度（72区間） | 3.29 秒 | 12.01 秒 |

各段階1回の経過時間で、Lean プロセスの起動を含む。ライブラリの通常ビルドでは 1,753 宣言の公理監査が通った。条件付き密度の保存証明も単独で再検査した。
