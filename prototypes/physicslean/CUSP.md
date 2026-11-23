# SASHIMI-C の prompt-cusp スペクトルモーメント

SASHIMI-C は CSV `powerspectrum31.txt` の各行を、波数 $k$ (Mpc⁻¹) と赤方偏移31での無次元スペクトル $\Delta^2_{31}$ として読む。各セルで $\log k$ に対して $\log\Delta^2_{31}$ を線形補間する。既定の積分域は $0.1\le k\le10^7$ Mpc⁻¹、$k_{\rm fs}=1.06\times10^6$ Mpc⁻¹。端を越える場合、上流の `interp1d` は両端それぞれの最外セルの傾きで外挿する。この処方を `--edges extrapolate` が再現する。`--edges inside` は表が積分域全体を含むことを要求する。[スペクトルとモーメントを読む実装](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L101-L133)・[カットオフの既定値](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/sashimi_c.py#L124-L130)

$$
\Delta^2(k,z)=\left(\frac{32}{1+z}\right)^{1.802}
\Delta^2_{31}(k)\exp[-(k/k_{\rm fs})^2],\qquad
\sigma_j^2(z)=\int_{k_{\min}}^{k_{\max}} k^{2j}\Delta^2(k,z)\,\frac{dk}{k}.
$$

$\sigma_j$ は $j=0,1,2$ について、それぞれの積分の平方根である。ガウス因子はパワーに1回掛ける。$R_*=\sqrt{3}\,\sigma_1/\sigma_2$ は Mpc、$\gamma=\sigma_1^2/(\sigma_0\sigma_2)$ は無次元。同じ帯域とカットオフなら、この赤方偏移モデルでは $R_*,\gamma$ は赤方偏移に依存しない。`log_measure` が $d\ln k$ と $dk/k$ の変数変換、`redshift_shape` が共通の成長因子の相殺、`peak_sound` が返った有理数区間の包含と許容幅を証明する。[Lean 実装](src/SashimiCusp.lean)

## 入力データと検査結果

family 固定版と公開 main を調べ、著者版 v1.2 の配布物も確認した。`data/` は Git の追跡対象外で、v1.2 の配布はソースアーカイブ2点。元の `powerspectrum31.txt` を取得できなかったため、元表・既定全帯域での数値は未確定である。[著者版の除外設定](https://github.com/shinichiroando/sashimi-c/blob/e07ef8fa13503bd80818995c0d366311104fb644/.gitignore)・[v1.2 配布物](https://github.com/shinichiroando/sashimi-c/releases/expanded_assets/v1.2)

独立の検査用に CAMB 1.6.4 で $z=31$ の線形 CDM 自己パワーを計算し、3点を入力に固定した。以下は $0.1\le k\le0.15$ Mpc⁻¹ における、この表の対数補間モデルに対する結果。区間端を外側へ丸めた。比較値は同じ表・帯域で実際の family 固定版 `prompt_cusps` を実行した浮動小数点の Simpson 積分（1,000点）。$z=0$ では `sigma_function`、$z=31$ では `powerspectrum` を使った。[入力 CSV](tests/data/camb_cdm31.csv)・[CAMB](https://camb.readthedocs.io/)

| 量 | $z=31$ 証明区間 | Simpson | $z=0$ 証明区間 | Simpson |
|---|---|---:|---|---:|
| $\sigma_0$ | [0.0216101, 0.0216856] | 0.02164779 | [0.4906785, 0.4923910] | 0.49153404 |
| $\sigma_1$ (Mpc⁻¹) | [0.00272245, 0.00274008] | 0.00273125 | [0.0618159, 0.0622161] | 0.06201578 |
| $\sigma_2$ (Mpc⁻²) | [0.000351708, 0.000354989] | 0.00035335 | [0.00798588, 0.00806037] | 0.00802308 |
| $R_*$ (Mpc) | [13.28331, 13.49398] | 13.3881900 | [13.28330, 13.49399] | 13.3881900 |
| $\gamma$ | [0.9628003, 0.9878356] | 0.9752363 | [0.9627991, 0.9878368] | 0.9752363 |

Lean は供給された表の有限バンド積分を包囲する。入力 CSV の生成に使うボルツマンコード、CDM の成長モデル、大波数での外挿にはそれぞれ物理仮定がある。ボルツマン計算の精度や、ここから求めるピークの分布・カタログ・消滅光度まで評価するには追加の証明が必要になる。

## 再現

`prototypes/physicslean` で実行する。数値生成の後、Python ランナーは積分と各出力の証明書を `decide +kernel` で検査する。`--tolerance` は選択した各量の単位での絶対幅。要求幅を満たさない量は `null` とし、全体は `unresolved` と終了コード2を返す。`--simpson` の比較には SciPy と NumPy を使用する。依存がビルド済みの環境で、証明生成・検査・比較を合わせて $z=0$ の8量に11.8秒、$z=31$ の5量に11.1秒を要した。

```bash
lake build +SashimiCusp
python3 scripts/cusp_moments.py tests/data/camb_cdm31.csv --edges inside --k-max .15 \
  --redshift 0 --depth 6 --bits 48 --degree 15 --tolerance .25 --simpson \
  --quantities variance0 variance1 variance2 sigma0 sigma1 sigma2 radius gamma \
  --certificate /tmp/CuspInput.lean
```

表の再生成には Python パッケージ `camb==1.6.4` を使う。下記の出力を十進の有理数として読み、以降の補間と積分の入力を確定する。

```python
import camb
import numpy as np

params = camb.CAMBparams()
params.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.1200, mnu=.06, omk=0, tau=.054)
params.InitPower.set_params(As=2.1e-9, ns=.965)
params.set_matter_power(redshifts=[31], kmax=1.0)
results = camb.get_results(params)
power = results.get_matter_power_interpolator(var1="delta_cdm", var2="delta_cdm",
    hubble_units=False, k_hunit=False, nonlinear=False)
for k in (.1, .125, .15):
    delta2 = k**3 * power.P(31, k) / (2 * np.pi**2)
    print(f"{k:.8g},{delta2:.14e}")
```
