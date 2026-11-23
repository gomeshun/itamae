# 標準模型プラズマの有効自由度

Saikawa・Shirai の2018年の計算を採用する。光子温度 $T$ は GeV、$g_{*\rho}$ と $g_{*s}$ は無次元。中心フィットの計算誤差と、状態方程式の不確かさを別々に検査できる。[論文 v2](https://arxiv.org/abs/1803.01038v2)・[付録 B と C、表 1](https://arxiv.org/html/1803.01038)

| $T$ | 中心フィット $g_{*\rho}$ | 表に基づく $g_{*\rho}$ の条件付き区間 | 中心フィット $g_{*s}$ | 表に基づく $g_{*s}$ の条件付き区間 |
|---|---|---|---|---|
| 30 MeV | [12.5681, 12.5684] | [12.3258, 12.8140] | [12.3551, 12.3553] | [12.1410, 12.5672] |
| 1 GeV | [69.7437, 69.7438] | [64.3112, 75.2118] | [68.7715, 68.7716] | [63.4548, 74.0203] |

各列の端点は外側へ丸めた。中心値の証明は付録 C の係数と式を実数式として固定し、FloatLib で包囲する。$0<T<0.12$ では式 C.3–9、$0.12\le T\le 10^{16}$ では式 C.1–2 を用いる。有理数の温度と許容幅を受ける `encloseFit` は定義域外、分母・超越関数の評価失敗、要求幅超過のいずれでも `none` を返す。[Lean 実装](src/SashimiThermal.lean)

## 誤差帯の出所と仮定

`standardmodel2018.dat` は温度、エネルギー自由度の中心と絶対誤差、エントロピー自由度の中心と絶対誤差の5列を持つ。10,000行の温度範囲は $1.9952623\times10^{-6}$ から $9.9738985\times10^{16}$ GeV。表の範囲をはみ出す問い合わせは拒否する。照合には固定版ミラーにある全ファイルを使い、SHA-256 `a4bf5f18dacb084a00e7202a29c689b4ae385950cfeb00be69c98aeb30e8014b` を確認した。組み込み例にはデータ行 1841・1842 と 2511・2512 を用いた。[著者の付随データ](https://arxiv.org/src/1803.01038v2/anc/standardmodel2018.dat)・[検査したミラー](https://raw.githubusercontent.com/carlzimmerman/zimmerman-formula/d7792c9d7033c60da3f93d77a9a60cfec26ab06a/reviews/collision_completion_2026_11_04/production/standardmodel2018.dat)・[取得元とハッシュの記録](https://github.com/carlzimmerman/zimmerman-formula/blob/d7792c9d7033c60da3f93d77a9a60cfec26ab06a/reviews/collision_completion_2026_11_04/production/SOURCE_RECORD.json)

この実装は隣接行の中心 $m$ と絶対誤差 $e$ をそれぞれ **温度 $T$ に対して線形に補間**する。著者の表は温度の対数上でサンプルされている。`SaikawaShiraiDoF cell gρ gs` は、物理的な両関数がセル全域で各補間帯 $[m-e,m+e]$ に含まれるという仮定を表す。表の端点から内点の被覆や特定の信頼水準を導くには追加情報を要する。この仮定のもとで、解析フィット $f$ に対しては $|g-f|\le e+|m-f|$ が成立する。したがって表と解析フィットの中心値がずれてもその差が誤差に入る。

## 物理量への伝播

`density_errors` は $\rho=\pi^2T^4g_{*\rho}/30$ と $s=2\pi^2T^3g_{*s}/45$ の誤差を自然単位系で与える。単位は順に GeV⁴、GeV³。`expansion_error` は $G\ge0$ のとき放射項 $H^2_{\rm rad}=8\pi G\rho/3$ へ伝播する。

共動エントロピー $sa^3$ を保存すると仮定すれば、`adiabatic_scale_cube` が2時刻のスケール因子比を与える。`radiation_transfer` は放射密度の希釈を $[\rho(T)/\rho(T_0)\,(a/a_0)^4]^3$ として導く。1 GeV と 30 MeV の両方の表の誤差を独立に包囲すると、この三乗の値は **[0.091, 0.350]** に入る。許容幅 0.001 の問い合わせは `none` になる。この背景放射の係数と、SASHIMI-F が物質スペクトルに掛ける伝達関数は異なる量である。[SASHIMI-F の物質伝達関数](https://github.com/gomeshun/sashimi-f/blob/f4893a3c24a7b6d3c34983cbc5099148ce770d91/sashimi_f.py#L495-L529)

`Cell.linear_derivative` は選んだセルでの補間関数の傾きを与える。物理的な $g_{*s}'(T)$ を $D$ として別に与えれば、`adiabatic_cooling` が $dT/dt=-3HTg_{*s}/(3g_{*s}+TD)$ を導く。表の点ごとの誤差からこの導関数の範囲は定まらない。SASHIMI-C の prompt-cusp 計算では 30 MeV と $a_{\rm kd}=5.3\times10^{-12}$ を個別に設定し、位相空間密度に $a_{\rm kd}^{-3}$ を使う。熱史からその値を導出するには温度の基準時刻とエントロピー保存の入力が要る。[C の初期値](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L96-L97)・[C の位相空間密度](https://github.com/gomeshun/sashimi-c/blob/9334d14622a0d66082375ea9b66ce21a6d8f1ad2/prompt_cusps.py#L420-L421)

## 再現

`prototypes/physicslean` で実行する。標準のビルドは上の2温度と希釈率、範囲外・許容幅超過を検査する。3番目のコマンドで、取得した全表を使って任意の表内温度の厳密な有理数区間を検査し、条件を明記した Lean 証明を保存できる。

```bash
lake build +SashimiThermal
python3 -m unittest discover -s tests -p test_thermal_table.py -v
python3 scripts/thermal_table.py /path/to/standardmodel2018.dat 0.03 --lean /tmp/ThermalInput.lean
```
