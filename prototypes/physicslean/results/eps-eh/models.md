# SASHIMI-C・Colossus 準拠の宇宙論計算で何を仮定し、何を証明したか

SASHIMI-C の Extended Press–Schechter（EPS）比較図に掲載した10量は、同じ膨張史から出発する。両処方で変わるのは**線形成長の関数と、現在の質量分散の決め方**である。この二つを変更すると、同じ崩壊障壁・EPS・ホスト質量履歴の式に異なる値が入る。

| 比較する量 | SASHIMI-C 準拠 Lean | Colossus 準拠 Lean |
|---|---|---|
| 膨張率 $H(z)$ | 平坦な物質＋Λの Friedmann 解 | 同じ定義・同じ証明区間 |
| 成長因子 $D$、微分 $D'$、障壁 $\delta_c$ | Carroll–Press–Turner（CPT）フィット、その赤方偏移微分、$\delta_c=1.686/D$ | Heath 積分、その赤方偏移微分、同じ障壁係数 |
| 質量ゆらぎ $\sigma$、分散 $S$、微分 $\partial_M S$ | Ludlow の質量フィットとその微分 | Eisenstein–Hu（EH98）の線形伝達関数、球形トップハット積分、$\sigma_8$ 規格化とその微分 |
| 条件付き密度 $F$ | Yang 式 (14) に自処方の障壁差・分散差を代入 | 同じ式に Heath と EH98 から得た差を代入 |
| ホスト $M_{200}$、質量変化率 | SASHIMI-C の形成史・濃度・ビリアル換算 | 同じ経験式の成長・分散を Heath・EH98 に置換 |

**保証の読み方。** フィットの実数値を $y_{\rm fit}$ とすると、Lean が検査する結論は $L\le y_{\rm fit}\le U$ である。フィットの係数や宇宙論パラメータは、採用値を固定して Lean 内で計算する。物理的な真値 $y_{\rm phys}$ に対して別途 $|y_{\rm phys}-y_{\rm fit}|\le e$ が裏付けられれば、

$$
 L-e\ \le\ y_{\rm phys}\ \le\ U+e
$$

まで言える。比較図の帯には、この $e$ や観測パラメータの不確かさは含めていない。多くの処方は実数関数の**定義**として書き、Physlib の解への接続や外部の物理量との一致は必要な定理の**仮定**に明示する。こうして自然界との一致を未証明の独自公理にせずに計算できる。[比較図](matched-comparison.pdf) · [各点の厳密な端点と検査区分](matched-samples.csv)

## 1. 入力、単位、Physlib から導く膨張率

SASHIMI-C 準拠モデルと Colossus 準拠モデルは、同じ膨張履歴のもとで成長・構造形成の処方を比較する。現在の物質密度パラメータを $m=\Omega_{m0}=0.315$、宇宙定数の密度パラメータを $\Omega_{\Lambda0}=0.685$ とし、空間曲率と放射の背景密度をゼロに固定する。無次元 Hubble パラメータは $h=0.674$ であり、$H_0$ を $\mathrm{km\,s^{-1}\,Mpc^{-1}}$ で表すと $100h=67.4$ である。Lean はこれらの小数を有理数として受け取る。赤方偏移は無次元で、$z=i/4$、$i=0,\ldots,28$ の29点を評価する。[共通入力と単位](../eps-plots/background.json) · [標本座標](../eps-plots/grid.json)

Colossus 準拠モデルのパワースペクトルでは、全物質に含まれるバリオンの密度 $\Omega_{b0}=0.049$、初期スペクトルの傾き $n_s=0.965$、振幅 $\sigma_8=0.811$、伝達関数に入る現在の CMB 温度 $T_{\rm CMB}=2.7255\,\mathrm K$ を指定する。膨張と成長の計算では両モデルとも放射項をゼロに保つ（外部 Colossus は `relspecies=False`）。SASHIMI-C 準拠モデルの Ludlow 分散フィットは質量と $h$ から評価し、その係数に振幅と形状が織り込まれている。[Colossus 準拠の入力](../../src/SashimiEisensteinIntegral.lean) · [Ludlow 分散フィット](../../src/SashimiSpectrumFit.lean)

**Friedmann 解との接続。** Physlib の `FriedmannSolution` は、時間区間上で第1・第2 Friedmann 方程式を満たす FLRW 解を表す。$a(t)$ をそのスケール因子、$t_0$ を基準時刻、$t$ を光の放射時刻とし、$H(t)=\dot a(t)/a(t)$ とする。接続定理は両時刻が解の区間に属すること、平坦性、$H(t_0)>0$、$H(t)\geq0$、放射流体が空であること、および次の物質密度のスケーリングを用いる。

$$
\rho_m(t)=\rho_m(t_0)\left(\frac{a(t_0)}{a(t)}\right)^3,
\qquad 1+z=\frac{a(t_0)}{a(t)}.
$$

前者は定理の仮定であり、後者は Physlib の赤方偏移の定義から従う。基準時の $\Omega_{m0}$ は物質密度と臨界密度の比である。Friedmann 方程式の閉包関係から $\Omega_{\Lambda0}=1-m$ が従う。規格化したスケール因子を $a=a(t)/a(t_0)$ と書くと、Physlib の `reducedHubble` は次の $E(a)$ を与える。

$$
E(a)=\sqrt{m a^{-3}+1-m},\qquad H(z)=H_0 E\!\left(\frac{1}{1+z}\right).
$$

数値モデルの領域は $H_0>0$、$0\le m\le1$、$z\ge0$ である。`HzFriedmann.hubbleConstant_eq_modelH` により、基準時の膨張率、密度、赤方偏移が入力と一致する解に、この等式で $H(z)$ の計算区間を移せる。図の実数式は入力だけで評価できる。これを特定の時空の膨張率として読めるのは、その時空が上の条件を満たす場合である。$H_0$ と $H$ は同一の逆時間単位をもつ実数として扱い、比較時に $\mathrm{km\,s^{-1}\,Mpc^{-1}}$ を用いる。両モデルは同じ区間を使い、例えば $z=1$ で $H(1)\simeq120.662943\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$ となる。[Friedmann 解との等式](../../src/HzFriedmann.lean) · [数値区間を解へ移す定理](../../src/SashimiBackground.lean) · [比較に使う共通区間](../../scripts/plot_eh_match.py)

## 2. 線形成長と崩壊障壁

$D(z)$ は線形物質密度ゆらぎの成長因子で、両モデルとも $D(0)=1$ に規格化する。$D'(z)$ は赤方偏移に関する微分である。以下では $0<m\leq1$、$z\geq0$ とし、$u=1+z$、$a=1/u$ を使う。[SASHIMI-C 準拠の成長則](../../src/SashimiGrowth.lean) · [Colossus 準拠の成長則](../../src/SashimiHeath.lean)

**SASHIMI-C 準拠の CPT フィット。** Carroll–Press–Turner の処方は、赤方偏移ごとの物質の割合 $x= m u^3/(1-m+m u^3)$ を使う。平坦な場合のフィットの分母を $P(x)$ と書くと、規格化した成長因子は次の式になる。

$$
P(x)=x^{4/7}-(1-x)+\left(1+\frac{x}{2}\right)\left(1+\frac{1-x}{70}\right),
$$

$$
D_{\rm CPT}(z)=\frac{P(m)}{m}\,\frac{x}{uP(x)}.
$$

$z=0$ で $x=m$ となるため $D_{\rm CPT}(0)=1$。$x'=3x(1-x)/u$ と $P'(x)=dP/dx$ を用いると、実装が計算する微分は次のように表せる。

$$
D'_{\rm CPT}(z)=D_{\rm CPT}(z)
\left[\left(\frac{1}{x}-\frac{P'(x)}{P(x)}\right)x'-\frac{1}{u}\right].
$$

Lean はフィットを微分してこの量を計算する解析式を証明し、指定した各点でその評価を有理数区間で囲む。この区間はフィット式を評価する際の数値誤差を表す。物理的な成長に対するフィットの精度の評価には、摂動計算などとの照合が要る。[規格化・解析微分・区間計算](../../src/SashimiGrowth.lean) · [各赤方偏移の証明区間](../eps-plots/background.json)

**Colossus 準拠の Heath 積分。** 平坦な物質＋宇宙定数の背景で、線形成長の増大モードに Heath の積分表示を採用する。第1節の $E(a)$ を使い、積分の規格化定数は $G(a)/G(1)$ で相殺される。

$$
G(a)=E(a)\int_0^a\frac{ds}{[sE(s)]^3},
\qquad D_{\rm H}(z)=\frac{G(a)}{G(1)}.
$$

$B(a)=m+(1-m)a^3$ と置く。積分の上端 $a=(1+z)^{-1}$ の変化を含めて微分すると、

$$
D'_{\rm H}(z)=\frac{3ma}{2B(a)}D_{\rm H}(z)-\frac{a^2}{B(a)G(1)}.
$$

Lean は $s=t^2$ と変数変換した実装の積分が上の Heath 積分に等しいこと、$D'_{\rm H}$ が規格化した関数の導関数であることを証明する。多項式近似の剰余を抑え、29点すべてで $D$、導関数、障壁の有理数端点をカーネル検査した。線形摂動の増大モードにこの積分を採用することは物理的な前提である。[Heath 積分と導関数の証明](../../src/SashimiHeath.lean) · [29点の具体的証明](growth.lean)

**EPS の崩壊障壁。** 両モデルは球対称崩壊に基づく処方として、現在の線形閾値を $\delta_{c0}=1.686$ に固定し、現在へ外挿した障壁 $\delta_c(z)=\delta_{c0}/D(z)$ を用いる。この定数はモデル入力である。例えば $z=1$ では、同じ $H(1)$ に対して次の値を得る。表は証明区間の中点を丸めたもの。[$z=1$ を含む CPT の端点](../eps-plots/background.json) · [Heath の端点](growth.json)

| $z=1$ のモデル | $D$ | $dD/dz$ | $\delta_c$ |
|---|---:|---:|---:|
| SASHIMI-C 準拠（CPT） | 0.607620 | −0.266632 | 2.774762 |
| Colossus 準拠（Heath） | 0.606754 | −0.265972 | 2.778722 |

## 3. SASHIMI-C 準拠の質量分散

拡張 Press–Schechter（EPS）理論では、質量 $M$ のスケールで平滑化した線形密度コントラストの標準偏差を $\sigma(M,z)$、その分散を $S(M,z)=\sigma(M,z)^2$ と書く。ここでは SASHIMI-C の `prompt_cusps=False` の経路に従い、現在の標準偏差 $\sigma_0(M)=\sigma(M,0)$ を Ludlow ほかのフィットから直接求める。振幅を含む係数はこの式に固定されている。[SASHIMI-C の実装](https://github.com/shinichiroando/sashimi-c/blob/e07ef8fa13503bd80818995c0d366311104fb644/sashimi_c.py#L139-L174)・[Ludlow ほか（2016）](https://arxiv.org/abs/1601.02624)

図の横軸は物理質量 $M$ を太陽質量で表した数値 $\mu=M/M_\odot$ である。フィットの基準質量は $10^{10}h^{-1}M_\odot$ なので、$h=0.674$ を用いて次の無次元量を作る。たとえば $M=10^{12}M_\odot$ なら、$h^{-1}M_\odot$ 単位での数値は $h\mu=6.74\times10^{11}$ となる。

$$
\xi=\frac{10^{10}h^{-1}M_\odot}{M}=\frac{10^{10}}{h\mu},
\qquad B(\xi)=1+1.53\,\xi^{0.275}+3.36\,\xi^{0.198}.
$$

フィットと第2節の成長因子 $D(z)$ を組み合わせると、

$$
\sigma_0(M)=\frac{22.26\,\xi^{0.292}}{B(\xi)},\qquad
\sigma(M,z)=D(z)\sigma_0(M),\qquad
S(M,z)=D(z)^2\sigma_0(M)^2.
$$

SASHIMI-C の `s_func(M)` と今回の質量分散の図は $z=0$ を用いる。このとき $D(0)=1$ である。条件付き質量分率の計算でも、親と子の分散をともに $z=0$ で評価し、時刻の違いを崩壊障壁 $\delta_c(z)=1.686/D(z)$ に入れる。[分散と赤方偏移の実装](../../src/SashimiEpsLudlow.lean)・[条件付き質量分率の実装](../../src/SashimiEpsKernel.lean)

$z$ を固定して物理質量で微分すると、$d\xi/dM=-\xi/M$ より

$$
\frac{\partial S(M,z)}{\partial M}
=-\frac{2S(M,z)}{M}
\left[
0.292-
\frac{0.275(1.53)\,\xi^{0.275}+0.198(3.36)\,\xi^{0.198}}{B(\xi)}
\right].
$$

この微分の単位は $M_\odot^{-1}$ で、$\sigma$ と $S$ は無次元である。$M>0$ では分数項は $0$、$0.275$、$0.198$ の加重平均なので、角括弧は正である。したがって分散は質量の増加とともに小さくなる。$M=10^{12}M_\odot$、$z=0$ での値は $\sigma\simeq2.2139104$、$S\simeq4.9013993$、$\partial S/\partial M\simeq-1.4581942\times10^{-12}\,M_\odot^{-1}$ となる。[解析微分の証明](../../src/SashimiSpectrumFit.lean)・[各点の厳密な有理数端点](../eps-plots/variance.json)

Lean の `variance_derivative` は、この微分が同じ分散関数の導関数であることを $M>0$ で証明する。公開した73質量すべてで、$\sigma$、$S$、$\partial S/\partial M$ の各区間がこれらの実数式の値を含むことをカーネルで検査した。物理的な前提は、採用した成長則と固定係数のフィットを対象質量まで適用できることである。最小質量 $10^{-6}M_\odot$ までの低質量側への外挿は、暗黒物質の自由流によるカットオフを省いたフィットをそのまま用いる。この外挿に伴うモデル誤差を定量化するには、対象とする暗黒物質モデルのパワースペクトルやシミュレーションとの比較が必要である。[数値包含の証明](../../src/SashimiEpsLudlow.lean)・[SASHIMI-C の質量分散の分岐](https://github.com/shinichiroando/sashimi-c/blob/e07ef8fa13503bd80818995c0d366311104fb644/sashimi_c.py#L125-L151)

## 4. Colossus 準拠の線形スペクトルと質量分散

Colossus 1.4.0 の `modelEisenstein98` に合わせ、赤方偏移 $z=0$ の冷たい暗黒物質（CDM）とバリオンを合わせた線形物質パワースペクトルを使う。入力は $h=0.674$、$\Omega_m=0.315$、$\Omega_b=0.049$、$T_{\rm CMB}=2.7255\,\mathrm K$、$n_s=0.965$、$\sigma_8=0.811$。EH98 は断熱初期ゆらぎを仮定し、初期宇宙の光子と質量ゼロのニュートリノ3種の寄与を物質と放射の密度が等しくなる時期の計算に含める。晩期の背景・成長では放射を無視する。

**伝達関数。** $k$ を $h/\mathrm{Mpc}$ 単位の共動波数、$\kappa=hk$ を $\mathrm{Mpc}^{-1}$ 単位の同じ波数とする。$f_b=\Omega_b/\Omega_m$、$\theta=T_{\rm CMB}/(2.7\,\mathrm K)$、$k_{\rm eq}=0.0746\,\Omega_m h^2\theta^{-2}\,\mathrm{Mpc}^{-1}$、$q=\kappa/(13.41k_{\rm eq})$ とおく。CDM とバリオンの**振幅**を加え、二乗して単位 $(\mathrm{Mpc}/h)^3$ のパワーを得る。

$$
T=(1-f_b)T_c+f_b T_b,\qquad P(k)=A k^{n_s}T(k)^2 .
$$

バリオンが光子の抗力から解放される時刻の音響地平線を $s$、$u=\kappa s$、$f=[1+(u/5.4)^4]^{-1}$ とする。$L_\beta=\ln(e+1.8\beta q)$、$C_\alpha=14.2/\alpha+386/(1+69.9q^{1.08})$、$U(\alpha,\beta)=L_\beta/(L_\beta+C_\alpha q^2)$ と書くと、正の波数で両成分は次の形になる。

$$
\begin{aligned}
T_c &= fU(1,\beta_c)+(1-f)U(\alpha_c,\beta_c),\\
T_b &= \left[\frac{U(1,1)}{1+(u/5.2)^2}
       +\frac{\alpha_b e^{-(\kappa/k_{\rm Silk})^{1.4}}}{1+(\beta_b/u)^3}\right]
       \frac{\sin y}{y},\qquad
y=\frac{u}{[1+(\beta_{\rm node}/u)^3]^{1/3}} .
\end{aligned}
$$

$s,k_{\rm Silk},\alpha_c,\beta_c,\alpha_b,\beta_b,\beta_{\rm node}$ は EH98 の式 (2)–(24) から得る。Colossus が CAMB と照合した範囲は約 $100\,h/\mathrm{Mpc}$ までで、小スケールではバリオン圧を無視している。最小質量 $10^{-6}M_\odot$ の代表波数 $1/R\simeq8.16\times10^5\,h/\mathrm{Mpc}$ には、粒子の自由流によるカットオフを加えずにこのフィットを外挿する。[EH98 原論文 §§2–3](https://arxiv.org/html/astro-ph/9709112v1) · [Colossus の仕様](https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html#cosmology.power_spectrum.modelEisenstein98) · [伝達関数の実装](../../src/SashimiEisensteinHu.lean)

**分散と質量微分。** 物理質量 $M$ について、第3節と同じく $\mu=M/M_\odot$ と書く。Colossus に渡す $M_\odot/h$ 単位の数値は $h\mu$ である。$\rho_*=277536627245.708$ を $h^2 M_\odot/\mathrm{Mpc}^3$ 単位の現在の臨界密度とする。半径 $R$ は $\mathrm{Mpc}/h$ 単位で、$x=kR$ は無次元である。

$$
R=\left(\frac{3h\mu}{4\pi\Omega_m\rho_*}\right)^{1/3},\qquad
W(x)=\frac{3(\sin x-x\cos x)}{x^3}\ (x>0),\qquad W(0)=1 .
$$

半径 $R_8=8\,\mathrm{Mpc}/h$ での標準偏差が $\sigma_8$ になるように振幅を決める。

$$
w(k)=k^{2+n_s}T(k)^2,\qquad
J(R)=\int_0^\infty w(k)W(kR)^2\,dk,\qquad A=\frac{2\pi^2\sigma_8^2}{J(8)} .
$$

$$
S(M,0)=\sigma_8^2\frac{J(R(M))}{J(8)},\qquad \sigma(M,0)=\sqrt{S(M,0)} .
$$

$dR/dM=R/(3M)$ と $xW'(x)/3=\sin(x)/x-W(x)$ より、

$$
\frac{dS}{dM}=\frac{\sigma_8^2}{M J(8)}
\int_0^\infty w(k)\,2W(kR)
\left[\frac{\sin(kR)}{kR}-W(kR)\right]dk .
$$

Lean は $M>0$、$0<r_{\min}<R(M)$、$Kr_{\min}\ge1$（$K=10^{10}\,h/\mathrm{Mpc}$）の下で積分と微分の交換を証明し、上式がこの実数モデルの質量微分であることを示す。Colossus は $kR<10^{-3}$ で窓を 1 に置き、既定値にスプラインを使う。補間を切った比較では数値積分と有限差分を使う。Lean の区間は指定したフィットの実数値を囲み、観測やシミュレーションに対するフィット誤差は別途評価を要する。[質量と分散の仕様](../../src/SashimiEhVariance.lean) · [全波数での質量微分](../../src/SashimiEhFullDerivative.lean) · [Colossus の規格化・窓・分散](https://bdiemer.bitbucket.io/colossus/cosmology_cosmology.html#cosmology.cosmology.Cosmology.filterFunction)

## 5. 両処方に共通する EPS とホスト

### 二時刻の条件付き質量分率

Yang et al. (2011) の式 (14) は、親ハローの質量のうち、早い時刻に分散 $S$ から $S+dS$ の範囲の前駆ハローに属した質量分率を $F\,dS$ と表す。親を $(M_0,z_0)=(10^{12}M_\odot,0)$、子を $(m,z_1)=(m,1)$ とし、$10^{-6}M_\odot\le m<M_0$ の72点を評価した。障壁を $B(z)=1.686/D(z)$、$d=B(z_1)-B(z_0)$、$s=S(m,0)-S(M_0,0)$ とすると、

$$
F(d,s)=\frac{d}{\sqrt{2\pi}\,s^{3/2}}\exp\!\left(-\frac{d^2}{2s}\right).
$$

両方の分散は $z=0$ にそろえる。$d$ と $s$ の区間下端がともに正であることを検査し、評価に成功すると、返った有理数区間に $F$ が入ることを証明する。親と同じ質量では $s=0$ となる。EPS の確率過程や現実のハロー集団がこの法則に従うかは数値証明の対象外である。質量変数での密度にはヤコビアン $|dS/dm|$ が必要になる。[EPS の実数仕様と定理](../../src/SashimiEpsKernel.lean)・[分散差の保証](../../src/SashimiEhYang.lean)

### ホストの質量履歴と変化率

$M_{200}(z)$ は平均密度が $200\rho_{\rm crit}(z)$ の球内の質量で、$M_0=M_{200}(0)$ である。$L_0=\log_{10}(M_0/M_\odot)$ として、形成赤方偏移のフィット $z_f=-0.0064L_0^2+0.0237L_0+1.8837$ から $q=4.137/z_f^{0.9476}$ を得る。今回の $z_f=1.2465$、$M_0/q\simeq2.97846\times10^{11}M_\odot$ である。初期赤方偏移 $z_i=0$、$D(0)=1$ の場合、$z_f>0$ と $s_q=S(M_0/q,0)-S(M_0,0)>0$ のもとで、

$$
\begin{aligned}
f&=s_q^{-1/2},\qquad \beta=-f,\qquad \alpha=f\bigl[1+1.686\sqrt{2/\pi}\,D'(0)\bigr],\\
M_{200}(z)&=M_0(1+z)^\alpha e^{\beta z},\qquad
\frac{d\ln M_{200}}{dz}=\frac{\alpha}{1+z}+\beta .
\end{aligned}
$$

図では $z=0,0.25,\ldots,7$ の29点を評価した。$D'(0)$ は各処方の成長関数を赤方偏移で微分した値である。$dM_{200}/dz$ は対数微分に $M_{200}$ を掛けて得る。図の `dM_dz` は SASHIMI-C が降着に用いる処方 $M_{\rm vir}^{\rm fit}(z)\,d\ln M_{200}/dz$ であり、単位は $M_\odot$／赤方偏移である。$M_{\rm vir}^{\rm fit}(z)$ 自体を微分すると濃度・過密度の変化も寄与する。Lean の評価対象は SASHIMI-C が定義するこの質量変化率である。[質量履歴と変化率](../../src/SashimiHalo.lean)・[成長率と分散差の伝播](../../src/SashimiEhHost.lean)

### ビリアル質量への換算

Bryan–Norman の過密度は臨界密度基準で $\Delta=18\pi^2+82y-39y^2$、$y=\Omega_m(z)-1$ である。Correa の濃度 $c=c_{200}$ は $v=1+z$、$L=\log_{10}[M_{200}(z)/M_\odot]$ を用いて $\log_{10}c=a+bL(1+gL^2)$ とする。

| 係数 | $0\le z\le4$ | $z>4$ |
|---|---|---|
| $a$ | $1.7543-0.2766v+0.02039v^2$ | $1.3081-0.1078v+0.00398v^2$ |
| $b$ | $0.2753+0.00351v-0.3038v^{0.0269}$ | $0.0223-0.0944v^{-0.3907}$ |
| $g$ | $-0.01537+0.02102v^{-0.1475}$ | $0$ |

Hu–Kravtsov の換算では $u(c)=\ln(1+c)-c/(1+c)$、$t=\Delta\,u(c)/(200c^3)>0$ を代入する。

$$
\begin{aligned}
p&=-0.4283-0.00313\ln t-0.0000352(\ln t)^2,\\
X&=\bigl[0.5116\,t^{2p}+(3/4)^2\bigr]^{-1/2}+2t,\\
M_{\rm vir}^{\rm fit}&=\frac{\Delta}{200}\frac{M_{200}}{(cX)^3}.
\end{aligned}
$$

ここで保証するのは、選んだ濃度・換算フィットと入力区間から求まる質量・変化率である。厳密な NFW 換算用には、$u(c_{\rm vir})/c_{\rm vir}^3=t$ の正の一意な根と、対応する質量 $M_{200}u(c_{\rm vir})/u(c)$ を囲む別の定理がある。経験式の物理的な較正誤差には独立の評価を要する。[換算フィット](../../src/SashimiHalo.lean)・[NFW 逆問題の定理](../../src/SashimiHaloRoots.lean)

## 6. 積分誤差をどう囲むか

有限帯域を $[L,K]=[10^{-8},10^{10}]\,h/\mathrm{Mpc}$ とする。モデルの $J$ の定義域は $0<k<\infty$ であり、帯域の外は解析的に囲む。有限積分ではスペクトルの重み $w(k)=k^{2+n_s}T(k)^2\ge0$ を利用する。各セル $C$ 上で $g(k)\in[g_-,g_+]$ なら

$$
 g_-\int_C w \ \le\ \int_C wg\ \le\ g_+\int_C w .
$$

これは窓 $W^2$ にも、符号を変えうる微分の窓にも使える。同じ7,680セルで区間を狭める際は、各セル幅 $\Delta k$、中点 $m$ において、微分可能な被積分関数 $v$ の導関数がセル全体で $d_-\le v'\le d_+$ を満たせば

$$
 \left|\int_C v(k)\,dk-\Delta k\,v(m)\right|
 \le \frac{(d_+-d_-)(\Delta k)^2}{8}
$$

となる一般定理を使う。導関数がほぼ一定なら中点の左右で一次の寄与が相殺し、同じ分割でも狭い区間になる。FloatLib の実数式評価と外向き有理数端点を使い、既存の包含とも交差させる。比較した両評価器でセルの両端、112ビット・次数16の設定、半径の区間を揃えた。定理自体は EH98 に固有の係数を要求しない。[中点積分定理](../../src/CertifiedCenteredIntegral.lean) · [重み付き積分定理](../../src/CertifiedWeightedIntegral.lean) · [EH98 への適用](../../src/SashimiEhFixedCells.lean)

帯域外では、採用した EH98 係数から $|T(k)|\le3$（$k>0$）、$|T(k)|\le c_T/k$（$k>K$、$c_T=10^{-9}$）を証明した。$0<r\le R$、$Kr\ge1$ のもとで、$J$ の非負の裾と、微分に使う符号付き積分の裾に、それぞれ次の上限を使う。

$$
 \varepsilon_J=9L^3+\frac{18c_T^2}{r^4 K^2},\qquad
 \varepsilon_{\partial} =36L^3+\frac{84c_T^2}{r^3 K}.
$$

前者は分子と $R=8$ の規格化分母に各々適用し、分母の下端が正であることも検査する。後者は微分式の積分部分に両側の誤差として適用し、その後で $\sigma_8^2/[MJ(8)]$ を掛ける。これにより、帯域外の寄与と有限積分の数値誤差の両方を包含に反映する。[伝達関数の裾の証明](../../src/SashimiEisensteinTail.lean) · [分散の裾](../../src/SashimiEhFullGrid.lean) · [符号付き積分の裾](../../src/SashimiEhSlopeTails.lean)

## 7. 公開した区間はどこまでカーネルで検査されたか

赤方偏移は $0,0.25,\ldots,7$ の29点、物理質量は $10^{-6}$ から $10^{12}\,M_\odot$ の73点（対数で0.25刻み、中間値は公表した有効数字12桁の十進数）。全465点は「物理量と評価座標」の組である。背景4量は $4\times29=116$、質量3量は $3\times73=219$、EPS は72、ホスト2量は $2\times29=58$。

| 対象 | 具体的に成立した検査 |
|---|---|
| SASHIMI-C 準拠の全465点 | 各点で入力から最終区間まで閉じた Lean 証明をカーネル検査。 |
| Colossus 準拠の背景116点 | 同様に閉じた証明。膨張率は共通、残り3量は Heath に基づく。 |
| Colossus 準拠の質量・EPS・ホスト349点、960セル | 全点で有限積分を含め、最終区間まで閉じた証明がある。 |
| 同じ349点、7,680セルの精細区間 | 正しさを証明した評価器を Lean 実行系で走らせ、各有限積分の区間を前提として、裾・規格化・微分・EPS・ホストへの伝播をカーネル検査。全7,680セルの評価を各点で再実行した閉じた証明は未作成。 |
| 7,680セルの追加検査 | $R=8$ の連続128セルは入力から部分積分まで閉じた証明。全349組の保存済み有理数端点について、改良後の区間が従来の7,680セル区間に含まれることを別途カーネル検査。348組で幅が厳密に縮小、$M_{200}(0)=10^{12}M_\odot$ は両方とも幅ゼロ。 |

したがって、**どちらの処方にも465点すべてに対する閉じた保証区間がある。図の最も狭い帯を使う場合、Colossus 準拠の積分については上表の実行系による有限和評価も信頼する。** 端点の入れ子の証明だけでは、内側の区間が実数積分を含むことは導けない。図の区間は各標本点に付く。点をつないだ任意の座標で同じ精度を保証するには、その座標か範囲の評価が要る。[960セルの証拠](full-closed.lean) · [精細区間の含意](full-tight.lean) · [128セルの証明](fixed-octave.lean) · [349組の端点比較](fixed-widths.lean) · [検証ログ](verification.txt)

一般に、評価器の健全性定理は「指定の条件で区間を返したなら真値を含む」という形を取る。どの入力でも希望幅の区間が得られるという主張はしていない。中点積分や重み付き積分の定理は別のスペクトルにも使えるが、今回の具体的な全波数の証明は固定した EH98 係数と裾の評価に依存する。宇宙論パラメータを変えると、係数・裾・規格化分母などの条件と各標本の区間を再検査する必要がある。

**区間と代表値の実例。** $M=10^{12}M_\odot$、$z=0$ の分散 $S$ を比較する。区間の十進表示は外向きに丸めた。

| 処方・方法 | 代表値（近似） | その実数モデルを囲む区間 |
|---|---:|---|
| SASHIMI-C 準拠 | $4.9013993307$ | $[4.9013993307,\ 4.9013993308]$（閉じた証明） |
| Colossus 準拠 | $5.0116001127$ | $[2.55279240,\ 9.73089505]$（960セル、閉じた証明） |
| Colossus 準拠、同じモデル | $5.0116001127$ | $[4.97334849,\ 5.05022745]$（7,680セル、上表の精細区間） |
| 外部 Python Colossus の既定補間 | $5.0120102537$ | 比較用の出力 |

この例の二つの処方の中心値の差は約 $0.1102$ で、同じ EH98 処方の Lean と外部 Colossus の差は約 $0.0004101$ である。全73質量で $S$ の両実装間の最大相対差は既定 Colossus と $0.000247$、同一帯域・同一の規格化手順で $4.70\times10^{-8}$。これらは実装を比較した測定値であり、それだけから計算誤差やモデル誤差の上限は導けない。[厳密な端点・比較数値](matched-samples.csv) · [設定別の比較](README.md)

**図の中心線。** SASHIMI-C 準拠の全量と、Colossus 準拠の背景4量では、基準値は証明区間の有理数中点。Colossus 準拠の $\sigma,S,\partial_M S$ の中心線は対数波数上の合成 Simpson 則・補償和による Lean の binary64 値で、$F$ とホスト2量の中心線はこの分散を下流へ通した値である。精度保証用のセル積分は別計算である。厳密な基準値と図示用の丸め値は CSV で区別できる。$y\in[L,U]$ と表示値 $x$ がわかれば $|x-y|\le\max(|x-L|,|U-x|)$ と言える。EPS のネイティブ計算全体の丸め誤差を全入力について証明したわけではない。例えば EH98 を直接検査した10波数では、Lean と Colossus の binary64 値は10点すべてで一致し、そのうち8点でその丸め値は実数 $T(k)$ の極めて狭い区間の外にある。[伝達関数の具体的な証明](transfer.lean) · [証明端点](transfer.json) · [数値積分の記録](lean.json)

**外部コードの役割。** SASHIMI-C 著者公開版を参照実行し、Colossus は Python 1.4.0 の既定スプライン、直接積分、同一帯域の再積分を比べた。Colossus 系列の $F$ とホストは、Colossus 由来の $S,D,D'$ に上記 SASHIMI-C の式を適用した比較である。Astropy は `FlatLambdaCDM` の $H$ と、その膨張史を使い SciPy で独立に積分した成長・障壁の比較に使う。外部パッケージのソース全体に対するプログラム等価性は証明対象外である。[参照実装と照合方法](../eps-plots/README.md) · [Colossus の実行記録](colossus.json)

**検査環境。** 対象実装 `819022e`、Lean 4.36.0、Physlib `db96e88`、Mathlib `fc4f59e`、FloatLib `d24844f`。ビルド時の公理監査と具体的証明の監査は `propext`、`Classical.choice`、`Quot.sound` のみを許す。主要定理の `#print axioms` でもこの3つを確認した。監査した結論に `sorryAx` や未証明の独自物理公理への依存はない。閉じた証明では Python の生成器が提案した証明書も Lean が検査する。実行系で評価した精細区間には上記の追加の信頼範囲がある。[依存バージョン](../../lake-manifest.json) · [公理監査コード](../../src/HzAudit.lean) · [ビルド・再検査記録](verification.txt)

## 8. 熱史・WDM と今回の比較の境界

今回の10量は冷たい暗黒物質（CDM）を用いる SASHIMI-C の計算である。以下の熱史と温かい暗黒物質（WDM）の仮定は、その10量の包含区間から独立している。

**熱史。** 光子温度 $T$（GeV）におけるエネルギー・エントロピーの有効自由度を $g_i(T)$（$i=\rho,s$）とする。Lean は $0<T\le 10^{16}$ で Saikawa–Shirai の解析フィット $f_i(T)$ を区間評価する。物理的な $g_i$ への適用には、表の隣接温度間で中心値と絶対誤差をそれぞれ $T$ について線形補間し、その結果を $m_i(T),e_i(T)$ とおく。仮定 `SaikawaShiraiDoF` はセル内の全温度で $|g_i(T)-m_i(T)|\le e_i(T)$ が成立することを要求する。このとき $|g_i(T)-f_i(T)|\le e_i(T)+|m_i(T)-f_i(T)|$ を証明できる。[熱史の出典・適用範囲](../../THERMAL.md) · [Lean の仮定と定理](../../src/SashimiThermal.lean)

**WDM。** 同じ $\alpha,k\ge0$、$\nu=1.12$ に対し、二つのパワー比は $R_q(k)=[1+(\alpha k)^{2\nu}]^{-q/\nu}$（$q=5,10$）である。SASHIMI-W の family 固定版の標準実装は $R_5$、公開 main の標準実装は $R_{10}=R_5^2$ を使う。この恒等式から文献の物理的な較正や採用すべき規約の正否までは決まらない。物理パワーの保証には、採用する規約と基準スペクトル $P_{\rm ref}$ を定め、較正帯域内で $|P_{\rm phys}(k)-P_{\rm ref}(k)R_q(k)|\le\epsilon(k)$ を満たすという仮定 `WdmFitAgreement` を要する。フィット値を $[L,U]$ に囲めば、物理パワーは $[L-\epsilon(k),U+\epsilon(k)]$ に入る。[版ごとの規約と一次資料](../../SASHIMI-COVERAGE.md) · [パワー比](../../src/SashimiSpectrum.lean) · [物理パワーへの伝播](../../src/SashimiWdm.lean)

SASHIMI-SI の散乱・コア進化、SASHIMI-F の波動暗黒物質の処方、C の追加カスプ模型にはそれぞれ固有の入力がある。サブハローの個数や消滅収量には集団形成・生存のモデル、観測尤度にはさらに観測データと誤差モデルが必要になる。[系列別の物理量と利用箇所](../../SASHIMI-COVERAGE.md) · [C のカスプ模型](../../CUSP.md)
