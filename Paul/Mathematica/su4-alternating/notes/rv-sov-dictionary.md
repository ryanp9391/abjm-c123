# RV-SoV formula dictionary for the alternating chain A F A F (2026-10-12)

Sources (TeX read locally): **RV1** = arXiv:1810.10996, **RV2** = 2002.12341, **GLMRV** = 1910.13442, **GLMR** = 2011.08229,
**BBLR** = 2603.26617. In parentheses: TeX label, and (RV1/RV2/BBLR) the equation number counted from the source (+-1, except RV2
(2.25), (4.12), (1.2), which match the README). Tags: **[Q]** quoted; **[D]** derived here from quoted formulas; **[I]** inference, untested.

## 1. Conventions and dictionary

**RV2 [Q].** Lax $L^\nu(u)=u-\hbar P^\nu$, $P^\nu=\sum_{ij}E_{ij}\otimes\pi^\nu(E_{ji})$ (below `YangR` 2.5); $T(u)=L_L(u-\theta_L)\cdots L_1(u-\theta_1)$ (`spinchain` 2.7), **site $L$ leftmost**;
twist on the right, $\mathbb T=TG$; $f^{[2n]}(u)=f(u+n\hbar)$; quantum minors $T[{}^I_J](u)=\sum_\sigma{\rm sgn}\sigma\,T_{i_{\sigma(1)}j_1}(u)T_{i_{\sigma(2)}j_2}(u-\hbar)\cdots$ (`quantumminor` 2.15: rows antisymmetrised, column $k$ at $u-(k-1)\hbar$);
$\mathbb T_{a,1}=\sum_{|I|=a}\mathbb T[{}^I_I]$; $\det(1-\mathbb Te^{-\hbar\partial_u})=\sum_a(-1)^a\mathbb T_{a,1}e^{-a\hbar\partial_u}$ (`talalaev` 2.13); CBR (2.16); genericity $\theta_\alpha-\theta_\beta\notin\hbar\mathbb Z$ (`thetageneric` 2.9);
companion twist $G_{ij}=(-1)^{j-1}\chi_j\delta_{i1}+\delta_{i,j+1}$, $\chi_j=e_j(z)$ (2.19) $=G_c(z)$ of the README; $\nu_j(u)=\prod_\alpha(u-\theta_\alpha-\hbar\nu_j^\alpha)$; GT vacuum $\langle0|$ (all $\mu^\alpha_{kj}=\nu^\alpha_{k+1}$): $\langle0|T_{j1}(u)=\delta_{j1}\nu_n(u)\langle0|$ (`lowering` 3.19).

| project | RV [D; agrees with tex §2] |
|---|---|
| $\mathcal L_F=u-\theta-hP$, $\rho_F(E_{ij})=-e_{ji}$ (Paul: "fundamental", built with $t_1^n$) | $\pi(E_{ij})=e_{ij}$, $\nu=(1,0,0,0)$, **$\theta^{\rm RV}=\theta_F$** |
| $\mathcal L_A=u-\theta+hK$, $\rho_A(E_{ij})=e_{ij}$ (Paul: "anti-fundamental", built with $t_a$) | $\pi(E_{ij})=-e_{ji}\cong(0,0,0,-1)\sim\nu=(1,1,1,0)$; $\pi^{(1,1,1,0)}(E_{ij})=\delta_{ij}-e_{ji}$ gives $L^{\rm RV}(u-\theta^{\rm RV})=u-\theta^{\rm RV}-\hbar+\hbar K$, so **$\theta^{\rm RV}_A=\theta_A-h$** |
| $T=\mathcal L_1\cdots\mathcal L_{2L}$, $T^G=TG$, $G=G_c(z)$ | identical operator: **our site $s$ = RV site $2L+1-s$**; $\hbar=h$ |
| $t_a(u)$, $a=0..4$ ($t_4$ = qdet) | $\mathbb T_{a,1}(u)$ |
| $V_1=\nu_4$, $V_2=V_3=\nu_2=\nu_3$, $V_4=\nu_1=\prod_A(u-\theta)\prod_F(u-\theta-h)$ [D] | all RV $\nu^\alpha_1=1$; $\nu_4(u)=\prod_A(u-\theta_A+h)\prod_F(u-\theta_F)$ |
| $\langle0|=e_1^T\otimes e_4^T\otimes\cdots$ | RV GT vacuum (unique solution, tex §12) |
| ABJM point $\theta_A-\theta_F=2h$ | $\theta^{\rm RV}_A-\theta^{\rm RV}_F=h\in\hbar\mathbb Z$: violates (2.9) |

**Other papers [Q].** RV1: $R^\nu=u-\hbar P^\nu$, $T=R_L\cdots R_1$, one rectangular rep $(S^A)$ on all sites, same $G$ (`twsitG` 2.24). GLMR/GLMRV/BBLR: $\hbar=i$, $f^{[a]}=f(u+ia/2)$ (**centred** Wronskians: $Q^{\rm BBLR}_I(u)\leftrightarrow Q^{\rm RV}_I(u+\hbar(|I|-1)/2)$
up to gauge); GLMR/GLMRV minors have ascending shifts ($\hbar\to-i$ relative to RV). BBLR: $\mathcal L_{jk}=u\delta_{jk}-i\mathbb E_{kj}$, site 1 leftmost like us, chain of $\mathbf 6=(1,1,0,0)$ ($A=2,S=1$), $Q_j=z_j^{-iu}q_j$;
on a site in the $\mathbf 4$ ($\mathbb E=e$) its Lax is $u-iP$, i.e. our F site with $h=i$. Re-derive shifts before importing any formula.

## 2. SoV basis for fundamental and anti-fundamental sites

**General [Q] (RV2).** $\mu_{kj}=\lambda_{n-k+j-1,j}$ (`mupattern` 3.13), $\bar\mu_{kj}=\mu_{kj}-\nu_{k+1}$ (3.15), $x^\alpha_{kj}=\theta_\alpha+\hbar(\mu^\alpha_{kj}-j+1)$, $B(u)\langle\Lambda^{\bf B}|=\prod_{\alpha,k}\prod_{j\le k}(u-x^\alpha_{kj})\langle\Lambda^{\bf B}|$ (4.3),
$$\langle\Lambda^{\bf B}|=\langle0|\prod_{\alpha}\prod_{k=1}^{n-1}\phi^{\,n-k-1}\big(\mathbb T_{\bar\mu^\alpha_k}(\theta_\alpha+\hbar\nu^\alpha_{k+1})\big)\ (4.12)\;=\;\langle0|\prod_{\alpha,k}\mathbb T^{(k)}_{\bar\mu^\alpha_k}(\theta_\alpha+\hbar\nu^\alpha_{k+1})\ (5.18),\quad
\mathbb T^{(k)}_\xi=\frac{\det_{i,j\le k}Q_i^{[2\hat\xi_j]}}{Q_{12..k}},\ \hat\xi_j=\xi_j-j+1\ (5.12)$$
($\phi:T_{ij}\mapsto T_{i+1,j+1}$ embedding morphism, 3.16-3.18; $\mathbb T^{(k)}$ Bäcklund transform; (1.2) is (5.14): $\mathbb T_{F_k+\bar\mu_k}/\mathbb T_{F_k}$ at $\theta_\alpha+\hbar\nu_n$, read with $\mu_k=\bar\mu_k+\nu_{k+1}$).

**Subtlety [D from (5.1), (5.14), (5.15)].** The embedded $\phi^{n-k-1}(\mathbb T)$ is needed in RV2's *proof* and equals the Bethe-algebra element $\mathbb T^{(k)}$ on $V_{(k)}$. For both our site types $\nu_4=0$ (after $\theta^{\rm RV}_A=\theta_A-h$) and the diagrams $F_k$ are empty for every
*dynamical* diagonal, so the generating word is a **plain gl(4) transfer matrix** at $\theta^{\rm RV}$: no embedding morphism, no inverse transfer matrices. For $\nu_{k+1}=0$ the Bäcklund transforms collapse, $\mathbb T^{(k)}_{(s)}=\mathbb T^{(1)}_{(s)}=Q_1^{[2s]}/Q_1$ (RV2 after 5.24).

**F site (RV $(1,0,0,0)$).** $\mu_{11}\ge\mu_{21}\ge\mu_{31}\in\{0,1\}$, $\bar\mu_k=(\mu_{k1})$; basis = **powers at one point, no shifted products**: $\langle k|_F=\langle0|\,t_1(\theta_F)^k$, $k=0..3$ ($k$ = number of $\mu_{j1}=1$).
Variables: $x_{j1}=\theta_F+h\mu_{j1}$ ($k$ of them at $\theta_F+h$); frozen $x_{22}=x_{32}=\theta_F-h$, $x_{33}=\theta_F-2h$. $B$: $(u-\theta_F)^{3-k}(u-\theta_F-h)^k(u-\theta_F+h)^2(u-\theta_F+2h)$.

**A site (RV $(1,1,1,0)$, $\theta^{\rm RV}=\theta_A-h$).** $\mu_1,\mu_2$ frozen, only $k=3$ dynamical, $\bar\mu_3=(1^a)$, $\mathbb T_{(1^a)}=t_a$: $\langle a|_A=\langle0|\,t_a(\theta_A-h)$, $a=0..3$. Variables: frozen $x_{11}=x_{21}=\theta_A$, $x_{22}=\theta_A-h$;
dynamical $x_{3j}=\theta_A-h+h(m_j+1-j)$, $m=(1^a,0^{3-a})$: $a=0:(\theta_A{-}h,\theta_A{-}2h,\theta_A{-}3h)$, $a=1:(\theta_A,\theta_A{-}2h,\theta_A{-}3h)$, $a=2:(\theta_A,\theta_A{-}h,\theta_A{-}3h)$, $a=3:(\theta_A,\theta_A{-}h,\theta_A{-}2h)$.

**Project words [D].** $\langle a_1,k_2,a_3,k_4|=\langle0|t_{a_1}(\theta_1-h)t_1(\theta_2)^{k_2}t_{a_3}(\theta_3-h)t_1(\theta_4)^{k_4}$ **equals** $\langle\Lambda^{\bf B}|$ exactly; RV's normalised basis differs by a state-independent scalar:
RV2 (5.21) $\langle x|=\langle\Lambda^{\bf B}|\prod_{\alpha,k,j}\Gamma[\nu_1(x^\alpha_{kj})]^{-1}$; RV1 (3.71) $\langle x|=\langle\Lambda|\prod\Phi(X)^{-1}$, $\Phi=\prod_{j\le A}\prod_{r\le m_j}Q_\theta(\theta-\hbar(S+j-r))$ (homogeneous chain).
RV1 (3.49), one rep on all sites: $\langle\Lambda|=\langle0|\prod_{k\le N-A}\prod_\alpha\mathbb T_{\mu^\alpha_k}(\theta_\alpha)$ ($A=1,S=1$: $t_1^n$; $A=3,S=1,N=4$: $t_m$). Shifted products $\prod_jt_1(\theta+jh)$ are the gl(2) Maillet-Niccoli form for $S>1/2$, replaced by $\mathbb T_{1,s}(\theta)$ (RV1 3.12, 3.14); they do not occur here.

**Wave functions [Q].** With $\langle\Omega_\sigma|\Psi\rangle=1$, $\langle\Omega_\sigma|:=\langle0|\prod_{\alpha,k}Q_{\sigma(I_k)}(\theta_\alpha+\hbar\nu^\alpha_{k+1})^{-1}$:
$$\langle x|\Psi\rangle=\prod_{\alpha,k}\det_{i,j\le k}\hat q_{\sigma(i)}(x^\alpha_{kj})\ (5.22,\ 1.1),\qquad\langle\Lambda^{\bf B}|\Psi\rangle=\prod_{\alpha,k}\det_{i,j\le k}Q_{\sigma(i)}(x^\alpha_{kj})\ (5.20),$$
$Q_i=\hat q_i\Gamma[\nu_1]$ (`gaugetr` 5.9), $\hat q_i=z_i^{u/\hbar}q_i$. (RV1 (4.1): $\hat q_I(\theta)^{-L(N-A)}\prod\det_{i\in I,j\le A}\hat q_i(x_{kj})$, any $A$-subset $I$.)
* **F site: single-index Q.** $\langle\Lambda|\Psi\rangle/\langle0|\Psi\rangle=\tau_1(\theta_F)^k$, $\tau_1(\theta_F)=Q_1(\theta_F+h)/Q_1(\theta_F)$, i.e. $\langle x|\Psi\rangle_F\propto\prod_{j\le3}\hat q_1(\theta_F+h\mu_{j1})$ up to $x$-independent vacuum factors; index-independent ($Q_1\to$ any $Q_i$, §3).
* **A site: one $3\times3$ Slater determinant** $\det_{i,j\le3}\hat q_i(x_{3j})$ (the $k=1,2$ factors are $a$-independent); Hodge-dual: $a=0$: $Q_{123}(\theta_A-h)$, $a=3$: $Q_{123}(\theta_A)$ (single dual $Q^4$); $a=1,2$: **gapped Casoratians**, no single-Q form (§3).

## 3. $B$-operator, Q-system, Baxter equation

**B [Q].** $B(u)=\sum_{J_1..J_{n-1}}T[{}^{J_1}_1]\,T^{[2]}[{}^{J_2}_{1\,J_1+1}]\,T^{[4]}[{}^{J_3}_{1\,J_2+1}]\cdots$ (`Bincomp` 2.25, untwisted $T$, $z$-independent) $=GT_1GT_2^{[2]}\cdots+{\rm nilpotent}$ (1.6). $n=4$: $\sum_{k_1\le3}\sum_{|K_2|=2}T_{k_11}(u)\,T[{}^{K_2}_{(1,k_1+1)}](u+h)\,T[{}^{123}_{(1,K_2+1)}](u+2h)$.
Eigenvalues: §2 (4.3), degree 6 per site; at generic $\theta$ all 256 words are $B$-eigen; at the ABJM point $B$ is degenerate.

**Q-system [Q].** $Q_\emptyset=1$; QQ (5.7) $Q_{Iij}Q_I^{[-2]}=Q_{Ii}Q_{Ij}^{[-2]}-Q_{Ij}Q_{Ii}^{[-2]}$, i.e. $Q_{12}(u)=Q_1(u)Q_2(u-h)-Q_2(u)Q_1(u-h)$, $Q_{12..k}(u)=\det_{i,j\le k}Q_i(u-(j-1)h)$;
$Q_I=N_I\hat q_I\prod_{j\le|I|}\Gamma[\hat\nu_j^{[2(1-|I|)]}]$ (5.8), $\hat q_{1234}=1$; Baxter (`qdefn` 5.6): $\det(1-\mathbb Te^{-\hbar\partial_u})Q_i^{[2]}=0$, i.e.
$$\textstyle\sum_{a=0}^4(-1)^at_a(u)Q_i(u+h-ah)=0,\qquad\frac{Q_{1234}(u+h)}{Q_{1234}(u)}=\det G\prod_{j=1}^4\nu_j(u-(4-j)h)=\text{qdet}.$$
$t_\xi=\sum_{\rm SSYT}\prod QEV_{\#(a,s)}(u+\hbar(s-a))$ (5.3), $QEV_k=\frac{Q_{I_{k-1}}^{[-2]}Q_{I_k}^{[2]}}{Q_{I_{k-1}}Q_{I_k}}=z_k\nu_k\frac{\hat q^{[-2]}_{I_{k-1}}\hat q^{[2]}_{I_k}}{\hat q_{I_{k-1}}\hat q_{I_k}}$ (5.10, 5.16), $QEV_r(\theta_\alpha+\hbar\nu^\alpha_r)=0$ (5.15).
Vacuum: $QEV_k=z_k\nu_k$, matching tex $\sum(-1)^a\tau_aD^a=\prod_i(1-V_iD)$ with $V_i=\nu_{5-i}$. Hodge dual $Q^i=\varepsilon^{\bar ii}Q_{\bar i}$ ($Q^4=Q_{123}$); quantisation (5.27)/(5.28).
**Gauge check [D].** The tex's polynomial Baxter equation $\sum_a(-1)^a\tau_aQ(u+h-ah)\prod_{b=a}^3V_4(u-bh)=0$ is exactly RV's gauge $\Phi=\Gamma[\nu_1]$, $\Phi(u+h)=V_4(u)\Phi(u)$; the quoted recursion $\Phi(u)=V_4(u)\Phi(u-h)$ is the same up to $u\to u+h$ (use the displayed equation).

**Baxter at the special points [D].** From (5.15): $t_{a\ge2}(\theta_F)=0$, $t_4(\theta_A-h)=0$ (tex: $g^{(a)}_F(\theta_F)=0$, $V_1^A(\theta_A-h)=0$). Hence for **all four** solutions (i) $Q_j(\theta_F+h)=\tau_1(\theta_F)Q_j(\theta_F)$;
(ii) at $\theta=\theta_A-h$: $Q_j(\theta+h)-\tau_1Q_j(\theta)+\tau_2Q_j(\theta-h)-\tau_3Q_j(\theta-2h)=0$, whose Cramer solution on $j=1,2,3$ is
$$\tau_a(\theta_A-h)=\frac{D_a}{Q_{123}(\theta_A-h)},\quad D_a=\det_{i\le3;\ b\in\{0..3\}\setminus a}Q_i\big(\theta_A-h+h(1-b)\big),\quad D_0=Q_{123}(\theta_A-h),\ D_3=Q_{123}(\theta_A),$$
= $\mathbb T^{(3)}_{(1^a)}$ of (5.12) with $\hat\xi=(0,-1,-2),(1,-1,-2),(1,0,-2),(1,0,-1)$.

## 4. Dual SoV basis, measure, scalar product

**What exists [Q].** GLMR (sl($N$), symmetric reps of spin $s$, compact fundamental $s=-\frac12$; GLMRV is $N=3$): $\mathbb B$ and $\mathbb C(u)=\sum T[{}^{1..N-1}_{j_{N-2}N}]\cdots T[{}^{j_2}_{j_1N}]T[{}^{j_1}_N]$ (same minors, reverse order); the $*$-map ($T_{ij}(u+a)\mapsto T_{ij}(u-a)$, products reversed) fixes all $\mathbb T_{a,1}$ and sends $\mathbb B\to\mathbb C$.
Left basis $\langle x|=\langle\bar0|\prod_{\alpha,j}\mathbb T_{1,n_{\alpha j}}(\theta_\alpha-\frac i2(2s+n_{\alpha j}-1))$ (`lwBeigen`, lowest-weight bra) or $\langle0|\prod_{\alpha,j}\mathbb T_{N-1,n_{\alpha j}}(\theta_\alpha+is-\frac i2(N-n_{\alpha j}-1))$ (`goodBtransfer`, highest-weight bra), $x_{\alpha a}=\theta_\alpha+i(s+n_{\alpha a})$, $n_{\alpha1}\ge..\ge n_{\alpha,N-1}\ge0$, $\langle x|\Psi\rangle=\prod_{\alpha,a}Q_1(x_{\alpha a})$.
**Right (dual) basis:** kets on the highest-weight ket, $|y\rangle=\prod_\alpha\mathbb T^*_{\mu^\alpha}(\theta_\alpha+\frac i2(2s+\mu^\alpha_1-\mu^{\alpha\prime}_1))|0\rangle$, $\mathbb T^*_\mu(u)=\det\mathbb T_{\mu'_j-j+k,1}(u+\frac i2(\mu'_1-\mu_1-\mu'_j+j+k-1))$ (`Ceigenstates`, `Tstardef`: CBR with the shift pattern reversed),
$\mathbb C\to\prod(u-y_{\alpha j})$, $y_{\alpha j}=\theta_\alpha+i(s+\mu^\alpha_j+1-j)$; $\langle\Psi|y\rangle=\prod_\alpha\det_{j,k\le N-1}Q^{j+1}(y_{\alpha k}+\frac i2(N-2))$ (GLMRV `SoVright`, $N=3$: $q^2(Y_1)q^3(Y_2)-q^3(Y_1)q^2(Y_2)$).
$\langle x|y\rangle$ is **not diagonal** but block-triangular in the SoV charge $\mathbf S=\sum n=\sum m$; completeness $1=\sum M_{y,x}|y\rangle\langle x|$, $M=\langle x|y\rangle^{-T}$, $\langle\Phi|\Phi'\rangle=\sum M_{y,x}\langle\Phi|y\rangle\langle x|\Phi'\rangle$; transfer-matrix eigenstates are orthogonal, $\langle\Psi_A|\Psi_B\rangle\propto\delta_{AB}$.
**Measure (`Msl3resNfinal`):** $M_{y,x}=\sum_k{\rm sign}(\sigma)\prod_{a=1}^{N-1}\frac{\Delta(x_{\sigma^{-1}(a)})}{\Delta(\theta)}\prod_{\alpha,a}\frac{r_{\alpha,n_{\alpha a}}}{r_{\alpha,0}}\big|_{\sigma_{\alpha a}=k_{\alpha a}-m_{\alpha a}+a}$, $r_{\alpha,n}=-\frac1{2\pi}\prod_\beta(n+1-i\theta_\alpha+i\theta_\beta)_{2s-1}$ (`resm`; $r_n/r_0=0$ for $n>-2s$, finite sums in the compact case); sl(2): $M_x=\frac{\Delta(x)}{\Delta(\theta)}\prod_\alpha\frac{r_{\alpha n_\alpha}}{r_{\alpha0}}$ (`mucompact`).
**Determinants (`mdef0`, `detproduct`, `detsov`, `PsiPhi`):** brackets $\langle\!\langle f\rangle\!\rangle_\alpha=\sum_nr_{\alpha n}f(x_{\alpha0}+in)=\frac1{2\pi i}\oint\mu_\alpha f$, $\mu_\alpha=\prod_{\beta\ne\alpha}(1-e^{2\pi(u-\theta_\beta-is)})/\prod_kQ_\theta(u+ik)$ (`muc`; fundamental: $\mu_\alpha\sim1/(Q_\theta^+Q_\theta^-)$).
Orthogonality $\det_{(a,\alpha),(b,\beta)}\langle\!\langle Q_1^Au^{\beta-1}\mathcal D^{N-2b}Q^{B,a+1}\rangle\!\rangle_\alpha\propto\delta_{AB}$ ($a,b\le N-1$, $\mathcal Df=f(u+\frac i2)$). If $\langle x|\Theta\rangle=\prod_{\alpha a}F_\alpha(x_{\alpha a})$ and $\langle\Phi|y\rangle=\prod_\alpha\det_{ab}G_{\alpha,a+1}(y_{\alpha b}+\frac i2(N-2))$ ("separable") then
$\langle\Phi|\Theta\rangle=N_0^{-1}\det_{(a,\alpha),(b,\beta)}\langle\!\langle F_\alpha u^{\beta-1}\mathcal D^{N-2b}G_{\alpha,a+1}\rangle\!\rangle_\alpha$ ($N_0=\Delta_\theta^2\prod r_{\alpha0}^2$ for sl(3)); $b_\beta(w),c_\beta(w)$ insertions act as $F_\alpha\to(w-u)^{\delta_{\alpha\beta}}F_\alpha$ (`FBBT`): off-shell Bethe states, different-twist overlaps (`overl`).

**Mixed 4/$\bar4$ chains: not covered [Q: absence].** RV2 gives the *left* basis for arbitrary irreps per site (§2); dual basis, $\mathbb C$, measure and determinants exist only for symmetric/fundamental reps (GLMR, GLMRV) and the $\mathbf 6$ (BBLR, §5).
**[I]** The antipode/$*$-map exchanges $4\leftrightarrow\bar4$ and $Q_i\leftrightarrow Q^i$, so the right basis of A F A F should be RV2's left construction for the mirrored chain F A F A: **A sites carry products of single dual $Q^i$, F sites a $3\times3$ Slater determinant of $Q^2,Q^3,Q^4$** (left-basis roles exchanged).
Counting: 3 dynamical variables per site and $2L_p$ sites match the $6L_p$ charges ($t_1,t_2,t_3$, $2L_p$ each), so the orthogonality determinant should be $6L_p\times6L_p$, rows (site, $Q^{a+1}$), columns ($\beta$, shift $b$). Pole lattices from §2: F $\{\theta_F,\theta_F+h\}$, A $\{\theta_A-3h,..,\theta_A\}$.
Before use: derive $\mu_s(u)$ from adjointness of the alternating-chain Baxter operators and test $\det M_{AB}\propto\delta_{AB}$ on `data/abjm-Qfunctions.wxf`.

## 5. BBLR (2603): final formulas, $\mathbf 6$ of su(4), $\hbar=i$

Basis $\langle x|=\langle0|\prod_\alpha t_{\mu_1^\alpha}(\theta_\alpha)t_{\mu_2^\alpha}(\theta_\alpha)$, $t_{[1,0]}=t_+/Q_\theta^{[-2]}$, $t_{[1,1]}=t_0^-/(Q_\theta^{[-2]}Q_\theta^{[-4]})$ (`eqn:xbasis`, `eqn:TsGT`), $\langle0|=\langle Z^L|e^{V_{13}\mathbb E_{13}}e^{V_{14}\mathbb E_{14}}e^{V_{23}\mathbb E_{23}}e^{V_{24}\mathbb E_{24}}$;
$\langle x|\Psi\rangle=\prod_\alpha\frac{\psi(x_{11},x_{12})\psi(x_{21},x_{22})}{Q_{12}(\theta_\alpha-\frac i2)^2}$, $\psi(x_1,x_2)=Q_1(x_1)Q_2(x_2)-Q_2(x_1)Q_1(x_2)$ (`eqn:xwavefunction`); $\langle\Psi|y\rangle=\prod_\alpha\frac{\psi(y_{11},y_{12})\psi(y_{21},y_{22})}{Q^{34}(\theta_\alpha+\frac i2)^2}$ (`eqn:ywave`), $|\bar y\rangle=\prod\bar t_{\mu_1}\bar t_{\mu_2}|1\rangle$.
Baxter $\mathcal B=Q_\theta^{[2]}Q_\theta D^{-4}-\tau_+Q_\theta^{[2]}D^{-2}+\tau_0^+-\tau_-^{[2]}Q_\theta^{[-2]}D^2+Q_\theta^{[-2]}Q_\theta D^4$ (`eqn:baxter`), $\mathcal B^\dagger$: $\tau_\pm\to\tau_\mp$; adjoint under $\llbracket f\rrbracket_\alpha=\oint\frac{du}{(-2\pi)^\alpha}\frac{e^{2\pi u(\alpha-1)}f}{Q_\theta^{[-2]}Q_\theta Q_\theta^{[2]}}$ (`eqn:Baxadjoint`);
SoV integrals: $\langle\!\langle f\rangle\!\rangle_\alpha=\oint\frac{du}{(-2\pi)^\alpha}\frac{e^{2\pi u(\alpha-1)}f}{Q_\theta^-Q_\theta^+}$ (`eqn:inhombracket`), homogeneous $\langle f\rangle_{\ell,\alpha}=\oint\frac{du}{(-2\pi)^\alpha}\frac{e^{2\pi u(\alpha-1)}f(u)}{(u+\frac i2)^\ell(u-\frac i2)^\ell}$ around $\pm\frac i2$ (`matrixelements` 16).

**Result (`CSoV` 11):** $C^0_\ell=\mathcal N_\ell\,\mathbb W_\omega\mathbb A_{\ell,\kappa}/\sqrt{\mathbb B}$; $\mathcal Q_j=Q_j/\sqrt{Q_{12}(-\frac i2)}$, $\mathcal Q^k=Q^k/\sqrt{Q^{34}(\frac i2)}$, $j\in\{1,2\}$, $k\in\{3,4\}$, $a,b=\pm1$:
$$\mathbb B=\det_{(\alpha,j,k),(\beta,a,b)}\langle\mathcal Q_j^{[a]}u^{\beta-1}\mathcal Q^{k[b]}\rangle_{L,\alpha}\,(12),\ \ \mathbb W_\omega=\det_{\alpha,\beta\le L}\langle\omega^{-iu}u^{\beta-1}\mathcal Q_{12}\rangle_{L,\alpha}\,(13),\ \ \mathbb A_{\ell,\kappa}=\det_{(\gamma,j),(\delta,a)}\langle\kappa^{-iu}u^{\delta-1}\mathcal Q_j^{[a]}\rangle_{\ell,\gamma}\,(14),$$
$$\mathcal N_\ell=\omega^{L/2}\kappa^\ell\,\frac{(1-z_1^2)^{L/2}}{(1-\kappa z_1)^\ell}\,\frac{(1-z_2^2)^{L/2}}{(1-\kappa z_2)^\ell}\,\frac{(1-z_1z_2)^L}{(1-\omega z_1z_2)^L}\ \ (\text{3ptNorm }17;\ z_1z_4=z_2z_3=1;\ \text{up to a phase}).$$
General twists (`app:FullResult`): $C_\ell=\mathcal N^{\rm gen}_\ell\mathbb W_\omega\mathbb A_{\ell,\kappa}/(\mathcal N\mathbb B)^{1/2}$, $\mathcal N=(z_{13}z_{23}z_{14}z_{24})^{-L}(-\frac{\nu_1}{\nu_1^*})^{K+M+\bar K}(\frac{\nu_2}{\nu_2^*})^{K+\bar K}$, $\mathcal N^{\rm gen}_\ell=\frac{\omega^{L/2}\kappa^\ell(1-\omega z_1z_2)^{-L}(1-\kappa z_1)^{-\ell}(1-\kappa z_2)^{-\ell}}{|p_{13}|^{2(K-\bar K)}|p_{14}|^{2\bar K}|p_{23}|^{2(M-K)}}$.
**Overlap (`eqn:finaloverlap`):** $\langle\mathcal W^\ell_{\omega,\kappa}\otimes\mathcal W_\omega^{L-\ell}|\Psi\rangle=\mathbb W_\omega\mathbb A_\ell\,\omega^{L/2}\kappa^\ell(1-\omega z_1z_2)^{-L}(1-\kappa z_1)^{-\ell}(1-\kappa z_2)^{-\ell}$ (unnormalised Q's: $Q_{12}(-\frac i2)^{-(L+\ell)}\det\det$, `eq:overlapHomogeneous`);
$\ell=0$: $\langle\mathcal W^L_\omega|\Psi\rangle=\frac{\omega^{L/2}}{(1-\omega z_1z_2)^L}\det\langle\omega^{-iu}u^{\beta-1}\mathcal Q_{12}\rangle$ (`eq:WPsiOverlap`); inner product $\langle\Psi_A|\Psi_B\rangle=\big(\frac{\mathcal N}{Q^A_{12}(-\frac i2)^2Q_B^{34}(\frac i2)^2}\big)^L\det\mathbb M_{AB}$, $\mathcal N=(z_{13}z_{23}z_{14}z_{24})^{-1}$, $\mathbb M_{AB}=\langle Q_B^{j[s_1]}u^{\beta-1}Q_k^{A[s_2]}\rangle_{L,\alpha}$ (`eqn:sovnorm`); $\|\Psi\|^2=\mathcal N^L\mathbb B\times$phases$\times|p|$-factors (`eqn:finalnorm`).
**Proof structure.** (i) $\langle W_\omega|$ defined by $\langle W_\omega|\Psi\rangle\propto\det\langle\!\langle u^{\beta-1}\omega^{-iu}Q_{12}\rangle\!\rangle$ is a product state $\langle\mathcal W_\omega|=\langle Z|e^{R_{13}\mathbb E_{13}}\cdots$, $R_{13}/V_{13}=\frac{1-\omega z_2z_3}{1-\omega z_1z_2}$ etc. (`eq:WL`; two GT patterns, `eqn:redGTpats`).
(ii) **Localisation** (`eqn:localization`): $\langle Z^L|K\prod_{\alpha\le\ell}t_{[1,0]}(\theta_\alpha)=\langle Z^L|K\prod_{\alpha\le\ell}t^{(\ell)}_{[1,0]}(\theta_\alpha)$ (also $t_{[1,1]}$), so $1-\kappa t_{[1,0]}(\theta_1)+\kappa^2t_{[1,1]}(\theta_1)$ (`eqn:loc1`) rotates site 1 by a group element and gives $K_\kappa(\tau)=\frac{\kappa^\ell}{(1-\kappa z_1)^\ell(1-\kappa z_2)^\ell}\det\langle\kappa^{-iu}u^{\delta-1}\mathcal Q_j^{[a]}\rangle_{\ell,\gamma}$ (`eq:Kkappa`), the su(3) anti-fundamental SoV determinant; $C_{ij}=\frac{1-\kappa z_j}{1-\kappa z_i}R_{ij}$.
Example (`eqn:twistedstructure` 35; $L=2$, $[1,0,1]$): $|C_{\ell=1}|=\frac{(1+\kappa)(z_2-\omega z_1)(1-z_1z_2)}{\sqrt{2z_2}(z_1-z_2)(1+z_2)(1-\omega z_1z_2)}$. $\omega$ enters only via $\omega^{-iu}$ in $\mathbb W$ and prefactors, $\kappa$ only via $\kappa^{-iu}$ in $\mathbb A$; duality invariant under $Q_{1,2}\to Q_{j,k}$.

## 6. How to use: A F A F in the project's conventions

Sites 1..4 = A F A F ($\theta_1,\theta_3$ A; $\theta_2,\theta_4$ F), $G=G_c(z)$, $\langle0|=e_1^T\otimes e_4^T\otimes e_1^T\otimes e_4^T$. $|\Psi\rangle$ has Q's solving $\sum_a(-1)^a\tau_a(u)Q(u+h-ah)=0$ (RV gauge $Q_j=\hat q_j\Gamma[V_4]$, $\hat q_j=z_j^{u/h}q_j$; polynomial gauge of the tex differs by $\Phi$, §3), $Q_{12}(u)=Q_1(u)Q_2(u-h)-Q_2(u)Q_1(u-h)$, $Q_{123}(u)=\det_{i,j\le3}Q_i(u-(j-1)h)$.

**[D] Wave function of the word basis:**
$$\boxed{\langle a_1,k_2,a_3,k_4|\Psi\rangle=\langle0|\Psi\rangle\;\tau_{a_1}(\theta_1-h)\,\tau_1(\theta_2)^{k_2}\,\tau_{a_3}(\theta_3-h)\,\tau_1(\theta_4)^{k_4}}$$
$$\text{F sites }(s=2,4):\ \tau_1(\theta_s)=\frac{Q_j(\theta_s+h)}{Q_j(\theta_s)}\ (\text{any }j=1..4)=V_4(\theta_s)\frac{z_jq_j(\theta_s+h)}{q_j(\theta_s)};\qquad
\text{A sites }(s=1,3):\ \tau_a(\theta_s-h)=\frac{D_a(\theta_s)}{Q_{123}(\theta_s-h)},$$
$D_a(\theta)=\det_{i\le3;\,b\in\{0,1,2,3\}\setminus\{a\}}Q_i(\theta-h+h(1-b))$: $\tau_0=1$; $\tau_3=Q_{123}(\theta_s)/Q_{123}(\theta_s-h)$ (single dual Q, $Q^4$); $\tau_1=\det[Q_i(\theta_s),Q_i(\theta_s{-}2h),Q_i(\theta_s{-}3h)]/Q_{123}(\theta_s-h)$; $\tau_2=\det[Q_i(\theta_s),Q_i(\theta_s{-}h),Q_i(\theta_s{-}3h)]/Q_{123}(\theta_s-h)$ ($i=1,2,3$, columns).
**F site = ratio of one single-index Q at $\theta_F+h$ and $\theta_F$; A site = $3\times3$ Casoratian of $Q_1,Q_2,Q_3$ with columns at $\theta_A-h+h\{1,0,-1,-2\}$ minus the column $1-a$ (dual-Q Wronskian data).**
Separated form (RV normalisation $\langle\Omega_\sigma|\Psi\rangle=1$): $\langle\Lambda|\Psi\rangle=\prod_s\prod_{k\le3}\det_{i,j\le k}Q_i(x^s_{kj})$ with $x^s_{kj}$ of §2, $\langle0|\Psi\rangle=\prod_{s\in F}Q_1Q_{12}Q_{123}(\theta_s)\cdot\prod_{s\in A}Q_1(\theta_s)Q_{12}(\theta_s)Q_{123}(\theta_s-h)$.

**[D] Overlap with a product bra.** $\langle B|=\sum_xc_x\langle x|$ gives $\langle B|\Psi\rangle=\langle0|\Psi\rangle\sum_xc_x\prod_s\tau$; F-only bras (9 words: $a_1=a_3=0$, $k\le2$): $\langle B|\Psi\rangle/\langle0|\Psi\rangle=P(r_2,r_4)$, $r_s=Q_j(\theta_s+h)/Q_j(\theta_s)$, degree $\le2$ in each; A-only: polynomial in the $D_a/D_0$ of sites 1, 3.
**[I] Expected determinant form.** For a bra separable in the dual basis, $\langle B|\Psi\rangle\propto\det$ of the $6L_p\times6L_p$ bracket matrix of §4 (F: single $Q_1$ left, $Q^{2,3,4}$ Slater right; A: mirror); by analogy with BBLR (13), (14) a rotated first $\ell$ pairs would give an $\mathbb A$-type $\kappa$-block and a $\mathbb W$-type $L\times L$ block $\langle\omega^{-iu}u^{\beta-1}Q_{(\cdot)}\rangle$.
Needs (i) the pair $\kappa$-family on F sites (NEXT.md item 1), (ii) the dual basis and measure of §4, (iii) the limit $\theta_A-\theta_F\to2h$ (RV genericity fails, $B$ degenerate).
