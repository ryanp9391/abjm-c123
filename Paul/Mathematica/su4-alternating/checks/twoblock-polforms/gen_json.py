import json, io, contextlib, collections, numpy as np, datetime, os
from fractions import Fraction

with contextlib.redirect_stdout(io.StringIO()):
    from symm import coeffs, res
    from polyverify import Fpoly, reduce_poly, NAMES
from core import S
TH = 1e-9
VN = ['v2_1', 'v2_2', 'v2_3', 'v2_4', 'u2_1', 'u2_2', 'u2_3', 'u2_4',
      'v3_1', 'v3_2', 'v3_3', 'v3_4', 'u3_1', 'u3_2', 'u3_3', 'u3_4']
OUT = r'C:\Users\pryan\git-projects\abjm-c123\Paul\Mathematica\su4-alternating\data\abjm-twoblock-forms.json'


def cx(z, nd=15):
    return [float('%.*g' % (nd, z.real)), float('%.*g' % (nd, z.imag))]


def cstr(z, nd=12):
    re_, im_ = float('%.*g' % (nd, z.real)), float('%.*g' % (nd, z.imag))
    if abs(im_) < 1e-13 * max(1, abs(re_)):
        return repr(re_)
    if abs(re_) < 1e-13 * max(1, abs(im_)):
        return repr(im_) + '*I'
    return '(%r%+r*I)' % (re_, im_)


def rat(x, maxden=5000, tol=1e-11):
    f = Fraction(x).limit_denominator(maxden)
    return f if abs(float(f) - x) < tol * max(1, abs(x)) else None


def invexpr(terms):
    parts = []
    for t in terms:
        mon = '*'.join(n if e == 1 else '%s**%d' % (n, e) for n, e in zip(NAMES, t['exp']) if e)
        parts.append(cstr(t['c']) + ('*' + mon if mon else ''))
    return ' + '.join(parts) if parts else '0'


def polystr(p, tol=1e-12, nd=12):
    terms = []
    for e, c in sorted(p.items(), reverse=True):
        if abs(c) < tol:
            continue
        mon = '*'.join(VN[i] if ee == 1 else '%s**%d' % (VN[i], ee) for i, ee in enumerate(e) if ee)
        terms.append(cstr(c, nd) + ('*' + mon if mon else ''))
    return ' + '.join(terms) if terms else '0'


def classify(supp, nm):
    if nm == 0:
        return 'kinematic_zero'
    if not supp:
        return 'dynamical_zero'
    s = supp[0][0] + supp[0][1]
    if s == 0:
        return 'no_pairing_dependence'
    e8min = min(a for a, b in supp)
    e9min = min(b for a, b in supp)
    if e8min >= 1 and e9min >= 1:
        return 'needs_both_pairings'
    if e8min == 0 and e9min >= 1:
        return 'survives_p23=0_only'
    if e9min == 0 and e8min >= 1:
        return 'survives_p32=0_only'
    return 'survives_either_pairing_zero'


DESC = (
    "Tree-level two-block overlaps F_Psi(v2,u2;v3,u3;M) = <(v2 x u2)^(L-M) (x) (v3 x u3)^M | Psi0> of the saved ABJM su(4) highest-weight cyclic states "
    "(data/abjm-hw-states.wxf), L = 2,3,4 pairs, M = 0..L. Bras are row vectors, bilinear (no conjugation), site 1 = A, A site carries v, F site carries u, pair k = sites (2k-1,2k); "
    "the first L-M pairs carry (v2,u2), the last M pairs (v3,u3). F is a polynomial of multidegree (L-M,L-M,M,M) in (v2,u2,v3,u3). "
    "HW trick: F is invariant under (v,u) -> (v n, u n^-T) for upper-unitriangular n (all four vectors) and torus covariant with weight w (A e_i: +eps_i, F e_j: -eps_j). "
    "N+ invariants: a2=v2_1, a3=v3_1, d2=u2_4, d3=u3_4, mv=v2_1*v3_2-v2_2*v3_1, mu=u2_3*u3_4-u2_4*u3_3, p22=v2.u2, p23=v2.u3, p32=v3.u2, p33=v3.u3. "
    "Key 'bps' = polarisations with p22=p33=0 imposed (monomials without p22,p33; unique and complete: the dimension of the N+-invariant space equals the number of monomials, checked for all weights, M, L<=4); "
    "key 'full' = no BPS condition (all ten invariants). Each entry has 'terms' = list of {exp: exponents of (a2,a3,d2,d3,mv,mu,p22,p23,p32,p33), c: [re,im]} and 'expr' = sympy-parsable string in the invariants. "
    "Reduced variables x=p23/(a2*d3), y=p32/(a3*d2), muv=mv/(a2*a3), muu=mu/(d2*d3): F/F_vac = sum_k c_k x^e8 y^e9 muv^e5 muu^e6 (bps) with F_vac = (a2*d2)^(L-M)*(a3*d3)^M; see 'reduced_bps'. "
    "Coefficients are numerical (states are accurate to about 1e-13; invariant-fit residual <= 3.4e-15 over all 184 states and all M); the phase of each state is the saved one (largest component real positive). "
    "'F16_full' and 'F16_bps_nf' (L=2,3 only) are sympy strings of the full 16-variable polynomial and of its normal form modulo v2.u2, v3.u3 (leading terms v2_1*u2_1, v3_1*u3_1 eliminated). "
    "'class' classifies the BPS form: kinematic_zero (no invariant monomial of that weight/degree), dynamical_zero (monomials exist, all coefficients vanish), "
    "no_pairing_dependence (p23,p32 absent, s=0), needs_both_pairings (every term has e8>=1 and e9>=1), survives_p23=0_only (some term has e8=0, none has e9=0), "
    "survives_p32=0_only, survives_either_pairing_zero. Generated 2026-10-12 by the python exploration behind twoblock-polarisation-forms.wb "
    "(the notebook recomputes L=2,3 symbolically and L=4 pointwise and checks against this file)."
)


def main():
    out = {'description': DESC, 'generated': str(datetime.date.today())}
    out['conventions'] = {
        'bra_action': 'A site v -> v g, F site u -> u g^{-T}; pairing v.u invariant',
        'blocks': 'pairs 1..L-M carry (v2,u2), pairs L-M+1..L carry (v3,u3)',
        'note_polarisation': {'v2': [1, 0, 0, 0], 'u2': [0, 0, 0, 1], 'v3': [1, 0, 0, -1], 'u3': [1, 0, 0, 1],
                              'x': 1, 'y': -1, 'muv': 0, 'muu': 0,
                              'C_M_equals': 'F/F_vac at this point = saved C[2(L-M)]'},
        'structure_count': 'N(L,M,[p,q,r]) = max(0, min(s,n2-e5,n3-e6) - max(0,s-n3+e5,s-n2+e6) + 1); e5=(2q+r-p)/4, e6=(p+2q-r)/4 (muv, muu powers), s=L-(p+2q+r)/2=e8+e9 (total pairing power), n2=L-M, n3=M',
        'exponents_bps': 'a2^e1 a3^e2 d2^e3 d3^e4 mv^e5 mu^e6 p23^e8 p32^e9 with e1=n2-e5-e8, e2=n3-e5-e9, e3=n2-e6-e9, e4=n3-e6-e8',
        'exponent_order': list(NAMES)}
    for L in (2, 3, 4):
        lst = []
        for k, s in enumerate(S[L]):
            ent = {'index': k, 'dyn': [int(x) for x in s['dyn']], 'w': [int(x) for x in s['w']],
                   'E': float(s['E']), 'C': [cx(complex(c)) for c in s['C']], 'M': {}}
            for M in range(L + 1):
                rec = {}
                for tag, bps in (('bps', True), ('full', False)):
                    r = res[(L, k, M, bps)]
                    terms = []
                    if r['nm']:
                        for m, c in zip(r['mons'], r['c']):
                            if abs(c) > TH:
                                terms.append({'exp': [int(m.get(n, 0)) for n in NAMES], 'c': c})
                    rec[tag] = terms
                supp = sorted((t['exp'][7], t['exp'][8]) for t in rec['bps'])
                rb = res[(L, k, M, True)]
                rf = res[(L, k, M, False)]
                e = {}
                e['class'] = classify(supp, rb['nm'])
                e['n_monomials_allowed_bps'] = int(rb['nm'])
                e['n_monomials_allowed_full'] = int(rf['nm'])
                e['fit_residual_bps'] = float(rb['res']) if rb['nm'] else None
                e['max_abs_F_if_no_monomial'] = None if rb['nm'] else float(rb['Fmax'])
                for tag in ('bps', 'full'):
                    e[tag] = {'terms': [{'exp': t['exp'], 'c': cx(t['c'])} for t in rec[tag]],
                              'expr': invexpr(rec[tag])}
                red = []
                for t in rec['bps']:
                    ex = t['exp']
                    rr = rat(abs(t['c']) ** 2)
                    red.append({'e5': ex[4], 'e6': ex[5], 'e8': ex[7], 'e9': ex[8], 'c': cx(t['c']),
                                'abs2_rational': (str(rr) if rr is not None else None)})
                e['reduced_bps'] = red
                if L <= 3:
                    F = Fpoly(s['psi'], L, M)
                    e['F16_full'] = polystr(F)
                    e['F16_bps_nf'] = polystr(reduce_poly(F))
                    e['n_terms_F16_full'] = len([1 for c in F.values() if abs(c) > 1e-12])
                ent['M'][str(M)] = e
            lst.append(ent)
        out['L%d' % L] = lst
    json.dump(out, open(OUT, 'w'), separators=(',', ':'))
    print(os.path.getsize(OUT) / 1e6, 'MB')


main()
