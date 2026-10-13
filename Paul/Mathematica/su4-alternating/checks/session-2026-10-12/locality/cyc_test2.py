from cyc import *
cy = Cyclic(4, 1, 0)
print('cyclic r', cy.r)
Lb = cy.loc_basis_coords()
print('loc dim', Lb.shape[0])
ops = [cy.op(f) for f in cy.F]
print('built ops', len(ops))
# check T membership in loc space: T coords in span of Lb
for k, t in cy.Tc.items():
    print(k, in_span(Lb, t, cy.p))
