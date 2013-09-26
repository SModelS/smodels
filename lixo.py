from Theory.element import Element

el1 = Element('[[[mu+]],[[mu-]]]')
el1.branches[0].masses = [300.,200.]
el1.branches[1].masses = [300.,200.]
el2 = Element('[[[mu+]],[[mu-]]]')
el2.branches[0].masses = [300.,200.]
el2.branches[1].masses = [300.,200.]

print el1 == el2

