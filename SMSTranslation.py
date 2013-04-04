SMSTranslation = {
  "T2": [ [ 2, [1, 0], ['jet'] ], [ 2, [1, 0], ['jet'] ] ],
  "T1": [ [ 2, [2, 0], ['jet','jet'] ], [ 2, [2, 0], ['jet','jet'] ] ],
  "T3w": [ [ 3, [2, 1, 0], ['jet', 'jet', 'W+'] ], [ 2, [2, 0], ['jet','jet'] ] ],
  "T5ww": [ [ 3, [2, 1, 0], ['jet', 'jet', 'W+'] ],[ 3, [2, 1, 0], ['jet', 'jet', 'W+'] ]],
  "T1gg": [ [3, [2,1, 0], ['jet', 'jet', 'photon']], [3, [2,1,0], ['jet', 'jet', 'photon']]], #correct name: T5gg
  "T1lnu": [[3, [2,2,0], ['jet','jet','nu','l+']], [3, [2,2,0], ['jet','jet','nu','l+']]], #T5lnu???
  "T1Lh":  [ [ 2, [2, 0], ['jet','jet'] ], [ 2, [2, 0], ['jet','jet'] ] ], #ist aber nicht einzige moeglichkeit?
  "T1bbbb": [ [ 2, [2,0], ['b','b'] ], [ 2, [2,0], ['b','b'] ] ],
  "T1ctct": [ [ 2, [2,0], ['c','t'] ], [ 2, [2,0], ['c','t'] ] ], #input-lhe falsch?!?!
  "T1tauh": [ [3, [2,2,0],['jet','jet','t+','t-']],[2, [2,0],['jet','jet']]], #lhe falsch + eigentlich T3?!
  "T1taunu":  [ [3, [2,2,0], ['jet', 'jet','nu','ta+']], [3, [2,2,0], ['jet', 'jet','nu','ta+']]], #T5
  "T1tttt": [ [ 2, [2, 0], ['t-','t+'] ], [ 2, [2, 0], ['t-','t+'] ] ],
  "T2bb": [ [ 2, [1,0], ['b'] ], [ 2, [1,0], ['b'] ] ],
  "T2blnu": [ [3, [1,2,0], ['b', 'nu','l+']], [3, [1,2,0], ['b', 'nu','l+']]], #T6blnu???
  "T2blsnu": [ [3, [1,1,0], ['b', 'l-']], [3, [1,1,0], ['b','l+']]], #fehler weil tau im event + T6
  "T2blsnu3":[ [ 2, [1,0], ['l+'] ], [ 2, [1,0], ['l-'] ] ], #lhe falsch
  "T2bw": [ [2, [1,0], ['W+']], [2, [1,0], ['b']]], #input.lhe falsch!!!, name auch?!?!
  "T5wz": [ [ 3, [2, 1,0], ['jet', 'jet', 'Z'] ],[ 3, [2, 1,0], ['jet', 'jet', 'W+'] ]],
  "T2tt": [[2, [1,0], ['t-']],[2, [1,0], ['t+']] ],
  "T2ttww": [[2, [1,0], ['W+']], [2, [1,0], ['W-']]], #stimmt jetzt ueberein mit lhe-file, aber lhe falsch, sollte ja T6ttww sein???
  "T2ttzz": [[2, [1,0], ['Z']], [2, [1,0], ['Z']]], #??? lhe falsch!!! 
  "T3wb": [ [ 3, [2, 1,0], ['jet', 'jet', 'W+'] ], [ 2, [2,0], ['b','b'] ] ],
  "T3z": [ [ 3, [2, 1, 0], ['jet', 'jet', 'Z'] ], [ 2, [2, 0], ['jet','jet'] ] ],
  "T5zz": [ [ 3, [2, 1, 0], ['jet', 'jet', 'Z'] ],[ 3, [2, 1, 0], ['jet', 'jet', 'Z'] ]], #gibt fehler weil bs im lhe
  "T5zzlnu": [[3, [2, 1, 0], ['jet', 'jet', 'Z']], [3, [2, 1, 0], ['b', 'b', 'Z']]], #name oder lhe falsch!!?
  "T5zzoff": [[3, [2, 2, 0], ['jet', 'jet', 'jet', 'jet']], [3, [2, 2, 0], ['b', 'b', 'jet','jet']]], # ??? im lhe nur ein teilchen am 2. vtx
  "T6bbzz": [ [3, [1,1,0], ['b', 'Z']], [3, [1,1,0], ['b', 'Z']] ],
  "TChiChipmSnuSlep": [ [3, [1,1,0], ['l+', 'l-']], [3, [1,1,0], ['l-', 'nu']] ], #lhe falsch, es werden zwei geladene chis erzeugt
  "TChiNuSlep": [ [3, [1,1,0], ['l+', 'l-']], [3, [1,1,0], ['l-', 'nu']] ], #lhe falsch, gibt keinen ll zerfall. welche topo??
  "TChipmSnuSlep": [ [3, [1,1,0], ['l+', 'nu']], [3, [1,1,0], ['nu', 'l-']] ], #gibt fehler weil l=tau im lhe
  "TChiStauStau": [ [3, [1,1,0], ['nu', 'ta+']], [3, [1,1,0], ['nu', 'ta-']] ], #ist das richtig? 
  "TChiww": [[2, [1,0], ['W+']], [2, [1,0], ['W-']]],
  "TChizz": [[2, [1,0], ['Z']], [2, [1,0], ['Z']]], #lhe falsch
  "TChiwz": [[2, [1,0], ['Z']], [2, [1,0], ['W-']]],
  "TGQ": [], #wie soll das ausschauen?
  "TSlepSlep": [ [2, [1,0], ['l-']], [2, [1,0], ['l+']]],
  "TSlepSnu": [ [2, [1,0], ['l+']], [2, [1,0], ['l-']]], #was ist der unterschied zu tslepslep?
}
