#include "pxs_gausq_3.h"

#define C0(a, b, c, d, e, f) C0i(cc0, a, b, c, d, e, f)
#define C1(a, b, c, d, e, f) C0i(cc1, a, b, c, d, e, f)
#define C2(a, b, c, d, e, f) C0i(cc2, a, b, c, d, e, f)
#define C00(a, b, c, d, e, f) C0i(cc00, a, b, c, d, e, f)
#define C11(a, b, c, d, e, f) C0i(cc11, a, b, c, d, e, f)
#define C12(a, b, c, d, e, f) C0i(cc12, a, b, c, d, e, f)
#define C22(a, b, c, d, e, f) C0i(cc22, a, b, c, d, e, f)
#define Power std::pow
//#define Conjugate conj

ComplexType ME_us_gGX_qQQ_single(POLE pIEPS, bool sc, bool uc, bool axial, int itsq, int itq,
                                 int color_flow, double Q2, double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // MG = -MG;

  double SS = sc, UU = uc, AXG = axial;
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){
  int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
  int iq = itq;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType iLp = conj(params->CHSQq[ch][isq][iq].R);
  ComplexType iRp = conj(params->CHSQq[ch][isq][iq].L);

  auto MQi = params->mSQ[isq];
  auto MQis = pow2(MQi);

  auto Mqi = params->mq[itq];
  auto Mqis = pow2(Mqi);

  ComplexType iLG = (params->GLSQq[isq][iq].L);
  ComplexType iRG = (params->GLSQq[isq][iq].R);
  ComplexType LGp = conj(params->GLSQq[sq][q].R);
  ComplexType RGp = conj(params->GLSQq[sq][q].L);
  if (color_flow) {
    auto tiLp = conj(iRp);
    auto tiRp = conj(iLp);
    auto tiLG = conj(iRG);
    auto tiRG = conj(iLG);
    iLp = tiLp;
    iRp = tiRp;
    iLG = tiLG;
    iRG = tiRG;
  }
  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C0(0, t, MXs, MQis, MQis, Mqis)
#define syFC2 C00(0, MXs, t, MQis, MQis, Mqis)
#define syFC3 C1(0, MXs, t, MQis, MQis, Mqis)
#define syFC4 C11(0, MXs, t, MQis, MQis, Mqis)
#define syFC5 C12(0, MXs, t, MQis, MQis, Mqis)
#define syFC6 C2(0, MXs, t, MQis, MQis, Mqis)
#define syFC7 C22(0, MXs, t, MQis, MQis, Mqis)

  _EPS0_(
      ret, auto x0 = L * RGp; auto x1 = iLG * iRp; auto x2 = x0 * x1; auto x3 = pow(MXs, 3);
      auto x4 = UU * x3; auto x5 = x2 * x4; auto x6 = 5. * s; auto x7 = iLp * iRG;
      auto x8 = LGp * R; auto x9 = x7 * x8; auto x10 = x4 * x9; auto x11 = AXG * UU;
      auto x12 = x11 * x2; auto x13 = pow(u, 4); auto x14 = 4. * x13; auto x15 = x12 * x14;
      auto x16 = AXG * SS; auto x17 = x16 * x2; auto x18 = 2. * x13; auto x19 = x17 * x18;
      auto x20 = pow(MXs, 2); auto x21 = UU * x20; auto x22 = x2 * x21; auto x23 = MUs * x22;
      auto x24 = s * x23; auto x25 = x21 * x9; auto x26 = s * x25; auto x27 = 10. * MUs;
      auto x28 = MXs * s; auto x29 = UU * x2; auto x30 = x28 * x29; auto x31 = pow(MUs, 2);
      auto x32 = 5. * x31; auto x33 = UU * x31; auto x34 = x28 * x9; auto x35 = x33 * x34;
      auto x36 = SS * x20; auto x37 = x2 * x36; auto x38 = MUs * x37; auto x39 = s * x38;
      auto x40 = SS * x2; auto x41 = x31 * x40; auto x42 = x28 * x41; auto x43 = x36 * x9;
      auto x44 = MUs * x43; auto x45 = s * x44; auto x46 = SS * x31; auto x47 = x34 * x46;
      auto x48 = x11 * x3; auto x49 = x2 * x48; auto x50 = x48 * x9; auto x51 = pow(MX, 4);
      auto x52 = UU * x51; auto x53 = x52 * x9; auto x54 = MUs * x53; auto x55 = 6. * s;
      auto x56 = x11 * x51; auto x57 = x56 * x9; auto x58 = MUs * x57; auto x59 = x17 * x20;
      auto x60 = MUs * x59; auto x61 = 4. * s; auto x62 = 4. * x31; auto x63 = x28 * x62;
      auto x64 = x16 * x9; auto x65 = x20 * x64; auto x66 = MUs * x65; auto x67 = x51 * x64;
      auto x68 = MUs * x67; auto x69 = 2. * s; auto x70 = x28 * x32; auto x71 = x11 * x9;
      auto x72 = x11 * x20; auto x73 = x2 * x72; auto x74 = MUs * x73; auto x75 = 10. * s;
      auto x76 = x72 * x9; auto x77 = MUs * x76; auto x78 = pow(s, 3); auto x79 = UU * x78;
      auto x80 = u * x0; auto x81 = x1 * x80; auto x82 = x79 * x81; auto x83 = u * x8;
      auto x84 = x7 * x83; auto x85 = x79 * x84; auto x86 = MXs * x2; auto x87 = pow(s, 2);
      auto x88 = pow(u, 2); auto x89 = SS * x88; auto x90 = x2 * x89; auto x91 = x87 * x90;
      auto x92 = MXs * x9; auto x93 = -x79 * x92; auto x94 = x89 * x9; auto x95 = x87 * x94;
      auto x96 = UU * x88; auto x97 = x2 * x96; auto x98 = 3. * x87; auto x99 = x9 * x96;
      auto x100 = x11 * x86; auto x101 = x11 * x92; auto x102 = Mqi * x0; auto x103 = iLG * iLp;
      auto x104 = x102 * x103; auto x105 = MX * x79; auto x106 = iRG * iRp; auto x107 = Mqi * x8;
      auto x108 = x106 * x107; auto x109 = x11 * x78; auto x110 = x109 * x81;
      auto x111 = x109 * x84; auto x112 = MG * MX; auto x113 = x0 * x7; auto x114 = x112 * x113;
      auto x115 = x114 * x79; auto x116 = x1 * x8; auto x117 = x112 * x116; auto x118 = x117 * x79;
      auto x119 = x11 * x88; auto x120 = x119 * x2; auto x121 = 7. * x87; auto x122 = x119 * x9;
      auto x123 = x17 * x88; auto x124 = x64 * x88; auto x125 = x109 * x114;
      auto x126 = x109 * x117; auto x127 = MX * x104; auto x128 = MX * x108; auto x129 = UU * x87;
      auto x130 = x7 * x80; auto x131 = x112 * x130; auto x132 = x129 * x131; auto x133 = 3. * x132;
      auto x134 = x11 * x87; auto x135 = x131 * x134; auto x136 = 7. * x135; auto x137 = 2. * x87;
      auto x138 = 6. * x87; auto x139 = x129 * x81; auto x140 = 5. * MXs; auto x141 = x129 * x84;
      auto x142 = MUs * UU; auto x143 = x142 * x81; auto x144 = 4. * x87; auto x145 = SS * x81;
      auto x146 = x137 * x145; auto x147 = SS * x84; auto x148 = x137 * x147; auto x149 = x16 * x81;
      auto x150 = MXs * x137; auto x151 = x16 * x84; auto x152 = x1 * x83; auto x153 = x112 * x152;
      auto x154 = x129 * x153; auto x155 = SS * x112; auto x156 = x130 * x155;
      auto x157 = x152 * x155; auto x158 = MUs * x11; auto x159 = x158 * x81;
      auto x160 = x134 * x81; auto x161 = 9. * MXs; auto x162 = x134 * x84; auto x163 = x131 * x16;
      auto x164 = x153 * x16; auto x165 = x134 * x153; auto x166 = pow(u, 3); auto x167 = s * x166;
      auto x168 = x12 * x167; auto x169 = x167 * x9; auto x170 = x11 * x169; auto x171 = x28 * x97;
      auto x172 = SS * x86; auto x173 = 2. * x172; auto x174 = MUs * x173; auto x175 = x174 * x87;
      auto x176 = SS * x92; auto x177 = 2. * x176; auto x178 = MUs * x177; auto x179 = x178 * x87;
      auto x180 = x167 * x17; auto x181 = x16 * x169; auto x182 = UU * x86; auto x183 = MUs * x144;
      auto x184 = x28 * x90; auto x185 = UU * x92; auto x186 = s * x96; auto x187 = x114 * x186;
      auto x188 = x117 * x186; auto x189 = x123 * x28; auto x190 = x16 * x86;
      auto x191 = MUs * x137; auto x192 = x16 * x92; auto x193 = x191 * x192; auto x194 = x69 * x89;
      auto x195 = x114 * x194; auto x196 = x117 * x194; auto x197 = x120 * x28;
      auto x198 = x114 * x16; auto x199 = x61 * x88; auto x200 = x117 * x16; auto x201 = s * x119;
      auto x202 = x114 * x201; auto x203 = x117 * x201; auto x204 = x167 * x29;
      auto x205 = UU * x169; auto x206 = x204 + x205; auto x207 = MUs * x176; auto x208 = s * x89;
      auto x209 = x102 * x106; auto x210 = MG * x209; auto x211 = x103 * x107;
      auto x212 = MG * x211; auto x213 = x114 * x208; auto x214 = s * x88; auto x215 = x16 * x210;
      auto x216 = x16 * x211; auto x217 = MG * x216; auto x218 = x167 * x40; auto x219 = SS * x169;
      auto x220 = x218 + x219; auto x221 = 2. * x185; auto x222 = MUs * x87; auto x223 = x69 * x96;
      auto x224 = x119 * x69; auto x225 = x69 * x88; auto x226 = x114 * x223 + x117 * x223;
      auto x227 = 8. * x87; auto x228 = 5. * x87; auto x229 = 5. * x182; auto x230 = 5. * x185;
      auto x231 = 3. * x172; auto x232 = MUs * x98; auto x233 = MUs * x228; auto x234 = x119 * x75;
      auto x235 = x6 * x88; auto x236 = 8. * MUs; auto x237 = 8. * x77; auto x238 = x12 * x62;
      auto x239 = 2. * x166; auto x240 = 3. * x60; auto x241 = 3. * x31; auto x242 = x17 * x241;
      auto x243 = 3. * x66; auto x244 = x241 * x64; auto x245 = MG * x11; auto x246 = x209 * x245;
      auto x247 = 4. * x166; auto x248 = x13 * x64; auto x249 = x11 * x114; auto x250 = x11 * x117;
      auto x251 =
          x14 * x71 + x15 - x19 + x198 * x239 + x200 * x239 - x247 * x249 - x247 * x250 - 2. * x248;
      auto x252 = x166 * x172; auto x253 = x166 * x176; auto x254 = MUs * x40;
      auto x255 = x166 * x254; auto x256 = x166 * x9; auto x257 = MUs * SS; auto x258 = x256 * x257;
      auto x259 = 12. * x166; auto x260 = x158 * x2; auto x261 = x166 * x260;
      auto x262 = x158 * x256; auto x263 = 3. * x9; auto x264 = x263 * x33; auto x265 = MUs * x17;
      auto x266 = x16 * x256; auto x267 = 4. * MUs; auto x268 = MG * x89; auto x269 = pow(MX, 3);
      auto x270 = x113 * x269; auto x271 = 4. * x270; auto x272 = 6. * x166; auto x273 = MG * x271;
      auto x274 = x11 * x31; auto x275 = x263 * x274; auto x276 = 2. * x31; auto x277 = x28 * x64;
      auto x278 = x147 * x51; auto x279 = SS * x51; auto x280 = x2 * x279; auto x281 = MUs * x280;
      auto x282 = x279 * x9; auto x283 = MUs * x282; auto x284 = x72 * x81; auto x285 = x72 * x84;
      auto x286 = x56 * x84; auto x287 = 4. * MXs; auto x288 = x158 * x84; auto x289 = pow(MX, 6);
      auto x290 = x11 * x289; auto x291 = x2 * x290; auto x292 = 2. * x291; auto x293 = x290 * x9;
      auto x294 = 2. * x293; auto x295 = -x292 - x294; auto x296 = x145 * x31;
      auto x297 = x147 * x31; auto x298 = x296 + x297; auto x299 = x2 * x56; auto x300 = 2. * MUs;
      auto x301 = x299 * x300; auto x302 = 2. * x56; auto x303 = x302 * x86; auto x304 = x302 * x92;
      auto x305 = x301 + x303 + x304; auto x306 = x17 * x51; auto x307 = MUs * x306;
      auto x308 = -x68; auto x309 = -x307 + x308; auto x310 = x300 * x57; auto x311 = x172 * x31;
      auto x312 = x176 * x31; auto x313 = x310 - x311 - x312 - x38 - x44; auto x314 = x36 * x81;
      auto x315 = x147 * x20; auto x316 = x151 * x51; auto x317 = x314 + x315 + x316;
      auto x318 = x149 * x300; auto x319 = x151 * x300; auto x320 = MXs * x318 + MXs * x319;
      auto x321 = x12 * x166; auto x322 = x11 * x256; auto x323 = x279 * x81; auto x324 = 2. * x256;
      auto x325 = x16 * x324; auto x326 = x142 * x2; auto x327 = x142 * x9; auto x328 = x190 * x88;
      auto x329 = x192 * x88; auto x330 = x149 * x51; auto x331 = 2. * x330; auto x332 = x56 * x81;
      auto x333 = 4. * x332; auto x334 = 8. * x119; auto x335 = 3. * x260; auto x336 = x86 * x89;
      auto x337 = MXs * x94; auto x338 = x123 * x300; auto x339 = x124 * x300;
      auto x340 = 2. * x336 + 2. * x337 - x338 - x339; auto x341 = x300 * x90;
      auto x342 = x300 * x94; auto x343 = -x341 - x342; auto x344 = 2. * syFC2;
      auto x345 = 4. * x51; auto x346 = 2. * x20; auto x347 = x123 * x346; auto x348 = x124 * x346;
      auto x349 = 6. * x20; auto x350 = x119 * x86;
      auto x351 = MUs * x325 + x239 * x265 - x252 - x253 + x255 + x258 - 4. * x261 - 4. * x262;
      auto x352 = 2. * syFC5; auto x353 = SS * x3; auto x354 = 2. * x289;
      auto x355 = 2. * pow(MX, 8); auto x356 = MXs * x296; auto x357 = MXs * x297;
      auto x358 = 4. * x291; auto x359 = 4. * x290; auto x360 = 4. * x48; auto x361 = x149 * x20;
      auto x362 = x151 * x20; auto x363 = pow(MUs, 3); auto x364 = x145 * x363 + x147 * x363;
      auto x365 = x274 * x81; auto x366 = x274 * x84; auto x367 = 4. * x366;
      auto x368 = x149 * x276; auto x369 = x151 * x276;
      auto x370 = -MXs * x367 + MXs * x368 + MXs * x369 - x287 * x365; auto x371 = x276 * x90;
      auto x372 = x276 * x94; auto x373 = x123 * x31; auto x374 = x124 * x31; auto x375 = 3. * MUs;
      auto x376 = 3. * MXs; auto x377 = MUs * x147; auto x378 = x120 * x276;
      auto x379 = x122 * x276; auto x380 = x119 * x92; auto x381 = x279 * x300;
      auto x382 = 2. * x290; auto x383 = MUs * x16; auto x384 = x209 * x383; auto x385 = MG * x384;
      auto x386 = x216 * x51; auto x387 = MG * MUs; auto x388 = x210 * x56; auto x389 = 2. * MXs;
      auto x390 = x212 * x56; auto x391 = 2. * x113; auto x392 = pow(MX, 7);
      auto x393 = x245 * x392; auto x394 = 2. * x116; auto x395 = pow(MX, 5);
      auto x396 = x257 * x395; auto x397 = MG * x396; auto x398 = x16 * x395; auto x399 = MG * x398;
      auto x400 = MUs * x399; auto x401 = x158 * x395; auto x402 = MG * x401;
      auto x403 = x245 * x395; auto x404 = MXs * x116; auto x405 = 2. * x403;
      auto x406 = x113 * x389 * x403 - x113 * x400 - x116 * x400 + 2. * x357 + x364 - x391 * x393 +
                  x391 * x397 + x391 * x402 - x393 * x394 + x394 * x397 + x394 * x402 + x404 * x405;
      auto x407 = 2. * syFC6; auto x408 = x2 * x353; auto x409 = MUs * x408;
      auto x410 = x172 * x363; auto x411 = x353 * x9; auto x412 = MUs * x411;
      auto x413 = x176 * x363; auto x414 = x276 * x37; auto x415 = x276 * x43;
      auto x416 = 2. * x280; auto x417 = x31 * x416; auto x418 = 2. * x282; auto x419 = x31 * x418;
      auto x420 = x269 * x46; auto x421 = x306 * x31; auto x422 = x31 * x67;
      auto x423 = x276 * x299; auto x424 = x276 * x57; auto x425 = 2. * x57; auto x426 = x20 * x425;
      auto x427 = x174 * x51; auto x428 = x279 * x92; auto x429 = x300 * x428;
      auto x430 = x104 * x269; auto x431 = x257 * x430; auto x432 = x108 * x269;
      auto x433 = x257 * x432; auto x434 = x190 * x51; auto x435 = MUs * x434;
      auto x436 = x192 * x51; auto x437 = MUs * x436; auto x438 = x267 * x56;
      auto x439 = x438 * x86; auto x440 = x438 * x92; auto x441 = 2. * x268;
      auto x442 = x116 * x269; auto x443 = 3. * x190; auto x444 = 3. * x192; auto x445 = x127 * x16;
      auto x446 = x128 * x16; auto x447 = 2. * x16; auto x448 = MG * x447; auto x449 = x442 * x448;
      auto x450 = 2. * x119; auto x451 = MG * x450; auto x452 = x11 * x127; auto x453 = x211 * x245;
      auto x454 = x11 * x128; auto x455 = MG * x46; auto x456 = x270 * x455;
      auto x457 = x442 * x455; auto x458 = x382 * x86; auto x459 = x257 * x270;
      auto x460 = MG * x459; auto x461 = MXs * x460; auto x462 = x257 * x442; auto x463 = MG * x462;
      auto x464 = MXs * x463; auto x465 = 2. * syFC7; auto x466 = x257 * x9;
      auto x467 = x114 * x129; auto x468 = x117 * x129; auto x469 = x114 * x134;
      auto x470 = x117 * x134; auto x471 = UU * x112; auto x472 = x130 * x471;
      auto x473 = x11 * x131; auto x474 = x142 * x87; auto x475 = x114 * x257;
      auto x476 = x117 * x257; auto x477 = x158 * x87; auto x478 = 2. * syFC4;
      auto x479 = 4. * x293; auto x480 = 2. * x299; auto x481 = -x20 * x480;
      auto x482 = x409 + x410 + x412 + x413 + x414 + x415 + x421 + x422 - x423 - x424 - x426 +
                  x435 + x437 - x439 - x440 + x481;
      auto x483 = 2. * x279; auto x484 = MG * Mqi; auto x485 = x106 * x80; auto x486 = x484 * x485;
      auto x487 = x16 * x486; auto x488 = x103 * x83; auto x489 = x484 * x488; auto x490 = 2. * SS;
      auto x491 = x395 * x490; auto x492 = MG * x491; auto x493 = 5. * x190;
      auto x494 = 5. * MUs * x329 + MUs * x493 * x88 - x130 * x399 - x130 * x405 + x130 * x492 -
                  x152 * x399 - x152 * x405 + x152 * x492 - x27 * x350 - x27 * x380 + x371 + x372 +
                  x373 + x374 - x378 - x379;
      auto x495 = 2. * x104; auto x496 = x11 * x392; auto x497 = 2. * x108; auto x498 = x104 * x398;
      auto x499 = x108 * x398; auto x500 = x209 * x36; auto x501 = MXs * SS;
      auto x502 = x209 * x501; auto x503 = MG * x502; auto x504 = x211 * x36;
      auto x505 = x211 * x501; auto x506 = MG * x505; auto x507 = x11 * x395;
      auto x508 = x389 * x507; auto x509 = MUs * x292 + MUs * x294 + x382 * x92;
      auto x510 = MG * x269; auto x511 = x130 * x257; auto x512 = x501 * x510;
      auto x513 = x130 * x510; auto x514 = x16 * x300; auto x515 = 4. * x158;
      auto x516 = x245 * x269; auto x517 = x287 * x516; auto x518 = x152 * x510;
      auto x519 = x182 * x31; auto x520 = MUs * x25; auto x521 = x185 * x31; auto x522 = x100 * x31;
      auto x523 = x101 * x31; auto x524 = x190 * x31; auto x525 = x192 * x31; auto x526 = MX * x142;
      auto x527 = x104 * x526; auto x528 = x108 * x526; auto x529 = MX * x257;
      auto x530 = x104 * x529; auto x531 = MXs * x530; auto x532 = x211 * x257;
      auto x533 = MG * x532; auto x534 = x108 * x529; auto x535 = MXs * x534;
      auto x536 = MUs * x445; auto x537 = MUs * x216; auto x538 = MG * MXs; auto x539 = MUs * x446;
      auto x540 = x127 * x158; auto x541 = x128 * x158; auto x542 = s * syFC1;
      auto x543 = x300 * x97; auto x544 = x86 * x96; auto x545 = x300 * x99; auto x546 = 2. * x190;
      auto x547 = x120 * x300; auto x548 = x122 * x300; auto x549 = x341 + x342;
      auto x550 = MXs * UU; auto x551 = x131 * x550; auto x552 = MXs * x16; auto x553 = MXs * x11;
      auto x554 = x131 * x553; auto x555 = -x131 * x501 + x131 * x552 + x551 - x554;
      auto x556 = MXs * x445; auto x557 = MXs * x446; auto x558 = x100 * x241 - 3. * x519;
      auto x559 = s * syFC3; auto x560 = x2 * x52; auto x561 = MUs * x560; auto x562 = x52 * x86;
      auto x563 = MXs * x53; auto x564 = MUs * x299; auto x565 = 6. * x56; auto x566 = MG * x142;
      auto x567 = 2. * x566; auto x568 = MG * UU; auto x569 = x270 * x568; auto x570 = x442 * x568;
      auto x571 = 2. * x307; auto x572 = MG * x270; auto x573 = 6. * x158; auto x574 = x245 * x270;
      auto x575 = 6. * MXs; auto x576 = MG * x442; auto x577 = x245 * x442; auto x578 = UU * x395;
      auto x579 = MG * x578; auto x580 = x116 * x579; auto x581 = 6. * x403;
      auto x582 = -x116 * x581 + 2. * x580; auto x583 = x16 * x267;
      auto x584 = 4. * x460 + 4. * x463 - x572 * x583 - x576 * x583; auto x585 = UU * x289;
      auto x586 = -x2 * x585 + x291 + x293 - x585 * x9; auto x587 = s * syFC5;
      auto x588 = x52 * x81; auto x589 = 4. * SS; auto x590 = x510 * x589; auto x591 = x130 * x590;
      auto x592 = x152 * x590; auto x593 = 16. * MUs; auto x594 = x112 * x511;
      auto x595 = x153 * x257; auto x596 = x131 * x142; auto x597 = x142 * x153;
      auto x598 = x153 * x550; auto x599 = 4. * x16; auto x600 = x513 * x599;
      auto x601 = x518 * x599; auto x602 = x153 * x553; auto x603 = x163 * x267;
      auto x604 = 2. * x561; auto x605 = 2. * x54; auto x606 = 3. * x142; auto x607 = 3. * UU;
      auto x608 = MXs * x607; auto x609 = 7. * x430; auto x610 = 7. * x432; auto x611 = 4. * x49;
      auto x612 = 4. * x50; auto x613 = 8. * x74; auto x614 = x611 + x612 + x613;
      auto x615 = 4. * x281; auto x616 = 4. * x283; auto x617 = 2. * x68;
      auto x618 = -x571 + x615 + x616 - x617; auto x619 = s * syFC6; auto x620 = x113 * x579;
      auto x621 = 2. * x620; auto x622 = x33 * x81; auto x623 = x33 * x84; auto x624 = MG * x52;
      auto x625 = x113 * x581; auto x626 = x442 * x566; auto x627 = 5. * x383;
      auto x628 = 8. * x158; auto x629 = 8. * MXs; auto x630 = MXs * x216;
      auto x631 = MG * x300 * x630; auto x632 = 3. * x21; auto x633 = 5. * x507;
      auto x634 = 3. * x72; auto x635 = x209 * x257; auto x636 = MG * x635; auto x637 = x389 * x636;
      auto x638 = x389 * x533; auto x639 = x142 * x209; auto x640 = MG * x376;
      auto x641 = x142 * x211; auto x642 = x158 * x209; auto x643 = x158 * x211;
      auto x644 = x209 * x552; auto x645 = x300 * x644; auto x646 = MG * x645;
      auto x647 = -x358 - x479; auto x648 = 2. * x562 + 2. * x563; auto x649 = 3. * x299;
      auto x650 = 3. * x56; auto x651 = 3. * x57; auto x652 = x300 * x59; auto x653 = x300 * x65;
      auto x654 = -x190 * x276 - x652 - x653; auto x655 = s * syFC7; auto x656 = 5. * x403;
      auto x657 = x11 * x84; auto x658 = 3. * x569; auto x659 = 3. * x570; auto x660 = 7. * x158;
      auto x661 = 7. * MXs; auto x662 = 5. * x142; auto x663 = 5. * x158; auto x664 = x200 * x267;
      auto x665 = -x100 * x375 + x182 * x375; auto x666 = 3. * x73; auto x667 = 3. * x76;
      auto x668 = 3. * x22 + 3. * x25 - x666 - x667; auto x669 = syFC5 * x87;
      auto x670 = x142 * x84; auto x671 = x550 * x81; auto x672 = x550 * x84;
      auto x673 = MUs * x145; auto x674 = x149 * x389; auto x675 = x151 * x389;
      auto x676 = 5. * x553 * x81; auto x677 = x140 * x657; auto x678 = 5. * x550;
      auto x679 = 3. * x147; auto x680 = x151 * x375; auto x681 = 4. * x114; auto x682 = 4. * x117;
      auto x683 = syFC6 * x87; auto x684 = x149 * x31; auto x685 = x151 * x31;
      auto x686 = MXs * x142; auto x687 = MXs * x475; auto x688 = MXs * x476;
      auto x689 = MUs * x198; auto x690 = MUs * x200; auto x691 = MXs * x158; auto x692 = MG * x389;
      auto x693 = x142 * x389; auto x694 = x158 * x389; auto x695 = x198 * x300;
      auto x696 = x200 * x300; auto x697 = x281 + x283; auto x698 = 3. * x29; auto x699 = MXs * x99;
      auto x700 = 15. * x119; auto x701 = 2. * x96; auto x702 = 8. * x329; auto x703 = 6. * x119;
      auto x704 = 4. * x88; auto x705 = x21 * x81; auto x706 = x21 * x84; auto x707 = x52 * x84;
      auto x708 = 12. * MXs; auto x709 = MUs * x149; auto x710 = MUs * x151; auto x711 = 20. * MXs;
      auto x712 = MUs * x99; auto x713 = 5. * SS; auto x714 = 5. * x16; auto x715 = x484 * x501;
      auto x716 = x485 * x715; auto x717 = MXs * Mqi; auto x718 = x245 * x717;
      auto x719 = x485 * x718; auto x720 = 4. * x516; auto x721 = x130 * x720 + x152 * x720;
      auto x722 = 12. * MUs; auto x723 = 6. * MUs;
      auto x724 = x122 - x124 - x156 + x163 + x472 - x473 + x94 - x99; auto x725 = syFC1 * x87;
      auto x726 = -x425; auto x727 = x209 * x568; auto x728 = x211 * x568; auto x729 = MXs * x728;
      auto x730 = MXs * x453; auto x731 = x300 * x446; auto x732 = -2. * x530;
      auto x733 = x300 * x445 + x732; auto x734 = x127 * x608 + 2. * x53 + 2. * x560;
      auto x735 = 2. * x534; auto x736 = x128 * x608 + 3. * x527 + 3. * x528 - x735;
      auto x737 = x114 * x607; auto x738 = x117 * x606; auto x739 = x117 * x607;
      auto x740 = 3. * x158; auto x741 = -x299 + x53 + x560 - x57; auto x742 = 2. * x476;
      auto x743 = x318 + x319 + 2. * x475 + x742; auto x744 = s * x478; auto x745 = x152 * x471;
      auto x746 = x11 * x153; auto x747 = x143 + x670; auto x748 = x478 * x87; auto x749 = 8. * x73;
      auto x750 = 2. * x11; auto x751 = 2. * x541;
      auto x752 = -x389 * x452 - x389 * x454 - 2. * x540; auto x753 = MX * SS;
      auto x754 = x104 * x753; auto x755 = x108 * x753;
      auto x756 = MXs * x754 + MXs * x755 + x530 + x534; auto x757 = x407 * x88;
      auto x758 = x173 * x51; auto x759 = x483 * x92; auto x760 = 2. * x38; auto x761 = SS * x114;
      auto x762 = SS * x117; auto x763 = x442 * x514; auto x764 = 2. * x577;
      auto x765 = x568 * x717; auto x766 = u * x407; auto x767 = 4. * x553; auto x768 = x217 * x389;
      auto x769 = 2. * MG; auto x770 = x643 * x769; auto x771 = 2. * x158;
      auto x772 = MXs * x761 + MXs * x762 - x114 * x771 - x117 * x771 - x249 * x389 - x250 * x389 +
                  x475 + x476 + x480 + x689 + x690;
      auto x773 = MG * x490; auto x774 = 2. * x59 + 2. * x65; auto x775 = x306 + x67;
      auto x776 = syFC3 * x87; auto x777 = MX * UU; auto x778 = x104 * x777;
      auto x779 = x108 * x777;
      auto x780 = UU * x117 - x100 - x101 - x172 + x182 + x185 + x190 + x192 + x200 - x250 - x445 -
                  x446 + x452 + x454 + x754 + x755 - x762 - x778 - x779;
      auto x781 = 2. * x452; auto x782 = 2. * x454; auto x783 = 3. * x762; auto x784 = 2. * x778;
      auto x785 = 2. * x779; auto x786 = 2. * x754; auto x787 = 2. * x755; auto x788 = 2. * x445;
      auto x789 = 2. * x446; auto x790 = 3. * x198 + 3. * x200 + x739 - 3. * x761 - x783 - x784 -
                                         x785 + x786 + x787 - x788 - x789;
      auto x791 = -2. * x100 + 2. * x182 - 3. * x250 + x546 + x781 + x782 + x790;
      auto x792 = x127 * x550; auto x793 = x128 * x550; auto x794 = MXs * x452;
      auto x795 = MXs * x454; auto x796 = u * x542; auto x797 = -x12 * x241; auto x798 = 3. * x59;
      auto x799 = 3. * x65; auto x800 = u * x559; auto x801 = x158 * x9;
      auto x802 = 2. * x217 - 6. * x326 - x786 - x787; auto x803 = 6. * x182; auto x804 = 6. * x185;
      auto x805 = 10. * x190; auto x806 = 10. * x192; auto x807 = 4. * x280; auto x808 = 4. * x282;
      auto x809 = 2. * x306; auto x810 = 2. * x67; auto x811 = 9. * MUs; auto x812 = 21. * MUs;
      auto x813 = x210 * x514 + x217 * x300; auto x814 = 4. * x11; auto x815 = 3. * MG;
      auto x816 = u * x619; auto x817 = 13. * MUs; auto x818 = 24. * MUs; auto x819 = x31 * x753;
      auto x820 = x113 * x529; auto x821 = MX * x116; auto x822 = x257 * x821;
      auto x823 = 4. * x442; auto x824 = 4. * x72; auto x825 = MX * x113; auto x826 = s * t;
      ,
      TR * TR * gs * gs * (Nc - 1) * (Nc + 1) *
          (-MG * x766 *
               (MXs * x820 + MXs * x822 + x113 * x819 + x116 * x819 - x158 * x271 - x158 * x823 +
                x209 * x824 + x211 * x483 + x211 * x824 + x270 * x501 + x270 * x514 - x271 * x553 +
                x287 * x642 + x287 * x643 - x386 + x442 * x501 - x459 - x462 - x500 - x504 -
                x553 * x823 - x645 + x763) +
           MG * x796 *
               (x113 * x526 + x116 * x526 - x158 * x821 - x158 * x825 - x209 * x550 + x209 * x553 -
                x211 * x550 + x211 * x553 + x383 * x821 + x383 * x825 - x384 - x404 * x753 + x502 +
                x505 + x532 - x537 + x550 * x821 + x552 * x821 - x553 * x821 - x630 + x635 - x639 -
                x641 + x642 + x643 - x644 - x820 - x822) +
           MUs * x725 * (UU * x114 + x147 - x151 + x198 - x249 - x761 + x780) +
           MUs * x776 *
               (3. * x145 - 3. * x151 - 3. * x249 - x607 * x84 + 3. * x657 + x679 + x737 + x791) +
           s * x344 *
               (-2. * x157 - 3. * x159 - x22 - 3. * x288 + 4. * x377 + x53 + x560 - 2. * x569 -
                2. * x570 + 2. * x574 - x649 - x651 + x666 + x667 - x671 - x672 + 4. * x673 + x674 +
                x675 - x676 - x677 + x743 + 2. * x745 - 2. * x746 + x747 + x764) -
           syFC1 * (-x117 * x208 - x127 * x186 + x127 * x201 + x168 + x170 - x180 - x181 -
                    x186 * x210 - x186 * x212 + x187 + x188 + x198 * x214 + x200 * x214 +
                    x201 * x210 + x201 * x212 - x202 - x203 - x204 - x205 + x207 * x87 +
                    x208 * x210 + x208 * x212 - x213 - x214 * x215 - x214 * x217 + x220) -
           4. * syFC2 *
               (-x159 * x287 - 2. * x278 + 2. * x281 + 2. * x283 - 4. * x284 - 4. * x285 +
                2. * x286 - x287 * x288 + x295 + x298 + x305 + x309 + x313 + x317 + x320) -
           syFC3 * (x101 * x191 - x114 * x224 - x117 * x224 + 3. * x168 + 3. * x170 + x175 + x179 -
                    3. * x180 - 3. * x181 - x193 + x194 * x210 - x195 - x196 + x198 * x225 +
                    x200 * x225 - 3. * x204 - 3. * x205 - x210 * x223 + x210 * x224 + x212 * x224 -
                    x215 * x225 + 3. * x218 + 3. * x219 - x221 * x222 + x226) -
           syFC5 *
               (x110 + x111 + x115 + x118 + x120 * x227 + x122 * x227 - x123 * x144 - x124 * x144 -
                x125 - x126 + x133 - x136 - x144 * x156 + x144 * x163 + x144 * x164 - x144 * x97 -
                x144 * x99 + 3. * x154 - 7. * x165 - x82 - x85 + 2. * x91 + 2. * x95) +
           syFC5 * (-x10 * x6 + x12 * x70 - x15 - x17 * x63 + x19 - 10. * x24 - x26 * x27 -
                    x30 * x32 - 5. * x35 + 2. * x39 + 2. * x42 + 2. * x45 + 2. * x47 + x49 * x6 -
                    x5 * x6 + x50 * x6 + x54 * x55 - x55 * x58 - x60 * x61 - x61 * x66 - x63 * x64 +
                    x68 * x69 + x70 * x71 + x74 * x75 + x75 * x77) -
           syFC6 * (-s * x237 + s * x240 + s * x243 - x215 * x239 + x236 * x26 - x238 * x28 +
                    8. * x24 + x242 * x28 + x244 * x28 + x246 * x247 + x251 + x30 * x62 + 4. * x35 -
                    x39 - x42 - x45 - x47 - x63 * x71) -
           syFC6 * (x100 * x233 + x101 * x233 - x114 * x234 - x117 * x234 + 10. * x168 +
                    10. * x170 - 5. * x180 - 5. * x181 - x190 * x232 - x192 * x232 + x198 * x235 +
                    x200 * x235 - 2. * x204 - 2. * x205 + x207 * x98 - 3. * x213 + x220 -
                    x222 * x229 - x222 * x230 + x222 * x231 + x226 - x228 * x25 + x228 * x76) +
           syFC6 * (x100 * x78 + x101 * x78 + x104 * x105 + x105 * x108 - x109 * x127 -
                    x109 * x128 - x110 - x111 - x115 - x118 - x120 * x121 - x121 * x122 +
                    x123 * x98 + x124 * x98 + x125 + x126 - x133 + x136 - x79 * x86 + x82 + x85 -
                    x91 + x93 - x95 + x97 * x98 + x98 * x99) +
           syFC7 * x79 *
               (AXG * x114 + AXG * x117 - AXG * x81 - AXG * x84 + AXG * x86 + AXG * x92 - x114 -
                x117 + x81 + x84 - x86) -
           syFC7 * x87 *
               (-MXs * x737 - MXs * x739 - x114 * x606 + x114 * x740 + x117 * x740 - 4. * x22 +
                x249 * x376 - 4. * x25 + x250 * x376 - 4. * x288 + 4. * x670 - x695 - x696 +
                4. * x73 - x738 + x741 + x743 + 4. * x76) -
           syFC7 * (-x100 * x259 - x101 * x259 + x119 * x273 - x16 * x273 * x88 + x190 * x272 +
                    x192 * x272 + x247 * x265 + x251 - 2. * x252 - 2. * x253 + 2. * x255 +
                    2. * x258 - 8. * x261 - 8. * x262 + x264 * x28 + x266 * x267 + x268 * x271 -
                    x275 * x28 + x276 * x277) +
           syFC7 * (MXs * x146 + MXs * x148 - x120 * x138 - x122 * x138 + x123 * x137 +
                    x124 * x137 - 2. * x132 + 6. * x135 + x137 * x156 + x137 * x157 - x137 * x163 -
                    x137 * x164 + x137 * x97 + x137 * x99 - x139 * x140 - x140 * x141 -
                    x143 * x144 + x144 * x159 - x149 * x150 - x150 * x151 - 2. * x154 +
                    x160 * x161 + x161 * x162 + 6. * x165 + x93) +
           syFC7 * (-x100 * x183 - x101 * x183 - 9. * x168 - 9. * x170 - 4. * x171 - x175 - x179 +
                    4. * x180 + 4. * x181 + x182 * x183 + x183 * x185 + 4. * x184 - x187 - x188 -
                    8. * x189 + x190 * x191 + x193 + x195 + x196 + 20. * x197 - x198 * x199 -
                    x199 * x200 + 9. * x202 + 9. * x203 + x206) -
           u * x465 *
               (MG * x763 - x237 + x292 + x294 + x31 * x761 + x31 * x762 + x408 + x411 + x434 +
                x436 - x463 + x571 - x611 - x612 - x613 - x615 - x616 + x617 + x652 + x653 + x687 +
                x688 - x758 - x759 + x760) +
           u * x655 *
               (MUs * x805 + MUs * x806 - x100 * x812 - x101 * x812 - x117 * x660 - x172 * x267 +
                x182 * x811 + x185 * x811 - 4. * x207 + 6. * x22 + 6. * x25 + 5. * x299 - x53 -
                x560 + 5. * x57 + x664 - 18. * x73 + x738 - x742 - 18. * x76 + x774 + x797 - x807 -
                x808 + x809 + x810) -
           u * x683 *
               (-10. * x100 - 10. * x101 - x231 + x246 - 7. * x250 + x443 + x444 + 6. * x452 +
                x453 + 6. * x454 - x727 - x728 + x790 + x803 + x804) -
           u * x725 * (-x176 - x254 - x260 + x265 + x326 + x327 + x780 - x801) -
           u * x776 *
               (-2. * x101 - x173 - x177 + 2. * x192 + x221 + 3. * x265 + 3. * x326 - x335 + x791) +
           x2 * x748 * x88 * (-SS + UU - x11 + x16) +
           x344 *
               (s * x97 + s * x99 + x120 * x6 + x122 * x6 - x123 * x69 - x124 * x69 + x137 * x254 +
                x137 * x466 + x139 + x141 - x146 - x148 - x156 * x69 + x160 + x162 + x467 + x468 -
                x469 - x470 + x472 * x69 - x473 * x69 - x69 * x90 - x69 * x94) -
           x344 * (x120 * x267 + x122 * x267 - x158 * x263 * x28 + x17 * x239 + x26 +
                   2. * x265 * x28 + x277 * x300 + x28 * x326 + x28 * x327 - x28 * x335 -
                   4. * x321 - 4. * x322 - 4. * x323 + x325 - 4. * x328 - 4. * x329 + x331 + x333 +
                   x334 * x86 + x334 * x92 + x340 + x343 + x41 * x69 + x46 * x69 * x9) +
           x352 * (-MUs * x231 * x51 - MUs * x289 * x64 + MUs * x479 - x241 * x280 - x241 * x282 +
                   x354 * x466 - x375 * x428 + x482) -
           x352 * (-x100 * x247 - x101 * x247 - x120 * x345 + x120 * x349 - x122 * x345 +
                   x122 * x349 + x123 * x51 + x124 * x51 + x18 * x71 + x190 * x239 + x192 * x239 +
                   x20 * x90 + x20 * x94 + x236 * x350 - x248 - x347 - x348 + x351) -
           x352 * (MUs * x314 + MUs * x315 - MUs * x358 + x12 * x355 - x236 * x284 - x236 * x285 -
                   x254 * x354 + x265 * x289 + x300 * x361 + x300 * x362 + x353 * x81 + x353 * x84 +
                   x355 * x71 + x356 + x357 - x359 * x86 - x359 * x92 - x360 * x81 - x360 * x84 +
                   x364 + x370) -
           x352 *
               (MUs * x336 + MUs * x337 + MXs * x316 + MXs * x330 + MXs * x333 + x145 * x354 +
                x147 * x354 - x149 * x289 - x151 * x289 + x236 * x380 + x267 * x286 - x267 * x328 -
                x267 * x329 + x267 * x332 - x278 * x376 + x286 * x287 - x323 * x375 - x323 * x376 -
                x371 - x372 - x373 - x374 - 3. * x377 * x51 + x378 + x379) +
           x407 * (-8. * x122 * x20 + x302 * x486 + x302 * x489 + x347 + x348 - x483 * x486 +
                   x487 * x51 + x494) -
           x407 * (-x210 * x381 + x210 * x382 - x212 * x381 + x212 * x382 - x300 * x388 -
                   x300 * x390 + x385 * x51 + x386 * x387 - x388 * x389 - x389 * x390 + x406) -
           x407 *
               (-x100 * x272 - x101 * x272 - x166 * x217 + x166 * x443 + x166 * x444 - x166 * x445 -
                x166 * x446 + x239 * x452 + x239 * x453 + x239 * x454 + x270 * x441 -
                x270 * x448 * x88 + x270 * x451 + x351 + x441 * x442 + x442 * x451 - x449 * x88) +
           x407 *
               (-MUs * x498 - MUs * x499 + x104 * x508 + x108 * x508 - x31 * x503 - x31 * x506 -
                x387 * x500 - x387 * x504 + x396 * x495 + x396 * x497 + x401 * x495 + x401 * x497 +
                x456 + x457 + x458 + x461 + x464 + x481 - x495 * x496 - x496 * x497 + x509) -
           x407 * (MXs * x431 + MXs * x433 + x104 * x420 + x108 * x420 - x409 - x410 - x412 - x413 -
                   x414 - x415 + x417 + x419 - x421 - x422 + x423 + x424 + x426 + x427 + x429 -
                   x435 - x437 + x439 + x440) +
           x465 * x88 *
               (-x416 - x418 + x425 - x442 * x773 + x449 - x749 - 8. * x76 - x764 + x772 + x774 +
                x775) +
           x465 * (-x417 - x419 - x427 - x429 + x482 + x509) +
           x465 * (-x130 * x512 + x130 * x517 - x152 * x512 + x152 * x517 + x494 + x510 * x511 -
                   x513 * x514 + x513 * x515 + x515 * x518) -
           x465 * (x300 * x315 + 2. * x356 + x370 + x406 - x456 - x457 - x458 - x461 - x464) +
           x478 * (MXs * x467 + MXs * x468 - MXs * x469 - MXs * x470 + x114 * x474 - x114 * x477 +
                   x117 * x474 - x117 * x477 + x122 * x28 - x168 - x170 - x171 + x180 + x181 +
                   x184 - x189 + x197 + x198 * x222 + x200 * x222 + x206 - x218 - x219 - x28 * x99 -
                   x475 * x87 - x476 * x87) +
           x542 * (-x119 * x128 - x127 * x89 - x128 * x89 + x128 * x96 - 2. * x329 + x340 -
                   x389 * x99 + x445 * x88 + x446 * x88 + x450 * x86 + x450 * x92 - x543 -
                   2. * x544 - x545 - x546 * x88 + x547 + x548 + x549 + x555) +
           x542 * (MXs * x527 + MXs * x528 - MXs * x533 + MXs * x536 + MXs * x539 - MXs * x540 -
                   MXs * x541 - x23 + x311 + x312 + x38 + x44 - x519 - x520 - x521 + x522 + x523 -
                   x524 - x525 - x531 - x535 + x537 * x538 - x60 - x66 + x74 + x77) -
           x542 * (MXs * x636 - MXs * x680 + MXs * x689 + MXs * x690 + x114 * x686 - x114 * x691 +
                   x117 * x686 - x117 * x691 + x288 * x376 + x298 + x365 + x366 + x376 * x377 -
                   x376 * x670 + x376 * x673 - x384 * x538 - x538 * x639 - x538 * x641 +
                   x538 * x642 + x538 * x643 - x622 - x623 - x684 - x685 - x687 - x688) +
           x559 * x88 *
               (5. * x100 + 5. * x101 + 5. * x172 + 5. * x176 - 5. * x192 - x211 * x773 - x229 -
                x230 + 6. * x254 + 6. * x260 - 6. * x265 - 6. * x327 - x493 + 2. * x728 - x781 -
                x782 + x784 + x785 + x788 + x789 + 6. * x801 + x802) -
           x559 *
               (MXs * x695 + MXs * x696 + x114 * x693 - x114 * x694 + x117 * x693 - x117 * x694 +
                3. * x297 + x309 - x389 * x475 - x389 * x476 - x54 - x561 + x564 + x58 - x631 +
                x637 + x638 - x639 * x692 - x641 * x692 + x642 * x692 + x643 * x692 - x646 + x697) +
           x559 * (x101 * x241 - x190 * x241 - x192 * x241 - 3. * x23 - x240 - x243 + x300 * x556 +
                   x300 * x557 + 3. * x311 + 3. * x312 + 3. * x38 + x389 * x527 + x389 * x528 -
                   x389 * x530 - x389 * x534 - x389 * x540 - x389 * x541 + 3. * x44 - 3. * x520 -
                   3. * x521 + x558 + 3. * x74 + 3. * x77) +
           2. * x559 *
               (-x124 * x375 - x131 * x158 + x131 * x383 - x153 * x158 + x153 * x383 - x153 * x501 +
                x153 * x552 + x158 * x486 + x375 * x94 - x485 * x765 - x486 * x552 + x488 * x715 +
                x488 * x718 - x488 * x765 - x489 * x552 + x555 - x594 - x595 + x596 + x597 + x598 -
                x602 + x716 + x719) +
           x587 * (-x270 * x567 - x389 * x569 - x389 * x570 - x442 * x567 + 6. * x561 + 6. * x562 +
                   6. * x563 - 6. * x564 - x565 * x86 - x565 * x92 + x571 + x572 * x573 +
                   x573 * x576 + x574 * x575 + x575 * x577 + x582 + x584 + x586) +
           x587 * (x122 * x593 - x123 * x27 - x124 * x27 - x131 * x573 - x153 * x573 + x164 * x267 -
                   x236 * x99 + 11. * x286 - x331 + 11. * x332 + x549 + 2. * x551 - 6. * x554 -
                   3. * x588 - x591 - x592 - 4. * x594 - 4. * x595 + 2. * x596 + 2. * x597 +
                   2. * x598 + x600 + x601 - 6. * x602 + x603) -
           x587 * (SS * x324 + x114 * x701 - x114 * x703 + x117 * x701 - x117 * x703 - x120 * x593 -
                   x166 * x698 - x17 * x272 + x198 * x704 + x200 * x704 + x236 * x97 + x239 * x40 -
                   x256 * x607 - 6. * x266 + 11. * x321 + 11. * x322 + 8. * x328 - 6. * x336 -
                   6. * x337 + 7. * x544 - x681 * x89 - x682 * x89 + 7. * x699 - x700 * x86 -
                   x700 * x92 + x702) -
           x587 * (-x143 * x708 + x159 * x711 + 15. * x284 + 15. * x285 + x288 * x711 + 2. * x314 +
                   2. * x315 + 2. * x316 - 4. * x361 - 4. * x362 + 5. * x365 + 5. * x366 +
                   x377 * x629 - x621 - 5. * x622 - 5. * x623 + x625 + x629 * x673 - x670 * x708 -
                   4. * x684 - 4. * x685 - 7. * x705 - 7. * x706 + 3. * x707 - x708 * x709 -
                   x708 * x710) +
           x619 * x88 *
               (22. * x100 + 22. * x101 + 6. * x172 + 6. * x176 + x210 * x447 - 5. * x246 +
                14. * x260 + 4. * x445 + 4. * x446 - 9. * x452 - 5. * x453 - 9. * x454 + x727 +
                x728 + x778 + x779 + x783 + x802 - x803 - x804 - x805 - x806) +
           x619 * (-4. * x10 - x158 * x609 - x158 * x610 + x310 + x430 * x583 + x430 * x606 +
                   x430 * x608 - 4. * x431 + x432 * x583 + x432 * x606 + x432 * x608 - 4. * x433 -
                   4. * x5 - x553 * x609 - x553 * x610 + x604 + x605 + x614 + x618) +
           x619 *
               (-x104 * x578 + x104 * x633 - x108 * x578 + x108 * x633 + x210 * x632 - x210 * x634 +
                x212 * x632 - x212 * x634 + x305 + 5. * x463 - x576 * x627 + x631 - x637 - x638 +
                x639 * x640 + x640 * x641 - x640 * x642 - x640 * x643 + x646 + x647 + x648) +
           x619 * (x149 * x241 + x151 * x241 - x209 * x624 - x211 * x624 - x271 * x566 -
                   x287 * x569 - x287 * x570 + x298 - x367 + x388 + x390 + 5. * x460 - x572 * x627 +
                   x572 * x628 + x574 * x629 + x576 * x628 + x577 * x629 + x582 + x621 + 4. * x622 +
                   4. * x623 - x625 - 4. * x626) -
           x619 * (-14. * MUs * x122 + Mqi * x376 * x485 * x568 + x123 * x236 + x124 * x236 -
                   x131 * x627 + x131 * x628 - x153 * x627 + x153 * x628 + x389 * x487 +
                   x513 * x713 - x513 * x714 + x518 * x713 - x518 * x714 - 4. * x551 + 8. * x554 +
                   3. * x594 + 3. * x595 - 4. * x596 - 4. * x597 - 4. * x598 + 8. * x602 +
                   6. * x712 - 2. * x716 - 7. * x719 + x721) +
           x655 * (MUs * x649 + MUs * x651 - 3. * x10 - 6. * x23 + 3. * x49 - 3. * x5 + 3. * x50 -
                   6. * x520 + x54 + x558 + x561 + x562 + x563 + x618 + x650 * x86 + x650 * x92 +
                   x654 + 6. * x74 + 6. * x77) +
           x655 * (-MXs * x658 - MXs * x659 - x113 * x656 - x116 * x656 - x241 * x657 -
                   3. * x270 * x566 + 2. * x296 + 2. * x297 + x368 + x369 + x572 * x660 +
                   x574 * x661 + x576 * x660 + x577 * x661 + x580 + x584 + x620 + 3. * x622 +
                   3. * x623 - 3. * x626 + x647) -
           x655 * (-x120 * x722 - x122 * x722 + x123 * x723 + x124 * x723 - x131 * x606 -
                   x131 * x608 + x131 * x660 - x153 * x608 + x267 * x97 - 4. * x337 - 20. * x380 +
                   x549 + 7. * x554 + x591 + x592 + 2. * x594 - x600 - x601 + 7. * x602 - x603 +
                   4. * x699 + x702 + 4. * x712 + x721) +
           x669 * (-x101 * x375 + x114 * x662 - x114 * x663 + x117 * x662 - x117 * x663 - x174 -
                   x178 + x185 * x375 + x190 * x300 + x192 * x300 + x198 * x267 - 4. * x475 -
                   4. * x476 - 3. * x53 - 3. * x560 + x649 + x651 + x664 + x665 + x668) +
           x669 * (x114 * x678 + x117 * x678 - x140 * x249 - x140 * x250 - 6. * x143 + x145 * x389 +
                   x147 * x389 - x149 * x267 - x151 * x267 + 4. * x157 + 6. * x159 + 6. * x288 +
                   2. * x377 + 3. * x574 + 3. * x577 - x658 - x659 - 6. * x670 - 5. * x671 -
                   5. * x672 + 2. * x673 - x674 - x675 + x676 + x677) -
           x683 *
               (-MXs * x246 + MXs * x727 - 5. * x22 - x376 * x452 - x376 * x454 + 3. * x476 - x480 -
                3. * x540 - 3. * x541 + x726 + x729 + 5. * x73 - x730 + x731 + x733 + x734 + x736) +
           x683 * (MXs * x679 + x142 * x681 + x142 * x682 - 5. * x143 - x149 * x375 - x158 * x681 -
                   x158 * x682 + 5. * x159 + x198 * x375 + x200 * x375 - x249 * x287 - x250 * x287 +
                   5. * x288 + x377 - 3. * x475 + x550 * x681 + x550 * x682 - x569 - x570 + x574 +
                   x577 - 5. * x670 + x673 - x680) -
           x725 * (x120 - x123 + x724 + x90 - x97) -
           x744 * (x10 + 2. * x23 - x300 * x73 - x300 * x76 + x308 + x313 - x49 + x5 - x50 + x519 +
                   2. * x520 + x521 - x522 - x523 + x524 + x525 + x60 - x605 + x66 + x697) +
           x744 * (x143 * x389 - x159 * x389 - x288 * x389 - x296 - x297 - x301 - x303 - x304 +
                   x307 + x320 - x365 - x366 - x377 * x389 + x389 * x670 - x389 * x673 + x586 +
                   x604 + x622 + x623 + x648 + x684 + x685) -
           x744 *
               (-x278 + x284 + x285 - x286 + x317 - x323 + x329 + x330 - x332 - x337 + x338 + x339 +
                x343 - x361 - x362 + x543 + x545 - x547 - x548 + x588 - x705 - x706 + x707) -
           x748 * (-x157 - x159 + x164 - x288 - x377 + x569 + x570 - x574 - x577 - x673 + x709 +
                   x710 + x724 + x745 - x746 + x747) +
           x757 * (-MG * x537 - x215 * x389 + x246 * x287 + x287 * x453 - x385 + x503 + x506 -
                   x533 - x636 + x642 * x769 - x768 + x770 + x772) -
           x757 *
               (-x306 + x416 + x418 + x430 * x447 - x430 * x490 - x430 * x750 + x432 * x447 -
                x432 * x490 - x432 * x750 + x536 + x539 - x67 + x726 + x749 - x751 + x752 + x756) +
           x766 * (-x104 * x491 - x108 * x491 + x209 * x455 + x211 * x455 + x295 + x430 * x501 -
                   x430 * x515 - x430 * x767 + x432 * x501 - x432 * x767 - x434 - x436 +
                   x495 * x507 + x497 * x507 + x498 + x499 + x618 + x631 + x758 + x759) +
           x766 * (-x192 * x276 + x237 + x31 * x754 + x31 * x755 - 2. * x311 - x408 - x411 +
                   x430 * x514 - x431 + x432 * x514 - x432 * x515 - x433 - 2. * x44 + 4. * x522 +
                   4. * x523 + x531 + x535 + x614 + x654 - x760) +
           3. * x776 * (-x120 - x122 + x123 + x124 - x472 + x473 - x90 - x94 + x97 + x99) +
           x796 * (x190 * x375 + x22 + x25 - x37 - x43 - x527 - x528 - x536 - x539 + x540 + x541 -
                   x556 - x557 + x59 + x65 + x665 - x73 + x756 - x76 - x792 - x793 + x794 + x795) -
           x800 * (-x280 - x282 + x389 * x445 + x389 * x446 - x389 * x754 - x389 * x755 +
                   2. * x527 - 2. * x533 - 2. * x636 + x639 * x769 + x641 * x769 + x733 + x741 +
                   x752 - x770 + x775 + 2. * x792 + 2. * x793 + x813) +
           x800 *
               (-x100 * x236 - x101 * x236 - x172 * x236 + x182 * x236 + x185 * x236 + x190 * x236 +
                x192 * x236 - 8. * x207 + x242 + x244 + x264 - x275 + x31 * x698 - 3. * x37 -
                3. * x41 - 3. * x43 - 2. * x528 + x668 - x731 + x735 + x751 + x797 + x798 + x799) -
           x816 * (-6. * x299 - x430 * x589 + x430 * x599 - x430 * x814 - x432 * x589 +
                   x432 * x599 - x432 * x814 - 2. * x506 - 6. * x57 + x639 * x815 + x641 * x815 -
                   x642 * x815 - x643 * x815 + 3. * x729 - 7. * x730 + x734 + x768 - 7. * x794 +
                   x807 + x808 - x809 - x810 + x813) -
           x816 * (7. * MUs * x172 + x100 * x818 + x101 * x818 - x182 * x722 - x185 * x722 -
                   x190 * x817 - x192 * x817 + 7. * x207 - 8. * x22 + x238 - 8. * x25 +
                   x267 * x445 + x267 * x446 + x37 + x43 - 7. * x540 - 7. * x541 + 20. * x73 +
                   x732 + x736 + 20. * x76 - 7. * x795 - x798 - x799)) *
          Denom(768 * (-MGs * x28 + MGs * x826 + MGs * x87 + MXs * x826 - s * pow(t, 2) - t * x87) *
                pow(Pi, 2)));
  // std::cout << (iLp * iRG * LGp * R + iLG * iRp * L * RGp) << std::endl;

  if (color_flow) {
    return -ret.real();
  }
  return ret.real();
}

ComplexType ME_us_gGX_qQQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;
  for (int invert_color_flow = 0; invert_color_flow < 2; invert_color_flow++) {
    for (int itq = 0; itq < 6; itq++) {
      for (int itsq = 0; itsq < 2; itsq++) {
        ret += ME_us_gGX_qQQ_single(pIEPS, sc, uc, axial, itsq, itq, invert_color_flow, Q2, P1K1,
                                    params);
      }
    }
  }
  return ret;
}

ComplexType ME_us_gGX_Qqq_single(POLE pIEPS, bool sc, bool uc, bool axial, int itsq, int itq,
                                 int color_flow, double Q2, double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // MG = -MG;
  double SS = sc, UU = uc, AXG = axial;
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){
  int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
  int iq = itq;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType iLp = conj(params->CHSQq[ch][isq][iq].R);
  ComplexType iRp = conj(params->CHSQq[ch][isq][iq].L);

  auto MQi = params->mSQ[isq];
  auto MQis = pow2(MQi);

  auto Mqi = params->mq[itq];
  // std::cout << "mq: " << Mqi << std::endl;
  auto Mqis = pow2(Mqi);

  ComplexType iLG = (params->GLSQq[isq][iq].L);
  ComplexType iRG = (params->GLSQq[isq][iq].R);
  ComplexType LGp = conj(params->GLSQq[sq][q].R);
  ComplexType RGp = conj(params->GLSQq[sq][q].L);

  if (color_flow) {
    auto tiLp = conj(iRp);
    auto tiRp = conj(iLp);
    auto tiLG = conj(iRG);
    auto tiRG = conj(iLG);

    iLp = tiLp;
    iRp = tiRp;
    iLG = tiLG;
    iRG = tiRG;
  }

  // std::cout << (iLp * iRG * LGp * R + iLG * iRp * L * RGp) << std::endl;
  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 B0(MXs, Mqis, MQis)
#define syFC2 C0(0, t, MXs, Mqis, Mqis, MQis)
#define syFC3 C00(0, MXs, t, Mqis, Mqis, MQis)
#define syFC4 C1(0, MXs, t, Mqis, Mqis, MQis)
#define syFC5 C11(0, MXs, t, Mqis, Mqis, MQis)
#define syFC6 C12(0, MXs, t, Mqis, Mqis, MQis)
#define syFC7 C2(0, MXs, t, Mqis, Mqis, MQis)
#define syFC8 C22(0, MXs, t, Mqis, Mqis, MQis)

  _EPS0_(
      ret, auto x0 = L * RGp; auto x1 = iLG * iRp; auto x2 = x0 * x1; auto x3 = pow(MXs, 2);
      auto x4 = UU * x3; auto x5 = x2 * x4; auto x6 = s * x5; auto x7 = iLp * iRG;
      auto x8 = LGp * R; auto x9 = x7 * x8; auto x10 = s * x9; auto x11 = AXG * UU;
      auto x12 = x11 * x3; auto x13 = 3. * x12; auto x14 = x10 * x4; auto x15 = x11 * x2;
      auto x16 = pow(u, 3); auto x17 = 4. * x16; auto x18 = x11 * x9; auto x19 = AXG * SS;
      auto x20 = 2. * x16; auto x21 = x19 * x20; auto x22 = s * x2; auto x23 = pow(MUs, 2);
      auto x24 = SS * x23; auto x25 = 2. * x24; auto x26 = MXs * s; auto x27 = MUs * UU;
      auto x28 = x26 * x27; auto x29 = pow(u, 2); auto x30 = x19 * x29; auto x31 = x2 * x30;
      auto x32 = MXs * x31; auto x33 = x30 * x9; auto x34 = MXs * x33; auto x35 = x11 * x29;
      auto x36 = x2 * x35; auto x37 = 4. * MUs; auto x38 = u * x2; auto x39 = pow(MX, 4);
      auto x40 = x11 * x39; auto x41 = x38 * x40; auto x42 = 4. * x41; auto x43 = x35 * x9;
      auto x44 = MXs * x36; auto x45 = MXs * x43; auto x46 = 3. * x2; auto x47 = MUs * x11;
      auto x48 = x26 * x47; auto x49 = 3. * x9; auto x50 = x19 * x2; auto x51 = 2. * MUs;
      auto x52 = x26 * x51; auto x53 = x19 * x9; auto x54 = SS * x29; auto x55 = x2 * x54;
      auto x56 = MXs * x55; auto x57 = x54 * x9; auto x58 = MXs * x57;
      auto x59 = 2. * x56 + 2. * x58; auto x60 = x51 * x55; auto x61 = x51 * x57;
      auto x62 = -x60 - x61; auto x63 = x31 * x51; auto x64 = x33 * x51; auto x65 = -x63 - x64;
      auto x66 = x10 * x25 + x14 - x15 * x17 - x17 * x18 + x2 * x21 + x2 * x28 + x21 * x9 +
                 x22 * x25 + x28 * x9 - 4. * x32 - 4. * x34 + x36 * x37 + x37 * x43 + x42 +
                 8. * x44 + 8. * x45 - x46 * x48 - x48 * x49 + x50 * x52 + x52 * x53 + x59 + x62 +
                 x65;
      auto x67 = pow(MXs, 3); auto x68 = UU * x67; auto x69 = x2 * x68; auto x70 = 5. * s;
      auto x71 = x68 * x9; auto x72 = pow(u, 4); auto x73 = x15 * x72; auto x74 = x19 * x72;
      auto x75 = x2 * x74; auto x76 = 10. * MUs; auto x77 = UU * x23; auto x78 = x26 * x77;
      auto x79 = x2 * x78; auto x80 = x78 * x9; auto x81 = MUs * x3; auto x82 = 2. * SS;
      auto x83 = x81 * x82; auto x84 = SS * x2; auto x85 = 2. * x84; auto x86 = x23 * x26;
      auto x87 = x85 * x86; auto x88 = SS * x9; auto x89 = x23 * x88; auto x90 = x26 * x89;
      auto x91 = 2. * x90; auto x92 = x11 * x67; auto x93 = 5. * x92; auto x94 = UU * x39;
      auto x95 = x9 * x94; auto x96 = MUs * x95; auto x97 = 6. * MUs; auto x98 = x19 * x81;
      auto x99 = 4. * x98; auto x100 = x19 * x23; auto x101 = x100 * x26; auto x102 = x101 * x9;
      auto x103 = MUs * x39; auto x104 = x103 * x19; auto x105 = x18 * x86;
      auto x106 = 10. * x11 * x81; auto x107 = pow(s, 2); auto x108 = 4. * x107;
      auto x109 = pow(s, 3); auto x110 = u * x109; auto x111 = x110 * x2; auto x112 = x110 * x9;
      auto x113 = MUs * x85; auto x114 = UU * x111; auto x115 = x51 * x88; auto x116 = UU * x112;
      auto x117 = UU * x29; auto x118 = x117 * x2; auto x119 = x108 * x118; auto x120 = iLG * iLp;
      auto x121 = Mqi * x0; auto x122 = x120 * x121; auto x123 = MX * x109; auto x124 = UU * x123;
      auto x125 = iRG * iRp; auto x126 = Mqi * x125; auto x127 = x126 * x8; auto x128 = x11 * x123;
      auto x129 = x11 * x16; auto x130 = x129 * x22; auto x131 = x10 * x129; auto x132 = x118 * x26;
      auto x133 = x117 * x9; auto x134 = 4. * x133; auto x135 = x107 * x113; auto x136 = MXs * x135;
      auto x137 = x107 * x115; auto x138 = MXs * x137; auto x139 = x17 * x19; auto x140 = x26 * x55;
      auto x141 = MXs * x37; auto x142 = UU * x107; auto x143 = x142 * x9; auto x144 = MG * x0;
      auto x145 = x144 * x7; auto x146 = MX * x117; auto x147 = s * x146; auto x148 = MG * x1;
      auto x149 = x148 * x8; auto x150 = x26 * x31; auto x151 = MUs * MXs; auto x152 = MXs * x51;
      auto x153 = x152 * x50; auto x154 = x152 * x53; auto x155 = MX * x145; auto x156 = 2. * s;
      auto x157 = x156 * x54; auto x158 = MX * x149; auto x159 = x26 * x36; auto x160 = x26 * x43;
      auto x161 = s * x30; auto x162 = 4. * x161; auto x163 = s * x35; auto x164 = 9. * x163;
      auto x165 = UU * x16; auto x166 = x165 * x22; auto x167 = x10 * x165; auto x168 = x166 + x167;
      auto x169 = u * x9; auto x170 = SS * x38; auto x171 = 2. * x107; auto x172 = SS * x169;
      auto x173 = x107 * x11; auto x174 = MX * x142; auto x175 = x145 * x174;
      auto x176 = x149 * x174; auto x177 = 2. * x30; auto x178 = 5. * x35; auto x179 = x155 * x173;
      auto x180 = x158 * x173;
      auto x181 = s * x118 + s * x133 - x10 * x177 + x10 * x178 + x135 + x137 + x142 * x169 +
                  x142 * x38 - x156 * x55 - x156 * x57 + x169 * x173 - x170 * x171 - x171 * x172 +
                  x173 * x38 + x175 + x176 - x177 * x22 + x178 * x22 - x179 - x180;
      auto x182 = x123 * x145; auto x183 = x123 * x149; auto x184 = 8. * x107; auto x185 = SS * u;
      auto x186 = MX * x185; auto x187 = x145 * x186; auto x188 = 3. * u; auto x189 = 7. * u;
      auto x190 = u * x19; auto x191 = x155 * x190; auto x192 = x158 * x190; auto x193 = x4 * x9;
      auto x194 = SS * x16; auto x195 = x194 * x22; auto x196 = 2. * x95; auto x197 = x107 * x23;
      auto x198 = x197 * x84; auto x199 = x197 * x88; auto x200 = 2. * Mqi; auto x201 = MUs * SS;
      auto x202 = x0 * x120; auto x203 = MX * x202; auto x204 = x201 * x203;
      auto x205 = x200 * x204; auto x206 = x125 * x8; auto x207 = MX * x206;
      auto x208 = x201 * x207; auto x209 = x107 * x208; auto x210 = 3. * x122;
      auto x211 = MUs * x174; auto x212 = 3. * MXs; auto x213 = x174 * x212; auto x214 = 3. * x127;
      auto x215 = MX * x122; auto x216 = x173 * x215; auto x217 = 3. * MUs; auto x218 = MX * x127;
      auto x219 = x173 * x218; auto x220 = MUs * x19; auto x221 = x203 * x220;
      auto x222 = MUs * x207; auto x223 = x19 * x222; auto x224 = x107 * x223;
      auto x225 = pow(MUs, 3); auto x226 = SS * x225; auto x227 = x22 * x226;
      auto x228 = x10 * x226; auto x229 = pow(MX, 3); auto x230 = x229 * x27;
      auto x231 = x125 * x144; auto x232 = x16 * x19; auto x233 = x200 * x232;
      auto x234 = x200 * x203; auto x235 = x120 * x8; auto x236 = MG * x235;
      auto x237 = x200 * x207; auto x238 = x121 * x125; auto x239 = MG * x238;
      auto x240 = x11 * x17; auto x241 = Mqi * x236; auto x242 = x201 * x229;
      auto x243 = 4. * s * x242; auto x244 = x19 * x229; auto x245 = x122 * x244;
      auto x246 = s * x37; auto x247 = x127 * x244; auto x248 = x11 * x229; auto x249 = x127 * x248;
      auto x250 = 7. * MUs; auto x251 = x103 * x84; auto x252 = 2. * x251; auto x253 = x103 * x88;
      auto x254 = 2. * x253; auto x255 = x104 * x2; auto x256 = x81 * x84; auto x257 = x81 * x88;
      auto x258 = x104 * x9; auto x259 = MXs * x23; auto x260 = x259 * x84; auto x261 = x40 * x9;
      auto x262 = -x259 * x88 - x260 + x261 * x51; auto x263 = -x256 - x257 - x258 + x262;
      auto x264 = x170 * x3; auto x265 = x172 * x3; auto x266 = x19 * x39; auto x267 = x169 * x266;
      auto x268 = x264 + x265 + x267; auto x269 = x170 * x23; auto x270 = x172 * x23;
      auto x271 = pow(MX, 6); auto x272 = x11 * x271; auto x273 = x2 * x272; auto x274 = 2. * x273;
      auto x275 = x272 * x9; auto x276 = 2. * x275; auto x277 = x172 * x39; auto x278 = x12 * x38;
      auto x279 = x12 * x169; auto x280 = x2 * x40; auto x281 = x280 * x51; auto x282 = 2. * MXs;
      auto x283 = x280 * x282; auto x284 = x261 * x282; auto x285 = x169 * x40;
      auto x286 = x11 * x38; auto x287 = x11 * x169; auto x288 = x19 * x38; auto x289 = x169 * x19;
      auto x290 = x152 * x288 + x152 * x289;
      auto x291 = -x141 * x286 - x141 * x287 + x269 + x270 - x274 - x276 - 2. * x277 - 4. * x278 -
                  4. * x279 + x281 + x283 + x284 + 2. * x285 + x290;
      auto x292 = x170 * x39; auto x293 = 2. * x39; auto x294 = x19 * x293; auto x295 = x294 * x38;
      auto x296 = 2. * syFC3; auto x297 = x3 * x55; auto x298 = x3 * x57; auto x299 = pow(MU, 4);
      auto x300 = 2. * x280; auto x301 = 2. * x261; auto x302 = 2. * x3; auto x303 = x19 * x299;
      auto x304 = 2. * x23; auto x305 = -x302 * x31 - x302 * x33; auto x306 = x23 * x31;
      auto x307 = x23 * x33; auto x308 = -x306 - x307; auto x309 = 2. * syFC4;
      auto x310 = SS * x299; auto x311 = x2 * x310; auto x312 = x310 * x9; auto x313 = x2 * x266;
      auto x314 = x266 * x9; auto x315 = MUs * x311; auto x316 = MUs * x312; auto x317 = MUs * x264;
      auto x318 = MXs * x269; auto x319 = MUs * x265; auto x320 = MXs * x270; auto x321 = 4. * x81;
      auto x322 = x19 * x259; auto x323 = 2. * x38; auto x324 = 2. * x169;
      auto x325 = x170 * x225 + x172 * x225 + x322 * x323 + x322 * x324; auto x326 = SS * x67;
      auto x327 = x326 * x9; auto x328 = 2. * pow(MX, 8); auto x329 = 4. * x273;
      auto x330 = 4. * x275; auto x331 = -MXs * x329 - MXs * x330 + u * x327 + x15 * x328 +
                                         x18 * x328 + x323 * x98 + x324 * x98;
      auto x332 = MXs * x16; auto x333 = 4. * x15; auto x334 = 4. * x29; auto x335 = 4. * x18;
      auto x336 = 2. * x50; auto x337 = 2. * x53; auto x338 = 6. * x43; auto x339 = 8. * MUs;
      auto x340 = MUs * x16; auto x341 = 2. * x18;
      auto x342 = -x332 * x84 - x332 * x88 - x333 * x340 - x335 * x340 + x336 * x340 + x337 * x340 +
                  x340 * x84 + x340 * x88 + x341 * x72 - x74 * x9;
      auto x343 = 2. * syFC6; auto x344 = x2 * x326; auto x345 = 2. * x271; auto x346 = x345 * x84;
      auto x347 = MUs * x346; auto x348 = x19 * x271; auto x349 = x2 * x348; auto x350 = MUs * x349;
      auto x351 = 4. * x92; auto x352 = 8. * x81; auto x353 = x11 * x259; auto x354 = 4. * x353;
      auto x355 = -x169 * x354 + x325 - x354 * x38; auto x356 = x304 * x55; auto x357 = x304 * x57;
      auto x358 = 3. * x103; auto x359 = x304 * x36; auto x360 = x304 * x43; auto x361 = x266 * x38;
      auto x362 = 4. * MXs; auto x363 = x229 * x24; auto x364 = x11 * pow(MX, 7);
      auto x365 = 2. * x145; auto x366 = 2. * x149; auto x367 = pow(MX, 5); auto x368 = x201 * x367;
      auto x369 = x220 * x367; auto x370 = MXs * x242; auto x371 = x11 * x367;
      auto x372 = x145 * x371; auto x373 = x149 * x371; auto x374 = 2. * syFC8;
      auto x375 = 6. * x332; auto x376 = x229 * x365; auto x377 = x19 * x332; auto x378 = x11 * x20;
      auto x379 = x145 * x229; auto x380 = x229 * x366; auto x381 = x2 * x98; auto x382 = x9 * x98;
      auto x383 = x268 + x361; auto x384 = UU * u; auto x385 = MX * x384; auto x386 = x145 * x385;
      auto x387 = u * x11; auto x388 = x155 * x387; auto x389 = MUs * x344; auto x390 = MXs * x225;
      auto x391 = x390 * x84; auto x392 = MUs * x327; auto x393 = x390 * x88; auto x394 = x3 * x85;
      auto x395 = x23 * x394; auto x396 = x302 * x89; auto x397 = x23 * x313;
      auto x398 = x23 * x314; auto x399 = -x280 * x302; auto x400 = -x261 * x302;
      auto x401 = MXs * x255; auto x402 = MXs * x258; auto x403 = 3. * x39; auto x404 = x403 * x84;
      auto x405 = x345 * x88; auto x406 = x348 * x9;
      auto x407 = MUs * x405 - MUs * x406 - x212 * x251 - x212 * x253 - x23 * x404 - x403 * x89;
      auto x408 = x273 * x51 + x275 * x51; auto x409 = x117 * x239; auto x410 = x122 * x146;
      auto x411 = x117 * x241; auto x412 = MX * x201; auto x413 = x127 * x412;
      auto x414 = x127 * x146; auto x415 = s * x54;
      auto x416 = -x10 * x194 + x10 * x232 - x130 - x131 + x168 - x195 + x22 * x232;
      auto x417 = x145 * x412; auto x418 = x149 * x412; auto x419 = x155 * x220;
      auto x420 = x158 * x220; auto x421 = 2. * syFC5;
      auto x422 = -x141 * x261 - x141 * x280 - x261 * x304 - x280 * x304 + x389 + x391 + x392 +
                  x393 + x395 + x396 + x397 + x398 + x399 + x400 + x401 + x402;
      auto x423 = x293 * x84; auto x424 = x293 * x88; auto x425 = x185 * x367;
      auto x426 = x145 * x242; auto x427 = x19 * x367; auto x428 = u * x427;
      auto x429 = x185 * x229; auto x430 = MXs * x429; auto x431 = 2. * x372; auto x432 = 2. * x373;
      auto x433 = 5. * MUs; auto x434 = x145 * x248; auto x435 = u * x434; auto x436 = MXs * x434;
      auto x437 = 4. * u; auto x438 = x149 * x248; auto x439 = u * x438; auto x440 = MXs * x438;
      auto x441 = UU * x259; auto x442 = x2 * x441; auto x443 = x2 * x92; auto x444 = x9 * x92;
      auto x445 = x2 * x322; auto x446 = x200 * x202; auto x447 = x200 * x206;
      auto x448 = x15 * x259; auto x449 = 7. * x81; auto x450 = MUs * x248; auto x451 = x220 * x229;
      auto x452 = -3. * x69 - 3. * x71; auto x453 = s * syFC4; auto x454 = UU * x271;
      auto x455 = x2 * x454; auto x456 = x454 * x9; auto x457 = UU * x367; auto x458 = x457 * x8;
      auto x459 = x2 * x94; auto x460 = MUs * x459; auto x461 = 10. * MXs; auto x462 = x202 * x371;
      auto x463 = x206 * x371; auto x464 = MXs * UU; auto x465 = x229 * x464;
      auto x466 = MXs * x248; auto x467 = MXs * x459; auto x468 = MXs * x95;
      auto x469 = 6. * x467 + 6. * x468; auto x470 = 4. * x251; auto x471 = 4. * x253;
      auto x472 = -x470 - x471; auto x473 = x148 * x458; auto x474 = x145 * x465;
      auto x475 = x149 * x465; auto x476 = 6. * MXs; auto x477 = 2. * x255;
      auto x478 = -x261 * x476 - x280 * x476 + x477; auto x479 = x149 * x242;
      auto x480 = x244 * x37; auto x481 = -x145 * x480 - x149 * x480 + 4. * x426 + 4. * x479;
      auto x482 = x273 + x275 - x455 - x456; auto x483 = s * syFC6; auto x484 = 3. * x459;
      auto x485 = x145 * x429; auto x486 = 4. * x485; auto x487 = x149 * x429;
      auto x488 = 4. * x487; auto x489 = 16. * MUs; auto x490 = 4. * x417; auto x491 = 4. * x418;
      auto x492 = 2. * x155; auto x493 = x27 * x492; auto x494 = x464 * x492;
      auto x495 = x158 * x27; auto x496 = 2. * x495; auto x497 = x158 * x464; auto x498 = 2. * x497;
      auto x499 = 4. * x190; auto x500 = x379 * x499; auto x501 = x149 * x229 * x499;
      auto x502 = x158 * x387; auto x503 = x191 * x37; auto x504 = x60 + x61; auto x505 = 2. * x96;
      auto x506 = MG * x4; auto x507 = x235 * x506; auto x508 = Mqi * x507; auto x509 = 5. * x371;
      auto x510 = x122 * x248; auto x511 = 7. * MXs; auto x512 = x201 * x231;
      auto x513 = MXs * x200; auto x514 = x201 * x236; auto x515 = MXs * x27;
      auto x516 = x239 * x515; auto x517 = x241 * x515; auto x518 = x239 * x47;
      auto x519 = x241 * x47; auto x520 = x220 * x231; auto x521 = x220 * x236;
      auto x522 = 2. * x460 + 2. * x467 + 2. * x468; auto x523 = s * syFC7; auto x524 = x20 * x88;
      auto x525 = MXs * x118; auto x526 = x203 * x54; auto x527 = 2. * x54; auto x528 = 9. * x35;
      auto x529 = x200 * x30; auto x530 = x158 * x35; auto x531 = 4. * x30; auto x532 = x218 * x30;
      auto x533 = x38 * x77; auto x534 = x169 * x77; auto x535 = x169 * x4; auto x536 = MUs * x170;
      auto x537 = 4. * x536; auto x538 = MUs * x172; auto x539 = 4. * x538; auto x540 = x238 * x506;
      auto x541 = x27 * x38; auto x542 = 6. * x541; auto x543 = x169 * x27; auto x544 = 6. * x543;
      auto x545 = x244 * x51; auto x546 = x149 * x545; auto x547 = x236 * x94;
      auto x548 = -Mqi * x547 + x239 * x40 - x239 * x94 + x241 * x40; auto x549 = 6. * x5;
      auto x550 = 6. * x193; auto x551 = 2. * x381; auto x552 = 2. * x258; auto x553 = 2. * x382;
      auto x554 = x40 * x46; auto x555 = x40 * x49; auto x556 = 6. * x81; auto x557 = s * syFC8;
      auto x558 = x11 * x23; auto x559 = x169 * x558; auto x560 = 3. * x230;
      auto x561 = x145 * x457; auto x562 = x289 * x304 + x473 + x561; auto x563 = 3. * x95;
      auto x564 = x12 * x46; auto x565 = x12 * x49; auto x566 = MXs * x47; auto x567 = 5. * x155;
      auto x568 = x11 * x433; auto x569 = x19 * x37; auto x570 = x158 * x569;
      auto x571 = syFC6 * x107; auto x572 = UU * x229; auto x573 = x145 * x572;
      auto x574 = 3. * x573; auto x575 = x149 * x572; auto x576 = 3. * x575; auto x577 = 2. * x536;
      auto x578 = 2. * x538; auto x579 = 3. * x434; auto x580 = 3. * x438; auto x581 = x149 * x186;
      auto x582 = 5. * MXs; auto x583 = x286 * x582; auto x584 = x287 * x582; auto x585 = MXs * x11;
      auto x586 = x11 * x158; auto x587 = x38 * x464; auto x588 = x169 * x464;
      auto x589 = x282 * x288; auto x590 = x282 * x289; auto x591 = x170 * x282 + x172 * x282;
      auto x592 = -5. * x587 - 5. * x588 - x589 - x590 + x591; auto x593 = 2. * x459;
      auto x594 = 4. * x543; auto x595 = 8. * x536; auto x596 = 8. * x538; auto x597 = x239 * x464;
      auto x598 = x241 * x464; auto x599 = 2. * x417; auto x600 = 2. * x418; auto x601 = -x600;
      auto x602 = x239 * x585; auto x603 = x241 * x585; auto x604 = x19 * x51;
      auto x605 = x155 * x604; auto x606 = x158 * x604; auto x607 = x11 * x51;
      auto x608 = x155 * x607; auto x609 = x11 * x282; auto x610 = x155 * x609;
      auto x611 = x158 * x609; auto x612 = -x158 * x607 - x608 - x610 - x611;
      auto x613 = 2. * x434 + 2. * x438 - 2. * x573 - 2. * x575; auto x614 = syFC7 * x107;
      auto x615 = 4. * x57; auto x616 = x236 * x384; auto x617 = x185 * x203;
      auto x618 = x200 * x617; auto x619 = x185 * x207; auto x620 = x200 * x619;
      auto x621 = x237 * x384; auto x622 = x239 * x387; auto x623 = x241 * x387;
      auto x624 = x215 * x387; auto x625 = x387 * x492; auto x626 = x190 * x203;
      auto x627 = x200 * x626; auto x628 = 2. * x502; auto x629 = x190 * x207;
      auto x630 = x200 * x629; auto x631 = 2. * x187; auto x632 = 2. * x581;
      auto x633 = -x190 * x492 - 2. * x192 - 4. * x541 + x631 + x632; auto x634 = 2. * x118;
      auto x635 = 2. * x133; auto x636 = x288 * x51; auto x637 = 2. * x386; auto x638 = x149 * x385;
      auto x639 = 2. * x638; auto x640 = 9. * MXs; auto x641 = syFC8 * x107; auto x642 = MXs * x109;
      auto x643 = x2 * x642; auto x644 = x642 * x9; auto x645 = 6. * x2; auto x646 = x38 * x47;
      auto x647 = x169 * x47; auto x648 = x289 * x51; auto x649 = x636 + x648;
      auto x650 = x459 + x599 + x600 + x95; auto x651 = x541 + x543;
      auto x652 = x537 + x539 - x554 - x555 + x564 - x583 - x584 - x587 - x588 + x589 + x590 -
                  x628 - x632 + x639 - 3. * x646 - 3. * x647 + x649 + x650 + x651;
      auto x653 = MXs * x133; auto x654 = 11. * x16; auto x655 = 4. * x54; auto x656 = 8. * x34;
      auto x657 = u * x5; auto x658 = 2. * x265; auto x659 = 12. * MXs; auto x660 = x100 * x38;
      auto x661 = 4. * x3; auto x662 = x100 * x169; auto x663 = x38 * x558; auto x664 = x220 * x38;
      auto x665 = x169 * x220; auto x666 = 20. * MXs; auto x667 = x185 * x231;
      auto x668 = x185 * x236; auto x669 = x239 * x27; auto x670 = x190 * x513;
      auto x671 = MUs * x134 + u * x599 + x118 * x37; auto x672 = 12. * MUs; auto x673 = 3. * x155;
      auto x674 = x27 * x673; auto x675 = x464 * x673; auto x676 = 3. * x495; auto x677 = 3. * x497;
      auto x678 = u * x155; auto x679 = x158 * x47; auto x680 = x207 * x387;
      auto x681 = x11 * x155 * x212 + x212 * x586 - x675 - x677; auto x682 = -4. * x5;
      auto x683 = x2 * x515; auto x684 = s * x117; auto x685 = u * x231; auto x686 = s * x27;
      auto x687 = x231 * x384; auto x688 = u * x236; auto x689 = u * x156; auto x690 = x203 * x387;
      auto x691 = syFC2 * x200; auto x692 = 2. * x202; auto x693 = 2. * x206; auto x694 = SS * x231;
      auto x695 = SS * x236; auto x696 = x103 * x82; auto x697 = x40 * x51;
      auto x698 = syFC7 * x200; auto x699 = x15 * x81; auto x700 = s * x421; auto x701 = MUs * x634;
      auto x702 = MUs * x635; auto x703 = x19 * x3; auto x704 = x36 * x51; auto x705 = x43 * x51;
      auto x706 = -x133 - x187 + x191 + x192 - x33 + x386 - x388 + x43 - x502 + x57 - x581 + x638;
      auto x707 = x107 * x309; auto x708 = x215 * x27; auto x709 = x215 * x464;
      auto x710 = x218 * x27; auto x711 = x218 * x464; auto x712 = x122 * x412;
      auto x713 = x215 * x220; auto x714 = x218 * x220; auto x715 = x215 * x47;
      auto x716 = x215 * x585; auto x717 = x218 * x47; auto x718 = x218 * x585;
      auto x719 = x417 + x418; auto x720 = x311 + x312; auto x721 = x107 * x421;
      auto x722 = 4. * x443; auto x723 = 4. * x444; auto x724 = MXs * x313; auto x725 = MXs * x314;
      auto x726 = MXs * x417 + MXs * x418 + x155 * x24 + x158 * x24 - x479;
      auto x727 = x274 + x276 + x552; auto x728 = u * x206; auto x729 = u * x24;
      auto x730 = MXs * x202; auto x731 = u * x202; auto x732 = 2. * x272; auto x733 = 2. * u;
      auto x734 = u * x204; auto x735 = u * x208; auto x736 = x248 * x37; auto x737 = x282 * x40;
      auto x738 = u * x451; auto x739 = x231 * x24; auto x740 = x236 * x24; auto x741 = 4. * x12;
      auto x742 = 2. * x40; auto x743 = x141 * x387; auto x744 = x152 * x190; auto x745 = x2 * x303;
      auto x746 = x303 * x9; auto x747 = x201 * x239; auto x748 = Mqi * x514;
      auto x749 = x220 * x239; auto x750 = x220 * x241; auto x751 = s * x309;
      auto x752 = -MXs * x577 - x152 * x286 - x152 * x287 - x282 * x538 + x282 * x541 +
                  x282 * x543 + x290 + x533 + x534 - x559 - x663;
      auto x753 = MXs * x88; auto x754 = 8. * x3; auto x755 = SS * x229; auto x756 = MXs * SS;
      auto x757 = x155 * x756 + x158 * x756; auto x758 = x19 * x302;
      auto x759 = x2 * x758 + x758 * x9; auto x760 = x301 + x759; auto x761 = x293 * x50;
      auto x762 = 4. * x39; auto x763 = 18. * x3; auto x764 = x294 * x9; auto x765 = 21. * x151;
      auto x766 = x220 * x461; auto x767 = x122 * x755; auto x768 = x127 * x755;
      auto x769 = x241 * x27; auto x770 = x236 * x756; auto x771 = MXs * x19; auto x772 = u * x751;
      auto x773 = x203 * x47; auto x774 = x11 * x222; auto x775 = x19 * x282; auto x776 = s * t;
      ,
      TR * TR * gs * gs * (Nc - 1) * (Nc + 1) *
          (-UU * syFC8 *
               (AXG * x107 * x29 * x645 + AXG * x111 + AXG * x112 - AXG * x182 - AXG * x183 -
                AXG * x643 - AXG * x644 - x111 - x112 + x182 + x183 + x643 + x644) -
           s * syFC1 * (-x494 - x498 + x610 + x611 - x625 - x631 + x637 + x652) +
           s * x296 * (-x5 + x565 + x613 + x652) +
           s * x691 *
               (-MXs * x204 - MXs * x208 + MXs * x221 + MXs * x223 - MXs * x773 - MXs * x774 -
                u * x221 + u * x773 + u * x774 - x12 * x203 - x12 * x207 - x190 * x222 + x231 * x4 +
                x231 * x515 - x231 * x94 + x236 * x515 + x462 + x463 + x507 - x547 + x734 + x735 -
                x739 - x740) -
           syFC1 * x181 + syFC1 * (-x10 * x13 + x6 + x66) +
           2. * syFC1 * (x256 + x257 + x262 + x291 - 2. * x292 - x381 - x382 + x383) -
           4. * syFC3 * (x252 + x254 - x255 + x263 + x268 + x291) -
           syFC4 * x107 *
               (x200 * x680 - x286 * x51 - x287 * x51 - x493 + 2. * x541 + 2. * x543 + x574 + x576 -
                x577 - x578 - x579 - x580 + x608 + x618 + x620 - x621 - x627 - x630 + x649 + x681) -
           syFC6 * (UU * x182 + UU * x183 - x107 * x134 - x108 * x187 + x108 * x191 + x108 * x192 -
                    x108 * x31 - x108 * x33 + x11 * x111 + x11 * x112 - x11 * x182 - x11 * x183 -
                    x114 - x116 - x119 + x171 * x55 + x171 * x57 + x175 * x188 + x176 * x188 -
                    x179 * x189 - x180 * x189 + x184 * x36 + x184 * x43) +
           syFC6 * (6. * s * x96 + 2. * x10 * x104 + x10 * x106 - x10 * x40 * x97 + x10 * x83 +
                    x10 * x93 - x10 * x99 - 4. * x101 * x2 - 4. * x102 + 5. * x105 + x106 * x22 -
                    x14 * x76 + 5. * x15 * x86 + x22 * x83 + x22 * x93 - x22 * x99 - x6 * x76 -
                    x69 * x70 - x70 * x71 - 4. * x73 + 2. * x75 - 5. * x79 - 5. * x80 + x87 + x91) +
           syFC7 * (-x108 * x55 + x109 * x113 + x109 * x115 - x111 * x82 - x112 * x82 + 2. * x114 +
                    2. * x116 + x119 + x122 * x124 - x122 * x128 + x124 * x127 - x127 * x128) -
           syFC7 *
               (x107 * x196 + x107 * x200 * x221 - x107 * x205 + x122 * x213 + x127 * x213 + x136 +
                x138 - x142 * x152 * x2 - x143 * x152 - 2. * x166 - 2. * x167 - x171 * x193 -
                x171 * x5 + 2. * x195 + 4. * x198 + 4. * x199 - x200 * x209 + x200 * x224 +
                x210 * x211 + x211 * x214 - x212 * x216 - x212 * x219 - x216 * x217 - x217 * x219) -
           syFC7 * (-8. * MXs * x239 * x35 - s * x214 * x230 + s * x249 * x250 + x122 * x243 +
                    x127 * x243 + x14 * x37 + x156 * x69 + x156 * x71 + x215 * x240 + x218 * x240 -
                    2. * x227 - 2. * x228 - x231 * x233 - x232 * x234 - x232 * x237 - x233 * x236 +
                    x239 * x240 + x240 * x241 - x245 * x246 - x246 * x247 + x37 * x6 + 2. * x79 +
                    2. * x80 - x87 - x91) +
           syFC8 *
               (x10 * x139 + x107 * x153 + x107 * x154 - x108 * x151 * x18 - 9. * x130 - 9. * x131 -
                4. * x132 - x134 * x26 - x136 - x138 + x139 * x22 + 4. * x140 + x141 * x143 -
                x145 * x147 - x147 * x149 - 8. * x150 + x155 * x157 - x155 * x162 + x155 * x164 +
                x157 * x158 - x158 * x162 + x158 * x164 + 20. * x159 + 20. * x160 + x168) +
           u * x309 *
               (-MXs * x311 - MXs * x312 + MXs * x404 - x251 - x253 - x261 * x37 - x280 * x37 +
                x282 * x746 + x315 + x316 - x344 - x346 + x349 + x403 * x753 - x405 + x406 + x478 +
                x722 + x723 - x724 - x725 + x727) -
           u * x374 *
               (-MXs * x423 - MXs * x424 - x18 * x352 + 2. * x256 + x327 + x344 + x472 + x477 +
                x546 + x551 + x553 - 8. * x699 - x722 - x723 + x724 + x725 + x726 + x727) -
           u * x523 *
               (x196 - x200 * x208 + x200 * x520 + x200 * x521 - x205 + x215 * x569 + x218 * x569 +
                4. * x245 + 4. * x247 - 4. * x249 - 4. * x510 - 3. * x519 + x593 + x682 +
                3. * x708 + 3. * x709 + 3. * x710 + 3. * x711 - 7. * x715 - 7. * x716 - 7. * x717 -
                7. * x718 - 4. * x767 - 4. * x768 + 3. * x769) +
           u * x557 *
               (-x141 * x84 - x15 * x763 - x15 * x765 - x18 * x763 - x18 * x765 + x2 * x766 +
                5. * x261 + 5. * x280 + x304 * x50 - x37 * x753 - x459 - x46 * x558 + x46 * x77 +
                9. * x515 * x9 + x549 + x550 + x570 + x601 + 9. * x683 + x759 + x761 - x762 * x84 -
                x762 * x88 + x764 + x766 * x9 - x95) +
           x107 * x691 *
               (-x203 * x27 + x204 - x221 + x512 + x514 + x616 - x617 - x619 + x626 + x629 - x667 -
                x668 - x680 + x687 - x690 + x773 + x774) +
           x2 * x29 * x721 * (-SS + UU - x11 + x19) +
           x29 * x374 *
               (-x15 * x754 - x18 * x754 + x244 * x366 + x300 + x313 + x314 - x366 * x755 + x419 +
                x420 - x423 - x424 + x612 + x719 + x757 + x760) +
           x29 * x698 *
               (x203 * x607 + x203 * x609 - x203 * x756 + x207 * x609 - x207 * x756 + x231 * x607 +
                x231 * x756 - x231 * x775 + 4. * x236 * x585 + x236 * x607 - x236 * x775 -
                x244 * x692 - x244 * x693 + x248 * x692 + x248 * x693 - x512 - x514 - x520 - x521 +
                x692 * x755 + x693 * x755 + x770) -
           x296 * (-4. * x292 + x295 + x66) +
           x296 * (-x156 * x187 + x156 * x386 - x156 * x388 + x181) -
           x309 * (MXs * x315 + MXs * x316 - x286 * x321 - x287 * x321 - x293 * x311 - x293 * x312 +
                   x299 * x313 + x299 * x314 + x3 * x311 + x3 * x312 - x317 - x318 - x319 - x320 +
                   x325 + x331) +
           x309 * (s * x409 + s * x410 + s * x411 + s * x414 + x107 * x413 + x161 * x215 +
                   x161 * x239 + x161 * x241 - x163 * x215 - x163 * x218 - x163 * x239 -
                   x163 * x241 + x198 + x199 - x215 * x415 - x239 * x415 - x241 * x415 + x416) +
           x309 * (-x152 * x261 - x152 * x280 + x347 - x350 + x389 + x391 + x392 + x393 + x395 +
                   x396 + x397 + x398 + x399 + x400 + x401 + x402 + x407 + x408) -
           x309 * (-x105 + x227 + x228 + x26 * x304 * x53 - x282 * x303 * x38 - x29 * x300 -
                   x29 * x301 + x293 * x31 + x293 * x33 - x293 * x55 - x293 * x57 + 2. * x297 +
                   2. * x298 + x299 * x31 + x299 * x33 + x302 * x36 + x302 * x43 + x305 + x308 +
                   x80 - x90) +
           x343 * (MUs * x330 + x407 + x422) -
           x343 * (-MUs * x329 + u * x344 - x169 * x351 - x286 * x352 - x287 * x352 + x317 + x318 +
                   x319 + x320 + x331 - x347 + x350 - x351 * x38 + x355) -
           x343 * (-x261 * x334 - x280 * x334 + x297 + x298 + x3 * x338 + 6. * x3 * x36 + x305 +
                   x31 * x39 + x33 * x39 - x332 * x333 - x332 * x335 + x332 * x336 + x332 * x337 +
                   x339 * x44 + x342) -
           x343 * (MUs * x56 + MUs * x58 + MXs * x267 + MXs * x361 + MXs * x42 - x169 * x348 +
                   x170 * x345 - x170 * x358 + x172 * x345 - x172 * x358 - x212 * x277 -
                   x212 * x292 + x285 * x362 + x285 * x37 + x308 - x32 * x37 + x339 * x45 -
                   x34 * x37 - x348 * x38 - x356 - x357 + x359 + x360 + x37 * x41) +
           x374 * (-MXs * x252 - MXs * x254 + MXs * x276 - x23 * x423 - x23 * x424 + x408 + x422) -
           x374 * (x102 - x15 * x375 + x155 * x232 - x155 * x378 + x158 * x232 - x158 * x378 -
                   x177 * x379 - x18 * x375 + x342 + x35 * x376 + x35 * x380 + x376 * x54 +
                   x377 * x46 + x377 * x49 + 2. * x73 - x75) -
           x374 * (-MXs * x274 - x145 * x363 - x145 * x369 - x145 * x370 - x149 * x363 -
                   x149 * x369 - x149 * x370 + x265 * x51 + x282 * x372 + x282 * x373 + 2. * x318 +
                   2. * x320 + x355 - x364 * x365 - x364 * x366 + x365 * x368 + x366 * x368 +
                   x372 * x51 + x373 * x51) +
           x374 * (u * x426 - u * x431 - u * x432 - x145 * x428 - x145 * x430 - x149 * x428 -
                   x149 * x430 - x190 * x379 * x51 + x306 + x307 + x32 * x433 + x34 * x433 + x356 +
                   x357 - x359 - x360 + x365 * x425 + x366 * x425 + x37 * x435 + x37 * x439 +
                   x436 * x437 + x437 * x440 - x44 * x76 - x45 * x76) +
           x421 * (MUs * x175 + MUs * x176 - MUs * x179 - MUs * x180 + MXs * x175 + MXs * x176 -
                   MXs * x179 - MXs * x180 - x107 * x417 - x107 * x418 + x107 * x419 + x107 * x420 -
                   x132 - x133 * x26 + x140 - x150 + x159 + x160 + x416) +
           x453 * (x15 * x449 + x18 * x449 - x193 * x433 + x230 * x446 + x230 * x447 - x242 * x446 -
                   x242 * x447 + 4. * x256 + 4. * x257 + 2. * x260 - 4. * x381 - 4. * x382 -
                   x433 * x5 - 2. * x442 + 5. * x443 + 5. * x444 - 4. * x445 - x446 * x450 +
                   x446 * x451 - x447 * x450 + x447 * x451 + 2. * x448 + x452) +
           x453 * (-x125 * x200 * x458 + x200 * x462 + x200 * x463 - x250 * x261 - x250 * x280 +
                   4. * x255 + 4. * x258 - x261 * x461 + 5. * x273 + 5. * x275 - x280 * x461 +
                   2. * x316 - x446 * x457 + x446 * x465 - x446 * x466 + x447 * x465 - x447 * x466 -
                   3. * x455 - 3. * x456 + 5. * x460 + x469 + x472 + 5. * x96) +
           x483 * (-x230 * x366 - x27 * x376 - x280 * x97 - 6. * x373 + x434 * x476 + x434 * x97 +
                   x438 * x476 + x438 * x97 + 6. * x460 + x469 + 2. * x473 - 2. * x474 - 2. * x475 +
                   x478 + x481 + x482) +
           x483 * (-u * x484 - u * x490 - u * x491 + u * x493 + u * x494 + u * x496 + u * x498 -
                   x133 * x339 + x192 * x37 + 11. * x285 - x295 - x31 * x76 - x33 * x76 -
                   x388 * x476 - x388 * x97 + 11. * x41 + x43 * x489 - x476 * x502 - x486 - x488 +
                   x500 + x501 - x502 * x97 + x503 + x504) -
           x483 * (MXs * x595 + MXs * x596 + u * x563 + x169 * x294 + 2. * x264 + 15. * x278 +
                   15. * x279 - x288 * x661 - x289 * x661 + 6. * x372 - 5. * x533 - 5. * x534 -
                   7. * x535 - x541 * x659 - x543 * x659 + 5. * x559 - 2. * x561 + x646 * x666 +
                   x647 * x666 - 7. * x657 + x658 - x659 * x664 - x659 * x665 - 4. * x660 -
                   4. * x662 + 5. * x663) -
           x483 * (x118 * x339 + x146 * x365 + x146 * x366 + x15 * x654 - 6. * x155 * x35 +
                   x155 * x531 - x155 * x655 + x158 * x531 - x158 * x655 - x165 * x46 - x165 * x49 +
                   x18 * x654 + x20 * x84 - x232 * x645 - 6. * x232 * x9 + 8. * x32 - x36 * x489 -
                   15. * x44 - 15. * x45 + x524 + 7. * x525 - 6. * x530 - 6. * x56 - 6. * x58 +
                   7. * x653 + x656) +
           x523 * (-MXs * x537 - MXs * x539 + MXs * x542 + MXs * x544 - x13 * x239 - x145 * x545 -
                   6. * x269 - 6. * x270 + x282 * x434 + x282 * x438 + 2. * x426 - x431 - x432 +
                   x434 * x51 + x438 * x51 + 2. * x479 + 2. * x533 + 2. * x534 + 4. * x535 +
                   3. * x540 - x546 + x548) +
           x523 * (-x122 * x457 + x122 * x509 - x126 * x458 + x127 * x509 - x13 * x241 +
                   x210 * x230 + x210 * x465 - x212 * x518 - x212 * x519 + x214 * x465 -
                   x249 * x511 - x250 * x510 + x505 + 3. * x508 - x510 * x511 - x512 * x513 -
                   x513 * x514 + x513 * x520 + x513 * x521 + 3. * x516 + 3. * x517 + x522) +
           x523 * (-MXs * x134 - x155 * x177 + x155 * x527 - x158 * x177 + x158 * x527 -
                   x178 * x239 - x178 * x241 - x200 * x526 - x215 * x528 + x215 * x531 -
                   x218 * x528 + x231 * x529 + x236 * x529 - x237 * x54 + x35 * x492 + x409 + x410 +
                   x411 + x414 - x524 - 4. * x525 + 2. * x530 + 4. * x532 + x59) -
           x523 * (u * x600 - x188 * x518 + x188 * x597 + x188 * x598 + x188 * x669 - x190 * x376 -
                   x190 * x380 - x191 * x51 - x192 * x51 + x231 * x670 + x236 * x670 + x282 * x388 +
                   x282 * x502 + x388 * x51 + 2. * x485 + 2. * x487 + x502 * x51 - x511 * x622 -
                   x511 * x623 - x513 * x667 - x513 * x668 - x55 * x97 - x57 * x97 + x671) +
           x557 * (MXs * x554 - x145 * x560 - x149 * x560 + x250 * x434 + x250 * x438 + 2. * x269 +
                   2. * x270 - x329 - x330 - 5. * x372 - 5. * x373 + 7. * x436 + 7. * x440 + x467 -
                   3. * x474 - 3. * x475 + x481 + 3. * x534 - 3. * x559 + x562) -
           x557 * (-MXs * x615 - u * x674 - u * x675 - u * x676 - u * x677 + x189 * x679 +
                   x31 * x97 + x33 * x97 - x36 * x672 + x388 * x511 - x43 * x672 + 4. * x435 +
                   4. * x439 + 7. * x47 * x678 + x486 + x488 - x500 - x501 + x502 * x511 - x503 +
                   x504 + x656 + x671) +
           x557 * (-MUs * x549 - MUs * x550 + MUs * x554 + MUs * x555 + MXs * x555 + x15 * x556 +
                   x18 * x556 + x353 * x46 + x353 * x49 - x441 * x46 - x441 * x49 - 2. * x445 +
                   x452 + x46 * x92 + x460 + x468 + x470 + x471 - x477 + x49 * x92 - x551 - x552 -
                   x553 + x96) +
           x571 * (x286 * x97 + x287 * x97 - x288 * x37 - x289 * x37 + x464 * x567 + 5. * x497 -
                   x542 - x544 - x567 * x585 - x574 - x576 + x577 + x578 + x579 + x580 + 4. * x581 -
                   x582 * x586 + x583 + x584 + x592) +
           x571 *
               (-MXs * x113 - MXs * x115 + x153 + x154 - x155 * x568 + x155 * x569 - x158 * x568 +
                x27 * x567 + x4 * x49 + x46 * x515 - x46 * x566 - x484 + x49 * x515 - x49 * x566 -
                x490 - x491 + 5. * x495 + 3. * x5 + x554 + x555 - x563 - x564 - x565 + x570) +
           x614 * (x493 + x494 + x496 + x498 - x593 - x594 + x595 + x596 - x597 - x598 - x599 +
                   x601 + x602 + x603 + x605 + x606 + x612 + x613) +
           x614 * (Mqi * x616 + x134 - 6. * x218 * x387 + x234 * x384 + x239 * x384 - 4. * x587 -
                   4. * x588 + x591 - x615 - x618 - x620 + x621 - x622 - x623 - 6. * x624 + x625 +
                   x627 + x628 + x630 + x633) +
           x641 * (x286 * x37 + x286 * x640 + x287 * x37 + x287 * x640 + 2. * x31 + 2. * x33 -
                   x338 + 6. * x388 + 6. * x502 + x592 + x633 + x634 + x635 - x636 - x637 - x639) -
           x641 *
               (x141 * x15 + x15 * x661 + x18 * x661 - 4. * x193 - x261 - x280 + x47 * x673 + x594 -
                x605 - x606 + x648 + x650 - x674 - x676 + 3. * x679 + x681 + x682 - 4. * x683) -
           x691 *
               (x142 * x222 - x161 * x203 - x161 * x207 + x163 * x203 + x163 * x207 + x203 * x415 +
                x207 * x415 - x209 + x224 + x231 * x415 - x231 * x684 + x236 * x415 - x236 * x684 +
                x26 * x616 - x26 * x617 - x26 * x619 + x26 * x626 + x26 * x629 - x26 * x680 +
                x26 * x687 - x26 * x690 - x512 * x689 - x514 * x689 + x685 * x686 + x686 * x688) +
           x698 *
               (MXs * x734 + MXs * x735 + x203 * x729 + x206 * x430 + x207 * x729 + x231 * x697 -
                x231 * x732 + x231 * x737 - x236 * x732 + x236 * x737 - x242 * x728 - x242 * x731 -
                x248 * x362 * x728 - x248 * x437 * x730 - x425 * x692 - x425 * x693 + x427 * x728 +
                x429 * x730 + x463 * x733 + x692 * x738 + x693 * x738 - x728 * x736 - x731 * x736) -
           x698 *
               (x104 * x231 + x104 * x236 + x202 * x363 + x202 * x369 + x202 * x370 + x206 * x363 +
                x206 * x369 + x206 * x370 - x231 * x696 - x236 * x696 - x236 * x697 + x259 * x694 +
                x259 * x695 - x282 * x462 - x282 * x463 + x364 * x692 + x364 * x693 - x368 * x692 -
                x368 * x693 - x462 * x51 - x463 * x51 + x694 * x81 + x695 * x81) +
           x698 *
               (-MUs * x203 * x30 - MUs * x526 + u * x739 + u * x740 + x207 * x35 * x51 -
                x222 * x30 - x222 * x54 - x231 * x743 + x231 * x744 - x236 * x743 + x236 * x744 +
                x266 * x685 + x266 * x688 - x293 * x667 - x293 * x668 + x3 * x667 + x3 * x668 +
                x427 * x731 + x462 * x733 - x685 * x741 + x685 * x742 - x688 * x741 + x688 * x742) +
           x700 * (x255 - x269 - x270 - x281 - x283 - x284 + x482 + x522 + x660 + x662 + x752) -
           x700 * (-x18 * x259 + x193 * x51 + x251 + x253 + x263 + x322 * x9 - x341 * x81 + x381 +
                   x382 + x441 * x9 + x442 - x443 - x444 + x445 - x448 + x5 * x51 - x505 + x69 -
                   2. * x699 + x71) -
           x700 * (u * x459 + u * x95 - x169 * x703 - x277 + x278 + x279 - x285 - x292 + x34 -
                   x38 * x703 + x383 - x41 - x535 - x58 + x62 + x63 + x64 - x657 + x701 + x702 -
                   x704 - x705) -
           x707 * (-x118 - x122 * x385 - x31 + x36 + x55 + x624 + x706) -
           x707 *
               (-x122 * x572 - x127 * x572 + x249 - x419 - x420 - x495 + x510 + x679 + x708 + x709 +
                x710 + x711 - x712 + x713 + x714 - x715 - x716 - x717 - x718 + x719 + x720) -
           x721 * (-x434 - x438 - x536 - x538 + x573 + x575 - x646 - x647 + x651 + x664 + x665 +
                   x706) +
           x751 * (-x155 * x310 - x158 * x310 + x288 * x304 - x372 - x373 + x436 + x548 + x562 -
                   x658 + x752) +
           x751 * (MXs * x745 + MXs * x746 - MXs * x747 - MXs * x748 + MXs * x749 + MXs * x750 -
                   x12 * x239 - x12 * x241 - x239 * x566 - x241 * x566 + x315 - x426 + x440 - x474 -
                   x475 + x508 + x516 + x517 + x540 + x726) +
           x751 * (-MXs * x388 - u * x573 - u * x575 - x218 * x54 - x32 - x34 + x435 + x439 + x44 +
                   x45 + x464 * x678 + x485 + x487 + x504 - x525 + x532 + x56 + x58 + x65 - x653 -
                   x701 - x702 + x704 + x705) +
           x772 * (-x15 * x302 - x245 - x247 - x3 * x341 - x394 + x413 + x423 + x424 - x708 - x709 -
                   x710 - x711 + x712 - x713 - x714 + x715 + x716 + x717 + x718 + x760 - x764 +
                   x767 + x768) -
           x772 * (-Mqi * x770 + x158 * x585 - x239 * x756 + x239 * x771 + x241 * x771 - x300 -
                   x497 - x518 - x519 + x597 + x598 - x602 - x603 + x669 + x720 + x745 + x746 -
                   x747 - x748 + x749 + x750 + x757 + x761 + x769)) *
          Denom(768 *
                (MGs * x107 - MGs * x26 + MGs * x776 + MXs * x776 - s * pow(t, 2) - t * x107) *
                pow(Pi, 2)));

  if (color_flow) {
    return -ret.real();
  }

  return ret.real();
}

ComplexType ME_us_gGX_Qqq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;
  for (int invert_color_flow = 0; invert_color_flow < 2; invert_color_flow++) {
    for (int itq = 0; itq < 6; itq++) {
      for (int itsq = 0; itsq < 2; itsq++) {
        ret += ME_us_gGX_Qqq_single(pIEPS, sc, uc, axial, itsq, itq, invert_color_flow, Q2, P1K1,
                                    params);
      }
    }
  }
  return ret;
}
