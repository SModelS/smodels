#include "pxs_gausq_3.h"

#define C0(a, b, c, d, e, f) C0i(cc0, a, b, c, d, e, f)
#define C1(a, b, c, d, e, f) C0i(cc1, a, b, c, d, e, f)
#define C2(a, b, c, d, e, f) C0i(cc2, a, b, c, d, e, f)
#define C00(a, b, c, d, e, f) C0i(cc00, a, b, c, d, e, f)
#define C11(a, b, c, d, e, f) C0i(cc11, a, b, c, d, e, f)
#define C12(a, b, c, d, e, f) C0i(cc12, a, b, c, d, e, f)
#define C22(a, b, c, d, e, f) C0i(cc22, a, b, c, d, e, f)

#define D0(a, b, c, d, e, f, g, h, i, j) D0i(dd0, a, b, c, d, e, f, g, h, i, j)
#define D1(a, b, c, d, e, f, g, h, i, j) D0i(dd1, a, b, c, d, e, f, g, h, i, j)
#define D2(a, b, c, d, e, f, g, h, i, j) D0i(dd2, a, b, c, d, e, f, g, h, i, j)
#define D3(a, b, c, d, e, f, g, h, i, j) D0i(dd3, a, b, c, d, e, f, g, h, i, j)
#define D00(a, b, c, d, e, f, g, h, i, j) D0i(dd00, a, b, c, d, e, f, g, h, i, j)
#define D11(a, b, c, d, e, f, g, h, i, j) D0i(dd11, a, b, c, d, e, f, g, h, i, j)
#define D12(a, b, c, d, e, f, g, h, i, j) D0i(dd12, a, b, c, d, e, f, g, h, i, j)
#define D22(a, b, c, d, e, f, g, h, i, j) D0i(dd22, a, b, c, d, e, f, g, h, i, j)
#define D13(a, b, c, d, e, f, g, h, i, j) D0i(dd13, a, b, c, d, e, f, g, h, i, j)
#define D23(a, b, c, d, e, f, g, h, i, j) D0i(dd23, a, b, c, d, e, f, g, h, i, j)
#define D33(a, b, c, d, e, f, g, h, i, j) D0i(dd33, a, b, c, d, e, f, g, h, i, j)

ComplexType ME_us_box_QGGq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  int itq = q;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int ftq = 0; ftq < 2; ftq++) {
    int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
    // int iq = (itq + ftq * 3) % 6;
    int iq = (sq - is_up_squark(sq) * 6) % 3 + is_up_squark(sq) * 3;
    ;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
    ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

    ComplexType kL = (params->CHSQq[ch][isq][iq].L);
    ComplexType kR = (params->CHSQq[ch][isq][iq].R);

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    ComplexType jLGp = conj(params->GLSQq[isq][q].R);
    ComplexType jRGp = conj(params->GLSQq[isq][q].L);
    ComplexType iLGp = conj(params->GLSQq[sq][iq].R);
    ComplexType iRGp = conj(params->GLSQq[sq][iq].L);

    auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C0(MUs, 0, u, Mqis, MGs, MGs)
#define syFC2 D0(MUs, 0, 0, MXs, u, s, Mqis, MGs, MGs, MQis)
#define syFC3 D1(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC4 D1(MUs, u, 0, s, 0, MXs, MGs, Mqis, MGs, MQis)
#define syFC5 C2(0, u, MXs, MQis, MGs, Mqis)
#define syFC6 D2(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC7 D2(MUs, u, 0, s, 0, MXs, MGs, Mqis, MGs, MQis)
#define syFC8 D3(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC9 D3(MUs, u, 0, s, 0, MXs, MGs, Mqis, MGs, MQis)
#define syFC10 D00(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC11 D12(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC12 D13(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC13 D22(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC14 D23(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)
#define syFC15 D33(0, s, MUs, u, 0, MXs, MGs, MQis, MGs, Mqis)


    _EPS0_(
        ret, auto x0 = pow(s, 2); auto x1 = u * x0; auto x2 = UU * x1; auto x3 = L * jRGp;
        auto x4 = iLGp * kR; auto x5 = x3 * x4; auto x6 = x2 * x5; auto x7 = MUs * SS;
        auto x8 = x0 * x7; auto x9 = R * jLGp; auto x10 = iRGp * x9; auto x11 = kL * x10;
        auto x12 = 4.0 * x11; auto x13 = pow(MXs, 2); auto x14 = UU * s; auto x15 = x13 * x14;
        auto x16 = SS * x1; auto x17 = x16 * x5; auto x18 = 4.0 * x17; auto x19 = 4.0 * s;
        auto x20 = pow(MUs, 2); auto x21 = SS * x20; auto x22 = x11 * x21; auto x23 = pow(u, 2);
        auto x24 = SS * x23; auto x25 = x10 * x24; auto x26 = 4.0 * x25; auto x27 = kL * s;
        auto x28 = pow(MX, 4); auto x29 = x14 * x28; auto x30 = AXG * x2; auto x31 = x30 * x5;
        auto x32 = x14 * x23; auto x33 = AXG * x32; auto x34 = MG * iRGp; auto x35 = x3 * x34;
        auto x36 = pow(MX, 3); auto x37 = kL * x36; auto x38 = 4.0 * x14; auto x39 = x37 * x38;
        auto x40 = x35 * x39; auto x41 = Mqi * iLGp; auto x42 = x3 * x41; auto x43 = x39 * x42;
        auto x44 = MXs * x14; auto x45 = MUs * x44; auto x46 = SS * u; auto x47 = MXs * s;
        auto x48 = x46 * x47; auto x49 = s * u; auto x50 = x49 * x7; auto x51 = 8.0 * x50;
        auto x52 = 2.0 * x0; auto x53 = UU * x3; auto x54 = x41 * x53; auto x55 = MX * x54;
        auto x56 = kL * x55; auto x57 = 4.0 * x7; auto x58 = x11 * x47; auto x59 = MUs * u;
        auto x60 = x14 * x59; auto x61 = MX * x35; auto x62 = 4.0 * x46; auto x63 = x27 * x62;
        auto x64 = MX * x42; auto x65 = u * x44; auto x66 = AXG * x65; auto x67 = x12 * x66;
        auto x68 = x27 * x57; auto x69 = kL * x35; auto x70 = u * x38; auto x71 = MX * x70;
        auto x72 = kL * x42; auto x73 = AXG * u; auto x74 = x38 * x73; auto x75 = MX * x74;
        auto x76 = UU * x23; auto x77 = kR * x0; auto x78 = iLGp * x3; auto x79 = x77 * x78;
        auto x80 = x24 * x78; auto x81 = AXG * x76; auto x82 = x79 * x81; auto x83 = MG * x9;
        auto x84 = iLGp * x83; auto x85 = UU * x84; auto x86 = x36 * x85; auto x87 = UU * x10;
        auto x88 = Mqi * x87; auto x89 = x36 * x88; auto x90 = x77 * x89; auto x91 = x24 * x5;
        auto x92 = AXG * x91; auto x93 = 2.0 * SS; auto x94 = MXs * x93; auto x95 = x20 * x94;
        auto x96 = s * x11; auto x97 = 2.0 * x11; auto x98 = x20 * x44 * x97; auto x99 = AXG * x86;
        auto x100 = x77 * x99; auto x101 = AXG * x89; auto x102 = x101 * x77;
        auto x103 = MUs * iLGp; auto x104 = UU * x103 * x83; auto x105 = MX * x104;
        auto x106 = MX * MXs; auto x107 = x106 * x85; auto x108 = x4 * x83; auto x109 = MX * x2;
        auto x110 = x108 * x109; auto x111 = x10 * x2; auto x112 = MX * Mqi;
        auto x113 = kR * x111 * x112; auto x114 = AXG * MUs; auto x115 = 2.0 * x16;
        auto x116 = x115 * x5; auto x117 = MX * x30; auto x118 = x108 * x117; auto x119 = Mqi * x10;
        auto x120 = kR * x119; auto x121 = x117 * x120; auto x122 = MX * x85;
        auto x123 = x114 * x122; auto x124 = AXG * MXs; auto x125 = x122 * x124;
        auto x126 = MX * x88; auto x127 = x126 * x77; auto x128 = 3.0 * x76; auto x129 = 6.0 * x28;
        auto x130 = pow(MXs, 3); auto x131 = 2.0 * x7; auto x132 = pow(MUs, 3);
        auto x133 = 4.0 * x13; auto x134 = MGs * x6; auto x135 = AXG * pow(u, 3);
        auto x136 = 4.0 * x135; auto x137 = x136 * x14; auto x138 = x11 * x7; auto x139 = AXG * x28;
        auto x140 = x139 * x93; auto x141 = 2.0 * x8; auto x142 = x141 * x5; auto x143 = MGs * x31;
        auto x144 = 4.0 * x28; auto x145 = UU * x144; auto x146 = x114 * x145;
        auto x147 = x11 * x146; auto x148 = x131 * x139; auto x149 = x11 * x148;
        auto x150 = 2.0 * x110; auto x151 = kR * s; auto x152 = x78 * x94; auto x153 = u * x152;
        auto x154 = x151 * x153; auto x155 = x15 * x5; auto x156 = x15 * x97; auto x157 = x20 * x78;
        auto x158 = x157 * x93; auto x159 = x29 * x5; auto x160 = s * x93; auto x161 = x160 * x20;
        auto x162 = x11 * x161; auto x163 = x33 * x5; auto x164 = x45 * x5; auto x165 = x45 * x97;
        auto x166 = x131 * x78; auto x167 = MXs * x166; auto x168 = x5 * x60;
        auto x169 = x131 * x47; auto x170 = x11 * x169; auto x171 = AXG * x127;
        auto x172 = 2.0 * MX; auto x173 = x172 * x44; auto x174 = x120 * x173; auto x175 = u * x93;
        auto x176 = MX * x119; auto x177 = x151 * x176; auto x178 = x5 * x66;
        auto x179 = x131 * x64; auto x180 = x179 * x27; auto x181 = u * x14;
        auto x182 = x120 * x172; auto x183 = x14 * x73; auto x184 = x78 * x93;
        auto x185 = x184 * x23; auto x186 = -x151 * x185; auto x187 = x57 * x78;
        auto x188 = u * x187; auto x189 = -x116 + x142 + x151 * x188 + x186 + x31 + x6;
        auto x190 = 2.0 * MGs; auto x191 = x5 * x8; auto x192 = -x143; auto x193 = Mqi * x35;
        auto x194 = kR * x193; auto x195 = 2.0 * x32; auto x196 = x194 * x2; auto x197 = MX * x195;
        auto x198 = x23 * x93; auto x199 = x193 * x198; auto x200 = MX * x84;
        auto x201 = pow(MX, 5); auto x202 = x14 * x201; auto x203 = 2.0 * x202;
        auto x204 = x203 * x69; auto x205 = x11 * x190; auto x206 = x205 * x29;
        auto x207 = 2.0 * x29; auto x208 = x41 * x83; auto x209 = kL * x208;
        auto x210 = x207 * x209; auto x211 = x195 * x208; auto x212 = AXG * x29;
        auto x213 = x11 * x212; auto x214 = 2.0 * x209; auto x215 = 2.0 * x33;
        auto x216 = 4.0 * x65; auto x217 = 4.0 * x66; auto x218 = 4.0 * UU; auto x219 = pow(MX, 6);
        auto x220 = AXG * x219; auto x221 = x218 * x220; auto x222 = x144 * x46;
        auto x223 = x175 * x20; auto x224 = x13 * x175; auto x225 = u * x166; auto x226 = AXG * x13;
        auto x227 = UU * u; auto x228 = 8.0 * x227; auto x229 = x226 * x228; auto x230 = 3.0 * x33;
        auto x231 = x139 * x175; auto x232 = AXG * x145; auto x233 = MXs * x232;
        auto x234 = 8.0 * x124; auto x235 = x59 * x87; auto x236 = 4.0 * x124; auto x237 = u * x7;
        auto x238 = x10 * x237; auto x239 = kL * x238; auto x240 = 3.0 * x66; auto x241 = AXG * x15;
        auto x242 = AXG * x60; auto x243 = -x212 * x5 + x241 * x5 - x242 * x5;
        auto x244 = x13 * x131; auto x245 = x144 * x7;
        auto x246 = x11 * x244 - x11 * x245 + x11 * x95 + x149 + x159 + x168 + x32 * x5 - x5 * x65;
        auto x247 = pow(s, 3); auto x248 = u * x247; auto x249 = x248 * x5; auto x250 = UU * x249;
        auto x251 = x247 * x5; auto x252 = x5 * x76; auto x253 = MX * x247; auto x254 = kR * x253;
        auto x255 = x254 * x88; auto x256 = x218 * x23; auto x257 = AXG * x256;
        auto x258 = x21 * x96; auto x259 = 2.0 * MUs; auto x260 = x35 * x37; auto x261 = x14 * x260;
        auto x262 = 2.0 * x260 * x44; auto x263 = x131 * x260; auto x264 = AXG * s;
        auto x265 = 2.0 * x114; auto x266 = MXs * x264; auto x267 = x131 * x266;
        auto x268 = AXG * x45; auto x269 = MUs * x10; auto x270 = x10 * x20; auto x271 = x270 * x44;
        auto x272 = x10 * x131; auto x273 = x272 * x28; auto x274 = x135 * x218;
        auto x275 = x36 * x42; auto x276 = x14 * x275; auto x277 = x275 * x44;
        auto x278 = x135 * x93; auto x279 = MX * x278; auto x280 = MX * x274;
        auto x281 = x10 * x244; auto x282 = s * x10; auto x283 = x35 * x36; auto x284 = 4.0 * x24;
        auto x285 = x283 * x284; auto x286 = x124 * x256;
        auto x287 = -AXG * x285 + s * x281 + x257 * x283 + x279 * x35 - x280 * x35 + x282 * x95 +
                    x285 + x286 * x61;
        auto x288 = kL * syFC15; auto x289 = x10 * x76; auto x290 = 3.0 * x289;
        auto x291 = AXG * x289; auto x292 = x115 * x269; auto x293 = x10 * x137;
        auto x294 = x10 * x23; auto x295 = s * x294; auto x296 = x295 * x57;
        auto x297 = x278 * x282; auto x298 = 3.0 * x32; auto x299 = x269 * x30;
        auto x300 = x131 * x23; auto x301 = x10 * x300; auto x302 = x23 * x94;
        auto x303 = x10 * x264; auto x304 = x302 * x303; auto x305 = x160 * x23;
        auto x306 = MX * x305; auto x307 = x306 * x35; auto x308 = u * x283; auto x309 = MXs * x10;
        auto x310 = x309 * x33; auto x311 = 4.0 * x310; auto x312 = MX * x298;
        auto x313 = 5.0 * x33; auto x314 = MX * x33; auto x315 = x314 * x35; auto x316 = x314 * x42;
        auto x317 = x10 * x305; auto x318 = x10 * x190; auto x319 = x318 * x32;
        auto x320 = x10 * x141; auto x321 = x24 * x264; auto x322 = x208 * x305;
        auto x323 = 2.0 * x10; auto x324 = x323 * x33; auto x325 = x195 * x61;
        auto x326 = x318 * x7; auto x327 = x172 * x35; auto x328 = x327 * x65; auto x329 = MX * x46;
        auto x330 = x329 * x42; auto x331 = x172 * x42; auto x332 = x181 * x331;
        auto x333 = 2.0 * x208; auto x334 = x266 * x46; auto x335 = x190 * x64;
        auto x336 = x49 * x57; auto x337 = x183 * x331; auto x338 = kL * syFC2; auto x339 = x0 * x3;
        auto x340 = x339 * x36; auto x341 = x34 * x340; auto x342 = UU * x341; auto x343 = MXs * x0;
        auto x344 = MGs * x87; auto x345 = x10 * x30; auto x346 = UU * x41; auto x347 = x346 * x83;
        auto x348 = x0 * x124; auto x349 = x327 * x33; auto x350 = x300 * x61;
        auto x351 = x302 * x61; auto x352 = x114 * x256; auto x353 = x352 * x61;
        auto x354 = x234 * x76; auto x355 = AXG * x300; auto x356 = x355 * x61;
        auto x357 = x236 * x24; auto x358 = kL * syFC8; auto x359 = AXG * x218;
        auto x360 = x359 * pow(MX, 8); auto x361 = x219 * x57; auto x362 = x201 * x35;
        auto x363 = x362 * x57; auto x364 = x20 * x93; auto x365 = x283 * x364;
        auto x366 = x221 * x269; auto x367 = 8.0 * MXs; auto x368 = x367 * x87;
        auto x369 = x218 * x362; auto x370 = x114 * x369; auto x371 = x124 * x369;
        auto x372 = AXG * x201; auto x373 = x35 * x372; auto x374 = x131 * x373;
        auto x375 = x131 * x283; auto x376 = MXs * x375; auto x377 = pow(MU, 4);
        auto x378 = x10 * x377; auto x379 = x13 * x93; auto x380 = SS * x378;
        auto x381 = x131 * x378; auto x382 = MXs * x381 + x140 * x378 - x144 * x380 + x378 * x379;
        auto x383 = x294 * x93; auto x384 = AXG * x383; auto x385 = x377 * x384;
        auto x386 = AXG * x25; auto x387 = x10 * x256; auto x388 = x20 * x383;
        auto x389 = x362 * x62; auto x390 = x175 * x373; auto x391 = x308 * x94;
        auto x392 = u * x359; auto x393 = x362 * x392; auto x394 = x106 * x131;
        auto x395 = x35 * x394; auto x396 = x283 * x57; auto x397 = u * x232; auto x398 = x34 * x53;
        auto x399 = x36 * x398; auto x400 = AXG * x399; auto x401 = x227 * x234;
        auto x402 = x283 * x401; auto x403 = x103 * x3; auto x404 = x15 * x403;
        auto x405 = x157 * x44; auto x406 = x29 * x403; auto x407 = x119 * x202;
        auto x408 = x36 * x84; auto x409 = x284 * x408; auto x410 = x166 * x28;
        auto x411 = Mqi * x36; auto x412 = x119 * x36; auto x413 = x412 * x44;
        auto x414 = x14 * x412; auto x415 = AXG * x413; auto x416 = x124 * x200;
        auto x417 = x244 * x78; auto x418 = x78 * x95;
        auto x419 = s * x417 + s * x418 - x200 * x274 + x200 * x278; auto x420 = kR * syFC15;
        auto x421 = iLGp * x339; auto x422 = x14 * x78; auto x423 = x136 * x422;
        auto x424 = s * x23; auto x425 = x187 * x424; auto x426 = x278 * x78; auto x427 = s * x426;
        auto x428 = x300 * x78; auto x429 = x302 * x78; auto x430 = x264 * x429;
        auto x431 = x38 * x408; auto x432 = u * x431; auto x433 = 4.0 * x78; auto x434 = x33 * x433;
        auto x435 = MXs * x434; auto x436 = x200 * x33; auto x437 = x176 * x33;
        auto x438 = x200 * x242; auto x439 = x200 * x66; auto x440 = MXs * x339;
        auto x441 = iLGp * x440; auto x442 = MGs * UU; auto x443 = x0 * x36; auto x444 = x443 * x84;
        auto x445 = x190 * x78; auto x446 = x203 * x84; auto x447 = Mqi * UU * x34;
        auto x448 = MGs * x78; auto x449 = x195 * x200; auto x450 = 2.0 * x200;
        auto x451 = x175 * x78; auto x452 = x201 * x57; auto x453 = x452 * x84;
        auto x454 = -u * x418; auto x455 = x359 * pow(MX, 7); auto x456 = x455 * x84;
        auto x457 = AXG * x130 * x228; auto x458 = MXs * x78; auto x459 = x129 * x46;
        auto x460 = x236 * x46; auto x461 = x201 * x84; auto x462 = x218 * x461;
        auto x463 = x114 * x462; auto x464 = x124 * x462; auto x465 = iLGp * x53;
        auto x466 = x465 * x59; auto x467 = 8.0 * x466; auto x468 = x131 * x372;
        auto x469 = x468 * x84; auto x470 = x245 * x78; auto x471 = kR * syFC8;
        auto x472 = 3.0 * x130; auto x473 = x422 * x472; auto x474 = 5.0 * x219;
        auto x475 = x367 * x78; auto x476 = 3.0 * x220; auto x477 = x412 * x57;
        auto x478 = 3.0 * x212; auto x479 = 3.0 * x241; auto x480 = x131 * x412;
        auto x481 = AXG * x185; auto x482 = AXG * x225; auto x483 = 3.0 * x268;
        auto x484 = 3.0 * x242; auto x485 = 5.0 * x66; auto x486 = MGs * syFC4;
        auto x487 = x108 * x253; auto x488 = x119 * x254; auto x489 = x37 * x52;
        auto x490 = x35 * x489; auto x491 = x42 * x489; auto x492 = x108 * x36 * x52;
        auto x493 = MXs * x1; auto x494 = kR * x412 * x52; auto x495 = MX * x1;
        auto x496 = x108 * x495; auto x497 = 2.0 * x496; auto x498 = x120 * x495;
        auto x499 = 2.0 * x498; auto x500 = x23 * x5 * x52; auto x501 = x1 * x5;
        auto x502 = 2.0 * MXs; auto x503 = 2.0 * x124;
        auto x504 = AXG * x249 + AXG * x500 - x249 - x500 + x501 * x502 - x501 * x503;
        auto x505 = UU * syFC11; auto x506 = 2.0 * x201; auto x507 = s * x506;
        auto x508 = kL * x507; auto x509 = 2.0 * x28; auto x510 = s * x509; auto x511 = kL * x269;
        auto x512 = x259 * x294; auto x513 = x264 * x506; auto x514 = x264 * x509;
        auto x515 = 4.0 * x49; auto x516 = x37 * x515; auto x517 = x47 * x59; auto x518 = u * x264;
        auto x519 = 4.0 * x518; auto x520 = x266 * x59; auto x521 = UU * syFC12;
        auto x522 = x10 * x115; auto x523 = MX * x346; auto x524 = x339 * x523;
        auto x525 = x49 * x94; auto x526 = x175 * x64; auto x527 = x340 * x41;
        auto x528 = UU * x527; auto x529 = s * x62; auto x530 = 4.0 * x321; auto x531 = 4.0 * x64;
        auto x532 = x264 * x62; auto x533 = 6.0 * x61; auto x534 = x19 * x24;
        auto x535 = x195 * x64 - x293 + x297 - x304 - x534 * x61 - x534 * x64;
        auto x536 = kL * syFC14; auto x537 = x35 * x455; auto x538 = u * x10;
        auto x539 = x538 * x95; auto x540 = x201 * x42; auto x541 = AXG * x131;
        auto x542 = x131 * x275; auto x543 = x10 * x227; auto x544 = x20 * x543;
        auto x545 = x218 * x540; auto x546 = x275 * x284; auto x547 = x236 * x7;
        auto x548 = x10 * x248; auto x549 = UU * x548; auto x550 = x10 * x247;
        auto x551 = x253 * x3; auto x552 = x41 * x551; auto x553 = UU * x552;
        auto x554 = x109 * x42; auto x555 = AXG * x64; auto x556 = kL * syFC6;
        auto x557 = x52 * x87; auto x558 = SS * x270; auto x559 = UU * x269; auto x560 = MXs * x52;
        auto x561 = MXs * x522; auto x562 = MXs * x320; auto x563 = x52 * x55;
        auto x564 = AXG * x305; auto x565 = x10 * x175; auto x566 = x10 * x62;
        auto x567 = x378 * x94; auto x568 = 8.0 * x235; auto x569 = x538 * x94;
        auto x570 = x538 * x57; auto x571 = MXs * x139; auto x572 = UU * x59;
        auto x573 = x234 * x572; auto x574 = x236 * x237; auto x575 = SS * x157;
        auto x576 = x122 * x52; auto x577 = x126 * x52; auto x578 = x38 * x412;
        auto x579 = x176 * x57; auto x580 = x112 * x566; auto x581 = 2.0 * x403;
        auto x582 = x112 * x559; auto x583 = x200 * x24; auto x584 = x112 * x26;
        auto x585 = x114 * x126; auto x586 = kR * syFC14; auto x587 = x377 * x78;
        auto x588 = x379 * x587; auto x589 = SS * x78; auto x590 = x144 * x589;
        auto x591 = x377 * x590; auto x592 = x587 * x93; auto x593 = x139 * x592;
        auto x594 = x131 * x587; auto x595 = MXs * x594; auto x596 = x157 * x227;
        auto x597 = x119 * x201; auto x598 = x218 * x597; auto x599 = x131 * x408;
        auto x600 = MXs * x599 - x146 * x458 + x221 * x403 + x364 * x408; auto x601 = AXG * x592;
        auto x602 = x26 * x411; auto x603 = x408 * x94; auto x604 = x359 * x461;
        auto x605 = x232 * x78; auto x606 = x458 * x7; auto x607 = x237 * x78;
        auto x608 = AXG * x80; auto x609 = x151 * x46; auto x610 = AXG * x607;
        auto x611 = x124 * x78;
        auto x612 = x190 * (x11 * x241 + x11 * x268 + x11 * x50 - x138 * x266 - x151 * x575 +
                            x151 * x606 + x151 * x607 - x151 * x608 + x151 * x610 -
                            x151 * x611 * x7 + x163 - x17 - x178 + x191 - x213 + x239 * x264 +
                            x243 - x258 + x268 * x5 - x458 * x609 + x58 * x7 + x6 + x609 * x611);
        auto x613 = x37 * x54; auto x614 = MX * x398; auto x615 = MUs * x614;
        auto x616 = MXs * x614; auto x617 = x108 * x329; auto x618 = MX * x7; auto x619 = kL * x618;
        auto x620 = x35 * x619; auto x621 = x42 * x619; auto x622 = MX * x227;
        auto x623 = x108 * x622; auto x624 = x112 * x543; auto x625 = AXG * x624;
        auto x626 = x114 * x614; auto x627 = x124 * x614; auto x628 = syFC13 * x52;
        auto x629 = x34 * x551; auto x630 = x35 * x495; auto x631 = 2.0 * x630;
        auto x632 = x42 * x495; auto x633 = 2.0 * x632; auto x634 = x294 * x52;
        auto x635 = -AXG * x548 - AXG * x634 + x1 * x124 * x323 + x548 + x634;
        auto x636 = x10 * x114; auto x637 = x172 * x424; auto x638 = x23 * x264;
        auto x639 = x201 * x398; auto x640 = 2.0 * x639; auto x641 = x506 * x54;
        auto x642 = 2.0 * x13; auto x643 = x114 * x87; auto x644 = x275 * x57;
        auto x645 = 2.0 * x139; auto x646 = AXG * x567; auto x647 = MXs * x399;
        auto x648 = x36 * x54; auto x649 = AXG * x54; auto x650 = x201 * x649;
        auto x651 = x114 * x399; auto x652 = 6.0 * x124; auto x653 = x114 * x648;
        auto x654 = AXG * x396; auto x655 = AXG * x644; auto x656 = s * x536;
        auto x657 = x175 * x283; auto x658 = 2.0 * x614; auto x659 = x59 * x658;
        auto x660 = x35 * x622; auto x661 = x208 * x572; auto x662 = 2.0 * x661;
        auto x663 = x208 * x227; auto x664 = u * x61; auto x665 = x131 * x664;
        auto x666 = x208 * x94; auto x667 = x131 * x208; auto x668 = AXG * x658;
        auto x669 = 2.0 * x660; auto x670 = 2.0 * x227; auto x671 = x208 * x670;
        auto x672 = MXs * x190; auto x673 = -u * x667 + x543 * x672; auto x674 = x175 * x378;
        auto x675 = AXG * x674; auto x676 = x318 * x46; auto x677 = u * x218;
        auto x678 = x275 * x677; auto x679 = 2.0 * x543; auto x680 = x59 * x614;
        auto x681 = 3.0 * MXs; auto x682 = x124 * x190; auto x683 = AXG * x680;
        auto x684 = s * x358; auto x685 = AXG * x679; auto x686 = x172 * x649;
        auto x687 = x42 * x622; auto x688 = AXG * x687; auto x689 = u * x272;
        auto x690 = AXG * x660; auto x691 = AXG * x235; auto x692 = -x235 - x291 + x691;
        auto x693 = 4.0 * SS; auto x694 = AXG * x93; auto x695 = x10 * x148; auto x696 = MX * x364;
        auto x697 = x35 * x696; auto x698 = 8.0 * x13; auto x699 = MGs * x543;
        auto x700 = AXG * x699; auto x701 = x28 * x451; auto x702 = 2.0 * x596;
        auto x703 = x201 * x85; auto x704 = 2.0 * x703; auto x705 = x377 * x451;
        auto x706 = x223 * x78; auto x707 = x224 * x78; auto x708 = x218 * x226;
        auto x709 = x708 * x78; auto x710 = u * x709; auto x711 = x397 * x78;
        auto x712 = MXs * x188; auto x713 = x411 * x566; auto x714 = AXG * x703;
        auto x715 = x466 * x503; auto x716 = 2.0 * x112; auto x717 = x112 * x679;
        auto x718 = u * x57; auto x719 = x200 * x718; auto x720 = -x719; auto x721 = x176 * x718;
        auto x722 = x112 * x691; auto x723 = s * x586; auto x724 = x190 * x575;
        auto x725 = x28 * x465; auto x726 = 2.0 * x725; auto x727 = x465 * x642;
        auto x728 = x190 * x465; auto x729 = x193 * x364; auto x730 = Mqi * x398;
        auto x731 = x509 * x730; auto x732 = x200 * x364; auto x733 = x642 * x730;
        auto x734 = x103 * x53; auto x735 = x502 * x734; auto x736 = MXs * x448;
        auto x737 = x445 * x7; auto x738 = x106 * x88; auto x739 = MQis * x176;
        auto x740 = MXs * x728; auto x741 = MXs * x259 * x730; auto x742 = x119 * x618;
        auto x743 = x126 * x502; auto x744 = 2.0 * MQis; auto x745 = x124 * x126;
        auto x746 = syFC2 * x151; auto x747 = -x445 * x76; auto x748 = x227 * x445;
        auto x749 = MXs * x748; auto x750 = x744 * x78; auto x751 = AXG * x466;
        auto x752 = x124 * x227; auto x753 = x59 * x730; auto x754 = x122 * x59;
        auto x755 = x622 * x84; auto x756 = x119 * x329; auto x757 = x445 * x46;
        auto x758 = x76 * x78; auto x759 = x200 * x76; auto x760 = MXs * x755;
        auto x761 = x193 * x81; auto x762 = x200 * x81; auto x763 = x200 * x94;
        auto x764 = x124 * x755; auto x765 = s * x471; auto x766 = x201 * x693;
        auto x767 = x372 * x93; auto x768 = x236 * x377 * x589; auto x769 = x114 * x89;
        auto x770 = x587 * x94; auto x771 = x408 * x57; auto x772 = AXG * x771;
        auto x773 = AXG * x477 - x770 + x772; auto x774 = -x599; auto x775 = x131 * x200;
        auto x776 = MXs * x775; auto x777 = x114 * x86; auto x778 = x732 + x774 + x776 - 8.0 * x777;
        auto x779 = SS * x193; auto x780 = x226 * x730; auto x781 = x114 * x730;
        auto x782 = x119 * x443; auto x783 = AXG * x639; auto x784 = x61 * x718;
        auto x785 = x55 * x59; auto x786 = x275 * x62; auto x787 = MX * x649;
        auto x788 = x59 * x787; auto x789 = x64 * x718; auto x790 = x10 * x223;
        auto x791 = x10 * x224; auto x792 = -x538 * x708;
        auto x793 = x10 * x397 - x28 * x565 + x674 + x790 + x791 + x792;
        auto x794 = MXs * x570 - x235 * x503 + x786 + x789; auto x795 = MUs * x399;
        auto x796 = u * x131; auto x797 = 5.0 * x124;
        auto x798 = x124 * x399 + x124 * x689 + x651 - x783; auto x799 = x190 * x87;
        auto x800 = -x697; auto x801 = MQis * x55; auto x802 = MXs * x347; auto x803 = x472 * x87;
        auto x804 = 3.0 * x13; auto x805 = 6.0 * x13; auto x806 = x218 * x59;
        auto x807 = AXG * x175; auto x808 = MQis * x543; auto x809 = x20 * x465;
        auto x810 = AXG * x770; auto x811 = x372 * x88; auto x812 = x104 * x36;
        auto x813 = x502 * x86 + 2.0 * x812; auto x814 = MXs * x86; auto x815 = x218 * x412;
        auto x816 = x124 * x225 + x124 * x86 + x706 + x707 - x714 + x777; auto x817 = x131 * x193;
        auto x818 = syFC6 * x151; auto x819 = MUs * x730;
        auto x820 = -x126 * x265 - x126 * x503 + x176 * x541 - x579; auto x821 = x408 * x693;
        auto x822 = x114 * x218; auto x823 = x329 * x35; auto x824 = x193 * x24;
        auto x825 = x408 * x46; auto x826 = MXs * x193; auto x827 = x200 * x237;
        auto x828 = x200 * x7; auto x829 = -x221; auto x830 = -x222; auto x831 = -x229;
        auto x832 = -x573; auto x833 = -x148 - x244 + x245 - x95;
        auto x834 = -x274 + x278 - x300 + x302 + x352 + x354 - x355 - x357 + x397;
        auto x835 = x15 + x45; auto x836 = x169 - x268 + x834 + x835;
        auto x837 =
            x146 + x161 + x223 + x224 + x231 + x233 + x574 + x829 + x830 + x831 + x832 + x834;
        auto x838 = x11 * x486;
        auto x839 = x205 * (x16 - x2 + x242 + x321 - x33 - x334 + x48 + x66 - x8);
        ,
        Nc * TR * TR * gs * gs * (Nc - 1) * (Nc + 1) *
            (-kL * syFC1 *
                 (-AXG * x331 * x44 + AXG * x524 + s * x526 + x10 * x336 + x10 * x525 + x111 -
                  x29 * x323 - x317 + x320 - x323 * x60 - x323 * x66 + x324 - x332 + x337 + x345 -
                  x522 - x524) +
             kL * syFC10 * x0 * (-x566 - x658 + x668 + x679 + x685 + x686) +
             kL * syFC3 *
                 (-AXG * x342 - MGs * x111 + MGs * x345 - x109 * x35 - x115 * x61 + x141 * x61 -
                  x2 * x208 + x208 * x30 + x30 * x61 + x342 + x343 * x344 + x343 * x347 -
                  x344 * x348 - x347 * x348 + x349) +
             kL * x505 *
                 (-AXG * x552 - AXG * x629 - AXG * x631 - AXG * x633 + x552 + x629 + x631 + x633 +
                  x635) +
             kL * x521 *
                 (-AXG * x341 - AXG * x527 + AXG * x630 + AXG * x632 - x1 * x269 + x1 * x636 +
                  x264 * x512 + x269 * x343 - x323 * x493 + x327 * x638 + x331 * x638 + x341 -
                  x343 * x636 - x35 * x637 - x42 * x637 + x527 - x630 - x632 + x635) +
             kL * x628 *
                 (-AXG * x238 + AXG * x330 + AXG * x823 + x238 - x25 + x289 - x330 + x386 + x399 -
                  x400 + x660 + x687 - x688 - x690 + x692 - x823) -
             kR * syFC10 *
                 (AXG * x431 - AXG * x576 - AXG * x577 + AXG * x578 + s * x200 * x57 + s * x579 -
                  s * x580 - x15 * x433 + x176 * x70 - x176 * x74 + x187 * x47 + x19 * x575 +
                  x19 * x80 - x200 * x529 + x200 * x70 - x200 * x74 + x29 * x433 - x431 -
                  x433 * x45 - x433 * x48 + x433 * x60 + x433 * x66 - x433 * x8 - x434 - x51 * x78 +
                  x576 + x577 - x578) +
             kR * syFC3 *
                 (-AXG * x432 + AXG * x446 + UU * x444 - x124 * x339 * x447 - x124 * x421 * x442 +
                  x141 * x200 - x193 * x195 - x193 * x207 + 2.0 * x193 * x212 + x193 * x215 +
                  x193 * x216 - x193 * x217 + x212 * x445 + x216 * x448 - x217 * x448 - x29 * x445 -
                  x32 * x445 + x33 * x445 + x33 * x450 + x432 + x440 * x447 + x441 * x442 - x446 -
                  x449) -
             kR * x521 *
                 (AXG * x444 + AXG * x782 - x103 * x440 + x114 * x441 + x119 * x507 - x119 * x513 -
                  x264 * x294 * x716 + x295 * x716 + x403 * x510 - x403 * x514 - x408 * x515 +
                  x408 * x519 - x412 * x515 + x412 * x519 + x424 * x450 + x424 * x581 -
                  x433 * x517 + x433 * x520 - x444 - x450 * x638 + x507 * x84 - x513 * x84 -
                  x581 * x638 - x782) +
             kR * x628 *
                 (-AXG * x742 + AXG * x756 - AXG * x828 - x101 - x105 - x107 + x123 + x125 - x466 -
                  x582 + x585 + x607 - x610 + x624 - x738 + x742 + x745 + x751 - x756 + x828 + x86 +
                  x89 - x99) -
             s * x288 *
                 (-AXG * x544 + MXs * x687 + MXs * x689 - x124 * x687 - x235 * x681 - x235 * x797 +
                  x544 + x639 + x64 * x796 - x647 - x650 + x665 - x678 + x785 - x788 + x793 - x795 +
                  x798) -
             s * x338 *
                 (MGs * x309 * x57 - MQis * x124 * x679 - MQis * x509 * x87 + MQis * x526 +
                  MQis * x569 + MQis * x570 - x124 * x326 - x139 * x799 + x190 * x42 * x618 -
                  x190 * x558 - x208 * x364 + x226 * x799 - x235 * x744 + x259 * x616 +
                  x259 * x802 - x347 * x509 + x347 * x642 + x502 * x801 - x503 * x801 - x55 * x672 +
                  x55 * x682 + x614 * x642 - x640 + x643 * x672 + x800) -
             s * x420 *
                 (-AXG * x596 + MXs * x225 + MXs * x624 + u * x775 - u * x815 + x112 * x235 -
                  x124 * x624 + x176 * x796 - x466 * x681 - x466 * x797 + x596 - x701 + x703 +
                  x705 - x710 + x711 - x722 - x811 - x812 - x814 + x816) +
             s * x556 *
                 (-AXG * x640 + AXG * x657 - AXG * x662 - MGs * x570 - u * x666 - x124 * x669 -
                  x124 * x671 + x131 * x61 * x73 + x190 * x235 + x190 * x25 - x190 * x289 +
                  x198 * x208 - x347 * x645 + x502 * x660 + x502 * x663 - x59 * x668 + x640 - x657 +
                  x659 + x662 - x665 + x666 * x73 + x667 * x73 + x673) -
             syFC1 * (-AXG * x174 - x127 - x131 * x177 - x151 * x158 - x151 * x167 + x154 +
                      2.0 * x155 + x156 - 2.0 * x159 - x162 + 2.0 * x163 + 2.0 * x164 + x165 -
                      2.0 * x168 - x170 + x171 + x173 * x72 + x174 + x175 * x177 - 2.0 * x178 -
                      x180 - x181 * x182 + x182 * x183 + x189) +
             syFC10 *
                 (-AXG * x40 - AXG * x43 + x11 * x51 + x12 * x15 - x12 * x29 + x12 * x33 +
                  x12 * x45 + x12 * x48 - x12 * x60 + x12 * x8 - x18 - x19 * x22 - x26 * x27 +
                  2.0 * x31 + x40 + x43 - x52 * x56 - x57 * x58 + 2.0 * x6 + x61 * x63 - x61 * x68 +
                  x63 * x64 - x64 * x68 - x67 - x69 * x71 + x69 * x75 - x71 * x72 + x72 * x75) +
             syFC14 * (AXG * x98 + MUs * x18 + MUs * x31 - MUs * x6 + x100 + x102 + x105 * x77 +
                       x107 * x77 + 3.0 * x110 + 3.0 * x113 - x114 * x116 + x118 + x121 -
                       x123 * x77 - x124 * x127 - x125 * x77 + x52 * x92 + 5.0 * x76 * x79 -
                       4.0 * x77 * x80 - x77 * x86 - x82 - x90 + x95 * x96 - x98) -
             syFC2 * (MGs * x127 - MGs * x171 - MQis * x116 - MQis * x127 + MQis * x142 +
                      MQis * x156 - MQis * x162 + MQis * x165 - MQis * x170 + MQis * x171 -
                      MQis * x180 + MQis * x31 + MQis * x6 + x108 * x197 - x115 * x194 +
                      3.0 * x134 + x141 * x194 + x150 - x151 * x198 * x200 - x151 * x199 -
                      x17 * x190 + x190 * x191 + x192 + x194 * x195 + 2.0 * x196) -
             syFC3 * (-AXG * x204 - MGs * x12 * x65 + MGs * x67 + MX * x108 * x115 + kL * x211 -
                      u * x40 + x100 + x110 - x118 + x134 - x190 * x213 + x192 - x194 * x30 + x196 +
                      x197 * x69 + x204 + x205 * x32 - x205 * x33 + x206 - x209 * x215 -
                      x209 * x216 + x209 * x217 + x210 - x212 * x214 + x40 * x73) +
             syFC5 * x11 *
                 (-x131 * x49 + x2 + x212 - x230 + x240 - x241 + x242 - x29 - x30 + x305 - x32 -
                  x525 - x60 + x65 + x836) +
             syFC5 * x5 *
                 (x146 + x223 + x224 + x231 + x233 + x574 + x829 + x830 + x831 + x832 + x833 +
                  x836) -
             syFC5 * (kL * x234 * x235 + x11 * x221 + x11 * x222 - x11 * x223 - x11 * x224 +
                      x11 * x229 - x11 * x231 - x11 * x233 - x147 + x151 * x225 + x154 + x186 +
                      x230 * x5 - x236 * x239 - x240 * x5 + x243 + x246 + x31 - x6) -
             syFC6 * x77 *
                 (x153 - x158 - x167 + x176 * x392 - x176 * x677 - x176 * x807 + x225 - x481 +
                  x482 + x580 + 2.0 * x582 - x611 * x670 - x726 + x727 + x735 + x737 + x743 + x748 -
                  2.0 * x751 - x757 + x820) -
             syFC6 *
                 (AXG * x250 + AXG * x255 - AXG * x262 - s * x263 + x131 * x251 + x15 * x205 +
                  x15 * x214 - x169 * x209 - x175 * x251 - x190 * x258 + x205 * x45 - x206 +
                  x209 * x267 - x210 - x214 * x241 - x214 * x268 + x214 * x45 + x250 - x252 * x52 -
                  x255 + x257 * x79 + x259 * x261 - x261 * x265 + x262 + x263 * x264) -
             syFC7 * x612 + syFC7 * x839 +
             syFC8 * (MUs * x116 + MXs * x116 + MXs * x129 * x138 - MXs * x142 + MXs * x147 -
                      MXs * x149 - MXs * x31 - MXs * x6 - x102 - x11 * x130 * x131 -
                      x11 * x132 * x94 - x11 * x140 * x20 + x113 + x121 + x128 * x79 + x129 * x22 -
                      x133 * x22 - x134 + x137 * x5 + x143 - x150 - x52 * x91 + x82 + x90) -
             syFC9 * x612 + syFC9 * x839 +
             u * x288 *
                 (x10 * x146 - x10 * x245 - x234 * x648 + x236 * x380 + x275 * x94 + x281 +
                  x359 * x540 - x375 + x381 + x394 * x42 + x395 + x42 * x696 - x540 * x693 +
                  x540 * x694 - x542 - x567 - x643 * x698 - 8.0 * x651 - 8.0 * x653 + x654 + x655 +
                  x695 + x697) +
             u * x420 *
                 (-x114 * x465 * x698 + x119 * x394 + x119 * x696 - x119 * x766 + x119 * x767 +
                  x146 * x78 + x148 * x78 - x234 * x89 + x359 * x597 + x412 * x94 + x417 - x470 -
                  x480 + x594 + x768 - 8.0 * x769 + x773 + x778) +
             u * x471 *
                 (x140 * x193 - x144 * x779 - x184 * x220 + x193 * x232 + x193 * x379 +
                  x193 * x547 - x221 * x78 - x234 * x86 - x367 * x781 - x594 + x603 + x604 + x729 -
                  x766 * x84 + x767 * x84 - x768 + x770 + x772 + x778 - 8.0 * x780) +
             u * x765 *
                 (x105 + x123 - x124 * x728 + x124 * x730 + x139 * x184 - x359 * x412 -
                  x412 * x693 + x412 * x694 - x465 * x645 + x465 * x805 + x589 * x672 + x601 +
                  x681 * x730 + x709 - 6.0 * x725 + x740 + x781 + x815 - x817 - x819 + x820 -
                  6.0 * x828) -
             x0 * x338 *
                 (AXG * x808 + MGs * x55 - MGs * x787 - MQis * x565 + MQis * x787 - x175 * x208 +
                  x326 + x667 + x669 + x671 - x676 + 3.0 * x699 - x700 - x801 + x808) +
             x0 * x358 * (x290 + x291 - x383 - x669 + x688 - x699 + x700) +
             x0 * x536 *
                 (-AXG * x689 - x26 + 5.0 * x289 + x36 * x649 + x384 - x399 + x400 + x570 + x615 +
                  x616 - x626 - x627 + 3.0 * x660 + 3.0 * x687 + x688 + x690 + x692) -
             x23 * x471 *
                 (AXG * x158 - AXG * x590 + AXG * x775 + AXG * x817 + AXG * x821 - x133 * x589 -
                  x193 * x822 - x193 * x94 - x200 * x822 - x218 * x416 + 4.0 * x226 * x589 -
                  x234 * x730 + x236 * x779 - x359 * x408 + x590 - x601 + x605 - x709 + x763 +
                  x775 + x817 - x821) +
             x288 * (-AXG * x271 - AXG * x277 + MUs * x276 - s * x273 - x114 * x276 + x135 * x272 +
                     x15 * x269 - x202 * x42 + x212 * x269 - x241 * x269 - x269 * x274 -
                     x269 * x29 + x271 + x277 + x279 * x42 - x280 * x42 + x287) -
             x288 * (MXs * x542 - x114 * x545 - x124 * x545 - x146 * x309 - x232 * x270 +
                     x234 * x544 + x275 * x364 - x363 + x365 + x366 - x370 - x371 + x374 + x376 +
                     x382 + x42 * x455 + x537 - x539 + x540 * x541 - x540 * x57) -
             x288 * (AXG * x546 - MXs * x301 - x114 * x289 * x367 - x257 * x270 - x257 * x275 -
                     x286 * x64 + x294 * x547 + x300 * x64 + x302 * x64 + x350 + x351 - x352 * x64 -
                     x353 + x355 * x64 + x356 + x385 + x388 + x389 - x390 - x391 - x393 + x402 -
                     x546) +
             x288 * (MUs * x111 + x0 * x290 + x0 * x291 + x242 * x61 - x25 * x52 + x264 * x301 -
                     x269 * x298 - x269 * x313 + x292 + x293 + x296 - x297 - x299 + x304 +
                     x306 * x42 + x307 + x308 * x38 - x311 - x312 * x35 - x312 * x42 - x315 - x316 -
                     x60 * x61 - x61 * x65 + x61 * x66) +
             x338 * (4.0 * MGs * x10 * x48 + MQis * x317 - MQis * x320 - MQis * x324 + MQis * x332 -
                     MQis * x337 + s * x190 * x330 - x181 * x335 + x183 * x335 -
                     x190 * x238 * x264 - x208 * x336 - x211 + x242 * x318 + x307 + x318 * x321 -
                     x318 * x334 + x318 * x65 - x319 + x322 - x325 - x326 * x49 + x327 * x60 +
                     x328 + x333 * x60 + x333 * x65 - x336 * x61) +
             x358 * (-x208 * x274 + x208 * x278 - x208 * x300 + x208 * x302 + x208 * x352 +
                     x208 * x354 - x208 * x355 - x208 * x357 + x264 * x281 + x287 + x303 * x95 -
                     x350 - x351 + x353 - x356) -
             x358 * (AXG * x528 + MXs * x111 + MXs * x345 + x208 * x298 + x208 * x33 - x282 * x302 -
                     x292 + 6.0 * x310 + x315 - 4.0 * x316 - x318 * x33 + x319 + 5.0 * x32 * x61 -
                     x322 - x528 + x535 - x554 - x561 + x562 + x564 * x64) +
             x358 *
                 (x10 * x13 * x232 + x10 * x360 - x10 * x361 + x146 * x208 - x148 * x208 -
                  x208 * x221 + x208 * x233 - x208 * x244 + x208 * x245 - x208 * x95 + x220 * x272 -
                  x220 * x368 + x363 - x365 - x366 + x370 + x371 - x374 - x376 + x382) +
             x358 * (-AXG * x388 - u * x375 + u * x395 + x13 * x26 - x144 * x25 - x144 * x291 +
                     x144 * x386 - x208 * x222 + x208 * x224 - x208 * x229 + x208 * x231 +
                     x208 * x397 + x223 * x61 - x226 * x26 + x226 * x387 + x385 - x389 + x390 +
                     x391 + x393 + x396 * x73 - 8.0 * x400 * x59 - x402) -
             x358 * (AXG * x245 * x538 - u * x273 + u * x281 + u * x381 - u * x567 + x10 * x457 -
                     x130 * x565 - x132 * x565 - x139 * x568 - x139 * x569 - x208 * x223 +
                     x208 * x573 - x208 * x574 - x219 * x566 + x220 * x565 + x221 * x538 +
                     x226 * x568 - x226 * x570 - x270 * x460 + x309 * x459 + x378 * x460 + x537 +
                     x539 - 12.0 * x543 * x571) -
             x420 * (MXs * x480 - x114 * x598 - x119 * x452 + x119 * x455 + x119 * x468 -
                     x124 * x598 - x157 * x232 + x234 * x596 + x364 * x412 - x453 + x454 + x456 -
                     x463 - x464 + x469 + x588 - x591 + x593 + x595 + x600) +
             x420 * (-AXG * x405 - AXG * x409 - s * x410 - x114 * x414 + x135 * x166 +
                     x14 * x269 * x411 - x176 * x274 + x176 * x278 + x212 * x403 - x241 * x403 +
                     x256 * x416 + x257 * x408 - x274 * x403 + x404 + x405 - x406 - x407 + x409 +
                     x413 - x415 + x419) -
             x420 * (AXG * x602 - MXs * x428 - u * x603 - u * x604 - x114 * x475 * x76 -
                     x157 * x257 + x158 * x23 - x175 * x372 * x84 - x176 * x286 + x176 * x300 +
                     x176 * x302 - x176 * x352 + x176 * x355 + x200 * x300 + x200 * x302 -
                     x200 * x352 + x200 * x355 + x23 * x547 * x78 + x23 * x601 - x257 * x412 +
                     x401 * x408 + x461 * x62 - x602) +
             x420 * (x115 * x403 + x128 * x421 - x176 * x298 + x176 * x305 + x2 * x403 -
                     x200 * x298 + x200 * x305 - x200 * x60 - x200 * x65 + x264 * x428 -
                     x298 * x403 - x30 * x403 - x313 * x403 + x421 * x81 + x423 + x425 - x427 +
                     x430 + x432 - x435 - x436 - x437 + x438 + x439 - x52 * x80) -
             x471 * (-x129 * x575 - x129 * x606 - x13 * x605 + x130 * x166 + x132 * x152 +
                     x133 * x575 + x139 * x158 + x148 * x458 - x166 * x220 + x193 * x244 -
                     x193 * x245 + x193 * x95 + x220 * x367 * x465 - x360 * x78 + x361 * x78 -
                     x588 + x591 - x593 - x595 + x600) +
             x471 * (-AXG * x473 + s * x477 - x193 * x274 + x193 * x278 + 6.0 * x212 * x458 -
                     x264 * x410 + x264 * x417 + x264 * x418 - x264 * x480 + x265 * x414 -
                     x29 * x475 + x403 * x478 - x403 * x479 + 3.0 * x404 - 5.0 * x406 - 2.0 * x407 +
                     2.0 * x415 + x419 + x422 * x474 - x422 * x476 + x473) +
             x471 * (u * x410 - u * x417 + x130 * x451 + x132 * x451 + x139 * x153 +
                     12.0 * x139 * x227 * x458 + x139 * x467 + x146 * x193 - x148 * x193 +
                     x157 * x460 + x188 * x226 - x193 * x221 + x193 * x233 + x219 * x62 * x78 -
                     x226 * x467 + x453 + x454 - x456 - x457 * x78 - x458 * x459 + x463 + x464 -
                     x469 - x470 * x73) +
             x486 * x5 * (x833 + x837) -
             x486 *
                 (AXG * x154 - x124 * x151 * x166 - x151 * x481 + x151 * x482 - x155 - x164 + x189 +
                  x246 + x313 * x5 - x478 * x5 + x479 * x5 + x483 * x5 - x484 * x5 - x485 * x5) -
             x505 * (AXG * x487 + AXG * x488 - AXG * x490 - AXG * x491 - AXG * x492 - AXG * x494 +
                     AXG * x497 + AXG * x499 - x487 - x488 + x490 + x491 + x492 + x493 * x97 +
                     x494 - x497 - x499 + x504) -
             x521 * (-AXG * x496 - AXG * x498 + MUs * x501 - x114 * x501 - x12 * x517 + x12 * x520 +
                     x260 * x519 + x27 * x512 + x35 * x508 - x35 * x516 + x37 * x42 * x519 +
                     x42 * x508 - x42 * x516 + x496 + x498 + x504 + x510 * x511 - x511 * x514 -
                     x513 * x69 - x513 * x72) -
             x536 *
                 (-MUs * x524 + x114 * x524 + x124 * x524 + x195 * x269 + x215 * x269 +
                  x283 * x529 - x283 * x532 - x296 + x311 - 6.0 * x315 - 6.0 * x316 + x321 * x531 +
                  x325 - x328 - x440 * x523 + x528 + x530 * x61 + x533 * x66 + x535) -
             x556 * (AXG * x0 * x387 + AXG * x549 + AXG * x553 + x111 * x190 - x115 * x555 +
                     x131 * x550 - x16 * x318 + x16 * x531 - x175 * x550 - x289 * x52 + x30 * x531 -
                     x386 * x52 + x549 - x553 - 4.0 * x554) -
             x556 * (MUs * x563 - MXs * x30 * x323 + x114 * x522 - x114 * x563 - x124 * x563 +
                     x13 * x557 + x141 * x555 - x208 * x215 + x208 * x564 + x211 - x28 * x557 +
                     x292 - 2.0 * x299 - x307 + x318 * x8 + x325 - x349 - x52 * x558 - x531 * x8 +
                     x55 * x560 + x559 * x560 + x561 - x562 + x564 * x61) -
             x586 * (-s * x584 - x0 * x582 + x0 * x585 - x126 * x343 + x176 * x195 + x176 * x530 -
                     x19 * x583 + x200 * x530 + x32 * x581 + x33 * x581 + x408 * x529 -
                     x408 * x532 - x423 - x425 + x427 - x430 + x435 - 6.0 * x436 - 6.0 * x437 +
                     6.0 * x438 + 6.0 * x439 + x449 - x450 * x60 - x450 * x65) -
             x628 * (AXG * x252 + AXG * x613 - AXG * x617 + AXG * x620 + AXG * x621 + AXG * x623 +
                     MUs * x56 + MXs * x56 + kL * x615 + kL * x616 - kL * x626 - kL * x627 +
                     kR * x625 - x114 * x56 - x124 * x56 - x252 - x613 + x617 - x620 - x621 - x623 +
                     x91 - x92) -
             x656 * (-AXG * x784 - AXG * x786 - AXG * x789 - x20 * x679 + x20 * x685 - x235 * x502 -
                     x502 * x687 + x652 * x687 - x659 - x675 + 6.0 * x683 + 6.0 * x783 + x784 -
                     2.0 * x785 + 6.0 * x788 + x793 + x794) +
             x656 * (-x259 * x399 - x259 * x648 - x273 + x281 + x396 + x399 * x652 - x502 * x648 +
                     x509 * x559 - x559 * x642 - x559 * x645 + x567 + x640 + x641 + x642 * x643 +
                     x644 - x646 - 2.0 * x647 + x648 * x652 - 6.0 * x650 + 6.0 * x651 + 6.0 * x653 -
                     x654 - x655) -
             x684 * (AXG * x790 + AXG * x791 - x139 * x347 - x179 * x73 - x275 * x807 + x28 * x347 -
                     x309 * x806 + x375 + x503 * x687 - x543 * x805 + x59 * x686 - x639 + x647 +
                     x790 + x791 + x792 + x794 + x795 + x798) -
             x684 * (AXG * x542 + AXG * x803 - MUs * x802 + MXs * x326 + x114 * x802 - x13 * x347 -
                     3.0 * x139 * x559 + x226 * x347 - x265 * x648 + x28 * x368 + 5.0 * x28 * x559 -
                     x395 - x474 * x87 + x476 * x87 - x503 * x648 - x559 * x804 - 6.0 * x571 * x87 +
                     x641 + x643 * x804 - x644 + x646 + x695 + x800 - x803) +
             x684 * (AXG * x661 + MXs * x660 + MXs * x676 + x10 * x231 + x10 * x355 + x124 * x660 +
                     x124 * x663 - x129 * x543 - x139 * x679 - x237 * x533 - x259 * x289 -
                     x265 * x289 - x275 * x392 + x283 * x677 + x301 - x543 * x682 + x657 - x661 +
                     x663 * x681 - x664 * x94 + x673 + x675 + x678 + x680 + x683) -
             x723 *
                 (-x114 * x727 + x410 + 2.0 * x411 * x559 - x417 - x418 - x477 + x502 * x809 +
                  x502 * x89 - x503 * x809 - x506 * x88 - x509 * x734 + x642 * x734 + x645 * x734 -
                  x652 * x89 - 6.0 * x769 - x771 + x773 - 6.0 * x777 + x810 + 6.0 * x811 + x813) +
             x723 * (-AXG * x702 + AXG * x705 + AXG * x713 + AXG * x719 + AXG * x721 + MXs * x717 +
                     x235 * x716 + x466 * x502 - x624 * x652 + x652 * x86 + x701 + x702 + x704 -
                     x705 - x706 - x707 + x710 - x711 - x712 - x713 - 6.0 * x714 + x715 + x720 -
                     x721 - 6.0 * x722) +
             x746 * (MQis * x158 + MQis * x167 + MQis * x726 - MQis * x727 - MQis * x735 -
                     MQis * x743 - x105 * x502 - x114 * x740 - x122 * x642 + x124 * x737 +
                     x131 * x739 + x139 * x728 + x190 * x738 - x190 * x742 - x190 * x745 -
                     x226 * x728 - x57 * x736 + x704 + x724 + x729 + x731 + x732 - x733 - x741 +
                     x744 * x745) +
             x746 * (-MQis * x112 * x685 - MQis * x153 + MQis * x185 - MQis * x188 + MQis * x717 -
                     x124 * x757 - x175 * x739 - x190 * x607 + x190 * x608 - x190 * x610 -
                     x190 * x624 + x190 * x625 + x190 * x751 + x190 * x756 + x193 * x227 * x502 -
                     x193 * x718 + x466 * x744 + x502 * x755 + x62 * x736 + x720 + x747 + x749 +
                     x750 * x752 - x750 * x81 + 2.0 * x753 + 2.0 * x754) -
             x765 * (AXG * x706 + AXG * x707 + MXs * x737 + MXs * x781 - MXs * x819 - x13 * x730 -
                     x139 * x730 + x28 * x730 - x458 * x806 + x599 - x703 + x712 - x715 - x732 -
                     x776 + x780 + x810 + x812 + x814 + x816) +
             x765 * (AXG * x429 - u * x763 - x112 * x384 - x128 * x193 + x175 * x408 + x176 * x257 +
                     x199 - x265 * x758 - x289 * x716 + x355 * x78 + x408 * x677 - x426 + x428 +
                     x429 + x445 * x81 - x581 * x76 + 4.0 * x583 + x584 - x652 * x758 + x747 -
                     5.0 * x759 + x760 - x761 - x762 + x764) -
             x818 * (AXG * x599 + AXG * x704 - MXs * x817 + u * x817 + x124 * x817 + x13 * x728 -
                     x190 * x466 - x190 * x725 - x265 * x86 + x448 * x718 - x502 * x781 -
                     x503 * x86 + x645 * x730 + x672 * x734 - x704 - x724 - x731 + x733 + x741 -
                     x749 + x774 - 2.0 * x780 + x813) +
             2.0 * x818 *
                 (AXG * x193 * x237 - AXG * x583 - AXG * x753 - AXG * x754 - AXG * x824 +
                  AXG * x825 + AXG * x827 + x124 * x193 * x46 - x193 * x752 - x193 * x76 +
                  x227 * x826 + x24 * x448 - x448 * x76 - x46 * x826 + x583 + x753 + x754 - x759 +
                  x760 + x761 + x762 - x764 + x824 - x825 - x827) +
             x838 * (x267 - x483 + x835 + x837) -
             x838 * (-x115 + x131 * x518 + x141 + x2 + x29 + x30 - x305 + x313 + x32 + x336 - x478 +
                     x479 - x484 - x485 + x518 * x94 - x564 + x60 - x65)) *
            Denom(768 * s * (MUs - u) * pow(Pi, 2)));

    //}
  }

  return ret.real();
}

#undef A0
#undef B0
#undef C0
#undef C1
#undef C2
#undef C00
#undef C11
#undef C12
#undef C22