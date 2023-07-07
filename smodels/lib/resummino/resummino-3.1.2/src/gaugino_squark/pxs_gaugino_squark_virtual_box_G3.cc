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

ComplexType ME_us_box_GQQq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;

  int itq = q;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int ftq = 0; ftq < 2; ftq++) {
    int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
    // int iq = (itq + ftq * 3) % 6;
    int iq = (sq - is_up_squark(sq) * 6) % 3 + is_up_squark(sq) * 3;

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


    //std::cout << "isq " << isq << " iq " << iq  << " q " << q << std::endl;

    auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C0(MXs, 0, MUs + MXs - s - u, Mqis, MQis, MQis)
#define syFC2 D0(MXs, 0, 0, MUs, MUs + MXs - s - u, s, Mqis, MQis, MQis, MGs)

#define syFC3 C1(MXs, MUs + MXs - s - u, 0, MQis, Mqis, MQis)
#define syFC4 D1(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC5 C2(MXs, MUs + MXs - s - u, 0, MQis, Mqis, MQis)
#define syFC6 D2(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC7 D3(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC8 D00(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC9 D12(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC10 D13(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC11 D22(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC12 D23(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

#define syFC13 D33(0, s, MXs, MUs + MXs - s - u, 0, MUs, MQis, MGs, MQis, Mqis)

    _EPS0_(
        ret, auto x0 = UU * u; auto x1 = pow(s, 3); auto x2 = R * jLGp; auto x3 = iRGp * x2;
        auto x4 = kL * x3; auto x5 = x1 * x4; auto x6 = AXG * x0; auto x7 = pow(s, 2);
        auto x8 = UU * x7; auto x9 = MXs * x8; auto x10 = 2.0 * x3; auto x11 = kL * x10;
        auto x12 = u * x11; auto x13 = pow(MXs, 2); auto x14 = L * jRGp; auto x15 = x13 * x14;
        auto x16 = iLGp * x15; auto x17 = kR * x16; auto x18 = MUs * SS; auto x19 = 4.0 * x18;
        auto x20 = s * x19; auto x21 = pow(MUs, 2); auto x22 = SS * x21; auto x23 = s * x22;
        auto x24 = 4.0 * MXs; auto x25 = iLGp * x14; auto x26 = kR * x25; auto x27 = x24 * x26;
        auto x28 = pow(u, 2); auto x29 = UU * x28; auto x30 = AXG * x29; auto x31 = x4 * x7;
        auto x32 = x30 * x31; auto x33 = 2.0 * x7; auto x34 = SS * x28; auto x35 = AXG * x34;
        auto x36 = x35 * x4; auto x37 = UU * s; auto x38 = x21 * x26;
        auto x39 = 2.0 * MXs * x37 * x38; auto x40 = SS * u; auto x41 = 4.0 * x40;
        auto x42 = iLGp * x7; auto x43 = x14 * x42; auto x44 = MX * Mqi; auto x45 = x43 * x44;
        auto x46 = kL * x41 * x45; auto x47 = 2.0 * s; auto x48 = MXs * x47; auto x49 = AXG * x26;
        auto x50 = AXG * x9; auto x51 = x14 * x41; auto x52 = MG * iRGp; auto x53 = x51 * x52;
        auto x54 = MX * x53; auto x55 = kL * x54 * x7; auto x56 = u * x8; auto x57 = Mqi * x25;
        auto x58 = MX * x57; auto x59 = kL * x58; auto x60 = x56 * x59; auto x61 = AXG * x56;
        auto x62 = x59 * x61; auto x63 = x14 * x52; auto x64 = MX * x63; auto x65 = kL * x64;
        auto x66 = x61 * x65; auto x67 = x29 * x4; auto x68 = 3.0 * kL; auto x69 = x64 * x68;
        auto x70 = x56 * x69 - x67 * x7; auto x71 = UU * x1; auto x72 = x64 * x71;
        auto x73 = kL * x72; auto x74 = x59 * x71; auto x75 = -AXG * x73 + AXG * x74 + x73 - x74;
        auto x76 = pow(MXs, 3); auto x77 = x25 * x76; auto x78 = 2.0 * SS; auto x79 = MUs * x78;
        auto x80 = kR * x79; auto x81 = pow(MUs, 3); auto x82 = MXs * x78; auto x83 = 4.0 * SS;
        auto x84 = pow(MX, 4); auto x85 = x25 * x84; auto x86 = 6.0 * x85; auto x87 = x4 * x56;
        auto x88 = MUs * x87; auto x89 = AXG * x84; auto x90 = x78 * x89; auto x91 = x21 * x78;
        auto x92 = pow(MX, 3); auto x93 = Mqi * iRGp; auto x94 = x2 * x93; auto x95 = x92 * x94;
        auto x96 = kR * x95; auto x97 = x4 * x61; auto x98 = MUs * x97; auto x99 = 2.0 * x14;
        auto x100 = x52 * x99; auto x101 = MX * x100; auto x102 = kL * x56; auto x103 = x33 * x40;
        auto x104 = x103 * x59; auto x105 = MXs * x95; auto x106 = x103 * x65;
        auto x107 = 2.0 * x57; auto x108 = MX * x107; auto x109 = MUs * x8; auto x110 = x50 * x69;
        auto x111 = 6.0 * x61; auto x112 = s * x29; auto x113 = x37 * x84; auto x114 = x13 * x79;
        auto x115 = x21 * x82; auto x116 = pow(MX, 6); auto x117 = 4.0 * UU; auto x118 = AXG * x117;
        auto x119 = x116 * x118; auto x120 = 4.0 * x84; auto x121 = x120 * x40; auto x122 = u * x78;
        auto x123 = x122 * x21; auto x124 = x122 * x13; auto x125 = x34 * x47;
        auto x126 = x120 * x18; auto x127 = AXG * x113; auto x128 = MUs * x37; auto x129 = u * x128;
        auto x130 = MXs * u; auto x131 = x130 * x37; auto x132 = x79 * x89; auto x133 = u * x18;
        auto x134 = x133 * x47; auto x135 = 8.0 * x6; auto x136 = AXG * x112;
        auto x137 = 3.0 * x136; auto x138 = x122 * x89; auto x139 = x118 * x84;
        auto x140 = MUs * x139; auto x141 = MXs * x139; auto x142 = AXG * x129;
        auto x143 = MUs * x0; auto x144 = AXG * x143; auto x145 = 8.0 * MXs;
        auto x146 = x144 * x145; auto x147 = SS * x130; auto x148 = AXG * x147;
        auto x149 = 4.0 * MUs; auto x150 = x148 * x149; auto x151 = 5.0 * AXG;
        auto x152 = x131 * x151; auto x153 = x34 * x4; auto x154 = x8 * x92; auto x155 = x154 * x63;
        auto x156 = 3.0 * x155; auto x157 = kL * x156; auto x158 = x13 * x37; auto x159 = MGs * x26;
        auto x160 = 2.0 * x159; auto x161 = AXG * x133; auto x162 = x40 * x45;
        auto x163 = kL * x162; auto x164 = MGs * MXs; auto x165 = x18 * x47;
        auto x166 = x164 * x165; auto x167 = MXs * x128; auto x168 = 3.0 * x7;
        auto x169 = x40 * x65; auto x170 = AXG * x167; auto x171 = AXG * x169;
        auto x172 = x109 * x65; auto x173 = x169 * x7; auto x174 = pow(MG, 2);
        auto x175 = x174 * x26; auto x176 = s * x18; auto x177 = MXs * x176;
        auto x178 = x175 * x177; auto x179 = MG * x93; auto x180 = x14 * x179;
        auto x181 = kR * x180; auto x182 = x177 * x181; auto x183 = MG * iLGp;
        auto x184 = x183 * x2; auto x185 = MX * x184; auto x186 = kR * x185;
        auto x187 = x177 * x186; auto x188 = x133 * x180; auto x189 = kR * s; auto x190 = s * x147;
        auto x191 = s * x181; auto x192 = x19 * x43; auto x193 = 4.0 * x26; auto x194 = s * x83;
        auto x195 = x25 * x34; auto x196 = kR * x195; auto x197 = 4.0 * s; auto x198 = x147 * x26;
        auto x199 = 4.0 * x37; auto x200 = x199 * x96; auto x201 = x133 * x25;
        auto x202 = kR * x201; auto x203 = MX * x94; auto x204 = kR * x203;
        auto x205 = 2.0 * x204 * x8; auto x206 = x184 * x92; auto x207 = kR * x199 * x206;
        auto x208 = 4.0 * u; auto x209 = x208 * x37; auto x210 = x186 * x209;
        auto x211 = x203 * x41; auto x212 = AXG * x131; auto x213 = s * x41;
        auto x214 = x204 * x209; auto x215 = MGs * u; auto x216 = x215 * x31;
        auto x217 = pow(MX, 5); auto x218 = x217 * x47; auto x219 = x184 * x218;
        auto x220 = kR * x219; auto x221 = MGs * x47; auto x222 = kR * x85; auto x223 = x28 * x47;
        auto x224 = x159 * x223; auto x225 = Mqi * x184; auto x226 = u * x7;
        auto x227 = kL * x225 * x226; auto x228 = x47 * x84; auto x229 = x181 * x223;
        auto x230 = MGs * s * x130 * x193; auto x231 = u * x168; auto x232 = x231 * x65;
        auto x233 = s * x208; auto x234 = x206 * x233; auto x235 = kR * x234; auto x236 = x47 * x89;
        auto x237 = x130 * x197; auto x238 = x181 * x237; auto x239 = AXG * x232;
        auto x240 = MX * x1; auto x241 = kL * x240; auto x242 = x241 * x63; auto x243 = AXG * x242;
        auto x244 = x242 - x243; auto x245 = UU * syFC4; auto x246 = u * x1; auto x247 = x246 * x4;
        auto x248 = x28 * x33; auto x249 = x248 * x4; auto x250 = x241 * x57;
        auto x251 = x130 * x33; auto x252 = x57 * x92; auto x253 = kL * x33;
        auto x254 = x252 * x253; auto x255 = x251 * x4; auto x256 = kR * x33;
        auto x257 = x256 * x95; auto x258 = x63 * x92; auto x259 = x253 * x258;
        auto x260 = x206 * x256; auto x261 = AXG * x250; auto x262 = u * x253 * x64;
        auto x263 = u * x33; auto x264 = x263 * x59; auto x265 = UU * syFC9; auto x266 = 4.0 * x3;
        auto x267 = AXG * pow(u, 3); auto x268 = x267 * x37; auto x269 = MUs * x3;
        auto x270 = 3.0 * x269; auto x271 = SS * x267 * x47; auto x272 = x154 * x57;
        auto x273 = 3.0 * x272; auto x274 = x161 * x3; auto x275 = x101 * x112;
        auto x276 = x19 * x7; auto x277 = x109 * x58; auto x278 = x58 * x9; auto x279 = x192 * x44;
        auto x280 = 5.0 * x64; auto x281 = AXG * x19; auto x282 = x136 * x64; auto x283 = x50 * x58;
        auto x284 = kL * syFC12; auto x285 = x3 * x76; auto x286 = x3 * x81; auto x287 = x116 * x3;
        auto x288 = x13 * x3; auto x289 = x288 * x83; auto x290 = x118 * pow(MX, 8);
        auto x291 = x3 * x84; auto x292 = 6.0 * x291; auto x293 = x21 * x3;
        auto x294 = x118 * pow(MX, 7); auto x295 = UU * x116; auto x296 = x295 * x3;
        auto x297 = AXG * x296; auto x298 = 8.0 * MUs; auto x299 = AXG * x287;
        auto x300 = x258 * x91; auto x301 = x217 * x57; auto x302 = MXs * x18;
        auto x303 = AXG * x79; auto x304 = MXs * x3; auto x305 = MXs * x252;
        auto x306 = x118 * x301; auto x307 = UU * x89; auto x308 = kL * syFC13;
        auto x309 = 2.0 * x34; auto x310 = x10 * x21; auto x311 = 4.0 * x6; auto x312 = MUs * x309;
        auto x313 = 4.0 * x30; auto x314 = x217 * x63; auto x315 = AXG * x122;
        auto x316 = x147 * x92; auto x317 = x311 * x314; auto x318 = MUs * x35;
        auto x319 = x3 * x318; auto x320 = MUs * x122; auto x321 = x145 * x30; auto x322 = MXs * x0;
        auto x323 = AXG * x322; auto x324 = 12.0 * x291; auto x325 = MUs * x147;
        auto x326 = 4.0 * x161; auto x327 = 8.0 * x144; auto x328 = 8.0 * x323;
        auto x329 = x37 * x76; auto x330 = x3 * x329; auto x331 = 3.0 * x330;
        auto x332 = x165 * x288; auto x333 = MXs * x293; auto x334 = SS * x333;
        auto x335 = x334 * x47; auto x336 = x117 * x267; auto x337 = x267 * x78;
        auto x338 = x333 * x37; auto x339 = 3.0 * x338; auto x340 = 6.0 * x269;
        auto x341 = x3 * x89; auto x342 = x19 * x252; auto x343 = AXG * x158;
        auto x344 = x128 * x252; auto x345 = x252 * x37; auto x346 = MXs * x345;
        auto x347 = AXG * x344; auto x348 = AXG * x346; auto x349 = x252 * x281;
        auto x350 = -x336 * x64 + x337 * x64; auto x351 = MGs * x3; auto x352 = 4.0 * x34;
        auto x353 = MXs * x309; auto x354 = MGs * x10; auto x355 = x24 * x35; auto x356 = 4.0 * x35;
        auto x357 = x149 * x30; auto x358 = 2.0 * x318; auto x359 = kL * syFC7;
        auto x360 = AXG * x272; auto x361 = s * x195; auto x362 = x112 * x225;
        auto x363 = x18 * x45; auto x364 = x112 * x58; auto x365 = 5.0 * x136;
        auto x366 = x35 * x351; auto x367 = s * x34; auto x368 = x367 * x64; auto x369 = AXG * x363;
        auto x370 = x136 * x58; auto x371 = AXG * x277; auto x372 = s * x35; auto x373 = x372 * x58;
        auto x374 = x35 * x47; auto x375 = x225 * x374; auto x376 = x176 * x288;
        auto x377 = s * x334; auto x378 = MUs * x10; auto x379 = MXs * x10; auto x380 = x176 * x252;
        auto x381 = x3 * x7; auto x382 = MX * x8; auto x383 = x107 * x382; auto x384 = x199 * x258;
        auto x385 = x24 * x3; auto x386 = x100 * x382; auto x387 = 4.0 * x345;
        auto x388 = x133 * x3; auto x389 = s * x388; auto x390 = x209 * x58; auto x391 = x209 * x64;
        auto x392 = x21 * x25; auto x393 = x25 * x35; auto x394 = MXs * x25; auto x395 = x25 * x30;
        auto x396 = x184 * x217; auto x397 = x25 * x318; auto x398 = x311 * x396;
        auto x399 = MUs * x206; auto x400 = MUs * x25; auto x401 = x147 * x185;
        auto x402 = kR * syFC13; auto x403 = MGs * x25; auto x404 = x206 * x352;
        auto x405 = x206 * x313; auto x406 = x185 * x312; auto x407 = x185 * x353;
        auto x408 = x206 * x356; auto x409 = x185 * x357; auto x410 = x185 * x30;
        auto x411 = x24 * x410; auto x412 = x179 * x99; auto x413 = x185 * x358;
        auto x414 = kR * syFC7; auto x415 = x154 * x94; auto x416 = AXG * x415;
        auto x417 = x109 * x203; auto x418 = x203 * x9; auto x419 = x203 * x367;
        auto x420 = x18 * x203; auto x421 = x420 * x7; auto x422 = x112 * x203;
        auto x423 = x112 * x185; auto x424 = 2.0 * x423; auto x425 = 3.0 * x185;
        auto x426 = AXG * x421; auto x427 = x136 * x203; auto x428 = AXG * x417;
        auto x429 = x203 * x50; auto x430 = x203 * x372; auto x431 = x180 * x35;
        auto x432 = x136 * x185; auto x433 = x183 * x71; auto x434 = MX * x2;
        auto x435 = x433 * x434; auto x436 = x154 * x184; auto x437 = x203 * x40;
        auto x438 = x437 * x7; auto x439 = AXG * x433; auto x440 = x434 * x439;
        auto x441 = AXG * x436; auto x442 = x18 * x185; auto x443 = x168 * x442;
        auto x444 = x425 * x56; auto x445 = x149 * x8; auto x446 = x185 * x445;
        auto x447 = 4.0 * x9; auto x448 = x185 * x40; auto x449 = x203 * x61;
        auto x450 = 4.0 * x185; auto x451 = AXG * x448; auto x452 = AXG * x443;
        auto x453 = 7.0 * x61; auto x454 = x185 * x453; auto x455 = x25 * x329;
        auto x456 = x16 * x176; auto x457 = x23 * x394; auto x458 = MXs * x392;
        auto x459 = x37 * x458; auto x460 = 2.0 * x113; auto x461 = 2.0 * x400;
        auto x462 = x25 * x89; auto x463 = x176 * x95; auto x464 = x128 * x95;
        auto x465 = x105 * x37; auto x466 = x185 * x337; auto x467 = x185 * x336;
        auto x468 = x25 * x29; auto x469 = x143 * x25; auto x470 = UU * kR; auto x471 = MUs * x470;
        auto x472 = x394 * x471; auto x473 = UU * x4; auto x474 = MUs * MXs;
        auto x475 = x473 * x474; auto x476 = x30 * x4; auto x477 = x26 * x302;
        auto x478 = x302 * x4;
        auto x479 = -AXG * x472 - AXG * x475 + AXG * x477 + AXG * x478 + x133 * x4 - x143 * x4 +
                    x144 * x4 + x147 * x4 - x148 * x4 - x153 - x161 * x4 + x202 - x322 * x4 +
                    x323 * x4 + x36 + x472 + x475 - x476 - x477 - x478 + x67;
        auto x480 = x13 * x31; auto x481 = MUs * u * x31; auto x482 = x31 * x474;
        auto x483 = x222 * x47; auto x484 = x168 * x258; auto x485 = kL * x484;
        auto x486 = Mqi * x43 * x68 * x92; auto x487 = u * x45 * x68; auto x488 = UU * syFC10;
        auto x489 = x3 * x8; auto x490 = x13 * x489; auto x491 = -x112 * x64; auto x492 = x18 * x58;
        auto x493 = x18 * x64; auto x494 = AXG * x493; auto x495 = AXG * x492;
        auto x496 = x372 * x64; auto x497 = 3.0 * x64; auto x498 = x109 * x497;
        auto x499 = x125 * x64 + x498; auto x500 = x122 * x291; auto x501 = 4.0 * x148;
        auto x502 = x19 * x217; auto x503 = x258 * x79; auto x504 = x118 * x314;
        auto x505 = -MUs * x504 + MXs * x503 - MXs * x504 + x294 * x63 + x303 * x314 - x502 * x63;
        auto x506 = x174 * x3; auto x507 = x493 * x7; auto x508 = x131 * x506;
        auto x509 = x136 * x225; auto x510 = s * x506; auto x511 = kL * syFC2;
        auto x512 = -x168 * x493 + x168 * x494; auto x513 = kL * syFC6; auto x514 = x445 * x64;
        auto x515 = x168 * x40 * x64; auto x516 = x0 * x25; auto x517 = x71 * x93;
        auto x518 = AXG * x517; auto x519 = x185 * x41 * x7; auto x520 = kR * syFC12;
        auto x521 = -x424; auto x522 = x206 * x41; auto x523 = s * x522; auto x524 = 2.0 * x185;
        auto x525 = x133 * x185; auto x526 = x161 * x185; auto x527 = 6.0 * x185;
        auto x528 = x394 * x47; auto x529 = x203 * x30; auto x530 = x16 * x8;
        auto x531 = x33 * x420; auto x532 = x33 * x442; auto x533 = x103 * x203;
        auto x534 = x109 * x425; auto x535 = -AXG * x534 - x425 * x50 + x425 * x9 + x534;
        auto x536 = x116 * x25; auto x537 = AXG * x536; auto x538 = x25 * x295;
        auto x539 = AXG * x538; auto x540 = x217 * x94; auto x541 = x118 * x540;
        auto x542 = x307 * x400; auto x543 = MXs * x206; auto x544 = x118 * x396;
        auto x545 = -MUs * x544 - MXs * x544 + x184 * x294 - x19 * x396 + x206 * x91 + x303 * x396 +
                    x543 * x79;
        auto x546 = x185 * x35; auto x547 = kR * syFC6; auto x548 = x16 * x470;
        auto x549 = x13 * x473; auto x550 = syFC5 * x47; auto x551 = x470 * x95;
        auto x552 = x18 * x186; auto x553 = x203 * x471; auto x554 = MXs * x470;
        auto x555 = x203 * x554; auto x556 = x185 * x471; auto x557 = x185 * x554;
        auto x558 = kR * x420; auto x559 = syFC11 * x33; auto x560 = x15 * x42;
        auto x561 = x43 * x84; auto x562 = x184 * x240; auto x563 = MUs * x43; auto x564 = u * x563;
        auto x565 = x240 * x94; auto x566 = x168 * x95; auto x567 = x168 * x206;
        auto x568 = AXG * x565; auto x569 = AXG * x562; auto x570 = AXG * x567;
        auto x571 = x185 * x231; auto x572 = x203 * x231; auto x573 = AXG * x571;
        auto x574 = kR * x488; auto x575 = x47 * x536; auto x576 = MXs * x563;
        auto x577 = x16 * x233; auto x578 = x223 * x400; auto x579 = x218 * x94;
        auto x580 = s * x25; auto x581 = u * x120 * x580; auto x582 = x185 * x223;
        auto x583 = x130 * x149 * x580; auto x584 = x233 * x95; auto x585 = x203 * x223;
        auto x586 = AXG * x582; auto x587 = 2.0 * x296; auto x588 = UU * x217;
        auto x589 = x100 * x588; auto x590 = UU * x92; auto x591 = x100 * x590;
        auto x592 = x288 * x311; auto x593 = x588 * x63; auto x594 = AXG * x593;
        auto x595 = 6.0 * x594; auto x596 = x149 * x3; auto x597 = x133 * x58;
        auto x598 = x590 * x63; auto x599 = MUs * x598; auto x600 = AXG * x599;
        auto x601 = MXs * x598; auto x602 = AXG * x601; auto x603 = x326 * x58;
        auto x604 = 6.0 * x58; auto x605 = x19 * x258; auto x606 = AXG * x605;
        auto x607 = AXG * x123; auto x608 = x124 * x3;
        auto x609 = AXG * x608 + x3 * x607 - x605 + x606; auto x610 = s * x284; auto x611 = -x593;
        auto x612 = x117 * x291; auto x613 = MUs * x612; auto x614 = MXs * x612;
        auto x615 = x0 * x293; auto x616 = 3.0 * UU; auto x617 = x258 * x616; auto x618 = x143 * x3;
        auto x619 = MXs * x618; auto x620 = x140 * x3; auto x621 = x141 * x3; auto x622 = x293 * x6;
        auto x623 = x57 * x588; auto x624 = x144 * x3; auto x625 = MXs * x624;
        auto x626 = -x296 + x297 + x623; auto x627 = s * x308; auto x628 = x302 * x506;
        auto x629 = MXs * x269; auto x630 = UU * x174 * x629; auto x631 = UU * x64;
        auto x632 = x474 * x631; auto x633 = x133 * x64; auto x634 = x133 * x225;
        auto x635 = x147 * x225; auto x636 = x302 * x64; auto x637 = x143 * x64;
        auto x638 = x143 * x225; auto x639 = x225 * x322; auto x640 = x144 * x64;
        auto x641 = x144 * x225; auto x642 = x225 * x323; auto x643 = x161 * x64;
        auto x644 = x161 * x225; auto x645 = x148 * x225; auto x646 = UU * x225;
        auto x647 = x474 * x646; auto x648 = x225 * x302;
        auto x649 = -AXG * x647 + AXG * x648 + x647 - x648; auto x650 = 8.0 * AXG;
        auto x651 = -x599 * x650; auto x652 = MGs * UU; auto x653 = x288 * x652;
        auto x654 = 3.0 * x653; auto x655 = AXG * x623; auto x656 = x164 * x3;
        auto x657 = x656 * x79; auto x658 = x117 * x258; auto x659 = x13 * x225;
        auto x660 = x616 * x659; auto x661 = UU * x10; auto x662 = MXs * x661;
        auto x663 = x629 * x652; auto x664 = 3.0 * x663; auto x665 = x18 * x258;
        auto x666 = x225 * x474 * x616; auto x667 = MXs * x79; auto x668 = MXs * x303;
        auto x669 = x291 * x652 - x341 * x652; auto x670 = s * x359; auto x671 = x120 * x3;
        auto x672 = SS * x671; auto x673 = MGs * x78; auto x674 = x217 * x83; auto x675 = AXG * x78;
        auto x676 = SS * x120; auto x677 = x64 * x79; auto x678 = AXG * x145;
        auto x679 = x652 * x678; auto x680 = AXG * x598; auto x681 = x145 * x680;
        auto x682 = UU * x650; auto x683 = AXG * x18 * x24; auto x684 = 2.0 * x538;
        auto x685 = 2.0 * UU; auto x686 = x685 * x77; auto x687 = x588 * x94;
        auto x688 = x117 * x85; auto x689 = MUs * x688; auto x690 = MXs * x688;
        auto x691 = MUs * x16; auto x692 = MUs * x685; auto x693 = x19 * x95;
        auto x694 = x114 * x25; auto x695 = AXG * x694; auto x696 = x140 * x25;
        auto x697 = x141 * x25; auto x698 = AXG * x687; auto x699 = 6.0 * MUs;
        auto x700 = x590 * x94; auto x701 = AXG * x700; auto x702 = 6.0 * MXs;
        auto x703 = x132 * x25; auto x704 = AXG * x693; auto x705 = x703 - x704;
        auto x706 = x19 * x206; auto x707 = AXG * x706; auto x708 = -x706 + x707;
        auto x709 = s * x520; auto x710 = x21 * x516; auto x711 = MXs * x469;
        auto x712 = x16 * x311; auto x713 = x138 * x25; auto x714 = x392 * x6;
        auto x715 = x41 * x95; auto x716 = x144 * x394; auto x717 = x143 * x203;
        auto x718 = x203 * x322; auto x719 = AXG * x715; auto x720 = x133 * x203;
        auto x721 = x184 * x590; auto x722 = MUs * x721; auto x723 = AXG * x722;
        auto x724 = MXs * x721; auto x725 = AXG * x724; auto x726 = x144 * x203;
        auto x727 = x203 * x323; auto x728 = x184 * x588; auto x729 = AXG * x728;
        auto x730 = -2.0 * x728 + 6.0 * x729; auto x731 = x124 * x25;
        auto x732 = AXG * x731 + x203 * x326 + x25 * x607; auto x733 = -x728;
        auto x734 = x206 * x616; auto x735 = x147 * x400; auto x736 = MUs * x203;
        auto x737 = s * x402; auto x738 = x174 * x25; auto x739 = x143 * x185;
        auto x740 = x185 * x322; auto x741 = x180 * x323; auto x742 = x144 * x185;
        auto x743 = x185 * x323;
        auto x744 = -x180 * x34 + x431 + x525 - x526 - x739 - x740 + x741 + x742 + x743;
        auto x745 = kR * syFC2; auto x746 = -x650 * x722; auto x747 = x16 * x652;
        auto x748 = 3.0 * x747; auto x749 = x164 * x25 * x79; auto x750 = x15 * x179;
        auto x751 = x616 * x750; auto x752 = x394 * x89; auto x753 = 3.0 * MUs * x394 * x652;
        auto x754 = x18 * x206; auto x755 = x180 * x474; auto x756 = x616 * x755;
        auto x757 = -x462 * x652 + x652 * x85; auto x758 = s * x414; auto x759 = AXG * x721;
        auto x760 = -x145 * x759; auto x761 = x21 * x40; auto x762 = x25 * x761;
        auto x763 = x16 * x40; auto x764 = x40 * x95; auto x765 = MUs * x315; auto x766 = UU * x84;
        auto x767 = x180 * x307 - x180 * x766; auto x768 = 2.0 * x25; auto x769 = 3.0 * x469;
        auto x770 = x616 * x95; auto x771 = x185 * x19; auto x772 = 5.0 * UU;
        auto x773 = MXs * x772; auto x774 = x144 * x25; auto x775 = x0 * x203;
        auto x776 = MUs * x185; auto x777 = UU * x151; auto x778 = MXs * x777;
        auto x779 = x79 * x85; auto x780 = UU * x678; auto x781 = UU * x691; auto x782 = x0 * x185;
        auto x783 = 3.0 * x448; auto x784 = x185 * x6;
        auto x785 = -AXG * x437 + x161 * x25 + x195 + x203 * x6 - x393 + x395 + x437 - x468 + x469 -
                    x774 - x775;
        auto x786 = x164 * x381; auto x787 = x223 * x351; auto x788 = x218 * x63;
        auto x789 = MXs * x7; auto x790 = x225 * x789; auto x791 = x233 * x258;
        auto x792 = x223 * x64; auto x793 = x223 * x225; auto x794 = x237 * x351;
        auto x795 = x225 * x237; auto x796 = x215 * x43; auto x797 = x164 * x43;
        auto x798 = x180 * x789; auto x799 = x180 * x226; auto x800 = -x562 + x569;
        auto x801 = x246 * x25; auto x802 = x248 * x25; auto x803 = AXG * x25;
        auto x804 = x203 * x263; auto x805 = x185 * x263; auto x806 = x14 * x174;
        auto x807 = MG * x434; auto x808 = x109 * x807; auto x809 = x7 * x807;
        auto x810 = x40 * x809; auto x811 = x18 * x809; auto x812 = x661 * x76;
        auto x813 = x21 * x662; auto x814 = x13 * x269; auto x815 = AXG * x3;
        auto x816 = x57 * x590; auto x817 = AXG * x816; auto x818 = x120 * x6;
        auto x819 = x3 * x818; auto x820 = x323 * x64; auto x821 = x252 * x41;
        auto x822 = x53 * x92;
        auto x823 = -AXG * x821 + AXG * x822 + x138 * x3 + x326 * x64 + x821 - x822;
        auto x824 = x269 * x30; auto x825 = x322 * x64; auto x826 = x143 * x58;
        auto x827 = x322 * x58; auto x828 = x144 * x58; auto x829 = x323 * x58;
        auto x830 = -x258 * x311 + x269 * x29; auto x831 = x143 * x351; auto x832 = x322 * x351;
        auto x833 = x258 * x40; auto x834 = x144 * x351; auto x835 = x323 * x351;
        auto x836 = x3 * x761; auto x837 = x288 * x40; auto x838 = x252 * x40;
        auto x839 = x225 * x766; auto x840 = x646 * x89; auto x841 = x252 * x83;
        auto x842 = x258 * x83; auto x843 = x58 * x79; auto x844 = x118 * x58;
        auto x845 = x118 * x64; auto x846 = x616 * x77; auto x847 = x458 * x616;
        auto x848 = 6.0 * x781; auto x849 = 7.0 * x701; auto x850 = MUs * x395;
        auto x851 = MUs * x468 - x206 * x311; auto x852 = x147 * x403; auto x853 = x322 * x403;
        auto x854 = x206 * x40; auto x855 = x144 * x403; auto x856 = x143 * x180;
        auto x857 = x180 * x322; auto x858 = x323 * x403; auto x859 = x148 * x403;
        auto x860 = x144 * x180; auto x861 = x13 * x646; auto x862 = x18 * x656;
        auto x863 = x47 * x513; auto x864 = UU * x750; auto x865 = UU * x755;
        auto x866 = x180 * x302; auto x867 = x47 * x547; auto x868 = MUs * x631;
        auto x869 = MXs * x631; auto x870 = UU * x58; auto x871 = MUs * x870;
        auto x872 = MXs * x870; auto x873 = x40 * x58; auto x874 = x10 * x116;
        auto x875 = x100 * x217; auto x876 = x208 * x288; auto x877 = x107 * x217;
        auto x878 = x28 * x378; auto x879 = u * x671; auto x880 = MX * x28; auto x881 = x208 * x252;
        auto x882 = x130 * x596; auto x883 = x107 * x880; auto x884 = x208 * x258;
        auto x885 = x0 * x99; auto x886 = x685 * x807; auto x887 = 2.0 * x56;
        auto x888 = x103 * x183;
        ,
        TR * TR * gs * gs * (Nc - 1) * (Nc + 1) *
            (-iLGp * x745 *
                 (-AXG * x808 - AXG * x810 + AXG * x811 + x112 * x806 - x136 * x806 - x367 * x806 +
                  x372 * x806 - x56 * x807 + x61 * x807 + x808 + x810 - x811) +
             kL * s * x488 *
                 (-AXG * x874 - AXG * x875 - AXG * x876 + AXG * x877 + AXG * x878 + AXG * x879 -
                  AXG * x881 - AXG * x882 + AXG * x883 + AXG * x884 + x100 * x880 - x378 * x84 +
                  x378 * x89 - x379 * x84 + x379 * x89 + x874 + x875 + x876 - x877 - x878 - x879 +
                  x881 + x882 - x883 - x884) +
             kL * syFC8 *
                 (-AXG * x383 - AXG * x384 + AXG * x386 + AXG * x387 - AXG * x390 + AXG * x391 +
                  s * x54 + x113 * x266 + x129 * x266 - x136 * x266 - x158 * x266 - x167 * x266 +
                  x176 * x385 - x19 * x381 - x190 * x266 + x194 * x293 + x20 * x58 - x20 * x64 +
                  x212 * x266 - x213 * x58 + x266 * x367 + x383 + x384 - x386 - x387 - 8.0 * x389 +
                  x390 - x391) -
             kL * x245 *
                 (-AXG * x484 + AXG * x786 - AXG * x787 + AXG * x788 - AXG * x790 - AXG * x791 +
                  AXG * x792 + AXG * x793 + AXG * x794 - AXG * x795 + x221 * x291 - x225 * x228 +
                  x225 * x236 - x236 * x351 + x484 - x786 + x787 - x788 + x790 + x791 - x792 -
                  x793 - x794 + x795) -
             kL * x559 *
                 (-AXG * x868 - AXG * x869 + AXG * x871 + AXG * x872 + AXG * x873 - x274 + x388 +
                  x44 * x516 + x492 - x493 + x494 - x495 - x598 - x618 + x624 + x680 + x816 - x817 +
                  x868 + x869 - x871 - x872 - x873) -
             kR * syFC8 * x42 * (AXG * x885 - AXG * x886 - x51 + x885 + x886) -
             kR * x245 *
                 (-AXG * x796 + AXG * x797 - AXG * x798 + AXG * x799 + x567 - x570 - x571 + x573 +
                  x586 + x796 - x797 + x798 - x799 + x800) -
             kR * x265 *
                 (-AXG * x801 - AXG * x802 - AXG * x804 + AXG * x805 + x251 * x803 + x565 - x568 +
                  x800 + x801 + x802 + x804 - x805) +
             kR * x559 * (-x201 - x448 + x451 + x721 - x759 + x782 - x784 + x785) +
             s * syFC1 *
                 (kR * x393 - kR * x395 + kR * x468 - kR * x469 + x144 * x26 - x148 * x26 -
                  x161 * x26 - x196 + x198 - x26 * x322 + x26 * x323 + x479) +
             s * x511 *
                 (-AXG * x628 + AXG * x630 - AXG * x632 + AXG * x636 + x147 * x64 - x148 * x64 +
                  x628 - x630 + x632 + x633 + x634 + x635 - x636 - x637 - x638 - x639 + x640 +
                  x641 + x642 - x643 - x644 - x645 + x649) +
             s * x745 *
                 (-x144 * x738 - x147 * x738 - x148 * x185 + x148 * x738 + x161 * x738 -
                  x174 * x201 + x174 * x469 + x180 * x29 - x180 * x30 + x185 * x29 - x185 * x34 +
                  x322 * x738 - x323 * x738 + x401 - x410 + x546 + x744) +
             syFC12 * (AXG * x39 - AXG * x46 + AXG * x55 - x0 * x5 - x12 * x50 + x12 * x9 +
                       x17 * x20 - x22 * x48 * x49 + x23 * x27 + 5.0 * x32 - x33 * x36 - x39 + x46 +
                       x5 * x6 - x55 - 3.0 * x60 + 7.0 * x62 - 7.0 * x66 + x70 + x75) +
             syFC13 * (-AXG * x104 + AXG * x106 + AXG * x109 * x69 - kR * x22 * x86 + x101 * x102 -
                       x102 * x108 + x104 + x105 * x80 - x106 + x110 + x111 * x59 - x111 * x65 +
                       x17 * x21 * x83 + x26 * x81 * x82 + x38 * x90 - x69 * x9 + x75 + x77 * x80 -
                       x88 + x91 * x96 + x98) -
             syFC2 *
                 (-AXG * x172 - AXG * x173 + AXG * x178 - AXG * x182 - AXG * x187 + x129 * x181 +
                  x131 * x181 - x142 * x181 + x148 * x191 + x161 * x191 + x167 * x175 -
                  x167 * x181 - x167 * x186 - x170 * x175 + x170 * x181 + x170 * x186 + x172 +
                  x173 - x178 - x181 * x190 + x182 + x187 - x188 * x189 - x56 * x65 + x66) +
             syFC3 * x26 *
                 (x112 - x113 - x125 + x127 + x129 + x131 + x134 + x137 - x142 - x152 - x158 -
                  x167 + x170 + x312 - x321 + x336 - x337 + x343 - x353 + x355 - x357 + x358 - x50 -
                  x56 + x61 - x818 + x9) -
             syFC3 * x4 *
                 (-x114 - x115 - x119 - x121 + x123 + x124 + x126 - x13 * x135 - x132 + x138 +
                  x140 + x141 - x146 + x150 + x158 + x167 - x170 - x312 + x321 - x336 + x337 -
                  x343 + x353 - x355 + x357 - x358 + x818) +
             syFC3 * (x112 * x4 - x113 * x4 + x114 * x26 + x115 * x26 + x119 * x26 + x121 * x26 -
                      x123 * x26 - x124 * x26 - x125 * x4 - x126 * x26 + x127 * x4 + x129 * x4 +
                      x131 * x4 + x132 * x26 + x134 * x4 + x135 * x17 + x137 * x4 - x138 * x26 -
                      x140 * x26 - x141 * x26 - x142 * x4 + x146 * x26 - x150 * x26 - x152 * x4 -
                      x4 * x50 + x4 * x9 - x87 + x97) +
             syFC6 * (-AXG * x157 - AXG * x163 + x110 - x133 * x31 + x153 * x7 + x157 -
                      x158 * x160 - x160 * x167 + x160 * x170 + x161 * x31 + x163 + x166 * x26 -
                      x166 * x49 - x168 * x169 + x168 * x171 - x31 * x35 + x32 - x60 - x61 * x69 +
                      x62 + x70 + x88 - x98) -
             syFC8 * (-AXG * x200 + AXG * x205 + AXG * x207 - AXG * x210 + AXG * x214 + kR * x192 +
                      8.0 * s * x202 + x11 * x56 + x11 * x61 - x113 * x193 - x129 * x193 +
                      x136 * x193 + x158 * x193 + x167 * x193 - x176 * x27 + x186 * x20 -
                      x186 * x213 + x189 * x211 - x193 * x212 - x194 * x38 - x196 * x197 +
                      x197 * x198 - x20 * x204 + x200 - x205 - x207 + x210 - x214 - x31 * x41) +
             u * x359 *
                 (MGs * x672 - MUs * x646 * x678 + MXs * x677 + x139 * x225 - x139 * x351 -
                  x225 * x676 + x225 * x683 + x225 * x90 + x225 * x91 + x258 * x82 + x269 * x679 -
                  x288 * x673 + x314 * x675 - x351 * x90 - x503 + x606 - x63 * x674 + x64 * x91 +
                  x650 * x653 + x651 - x659 * x682 + x659 * x78 - x681) +
             u * x402 *
                 (16.0 * AXG * x781 + MUs * x650 * x700 + MXs * SS * x86 - 12.0 * UU * x752 +
                  x119 * x25 - x16 * x281 - x203 * x667 - x203 * x91 + x392 * x780 - x462 * x82 -
                  x536 * x83 + x537 * x78 - x540 * x675 - x541 - 12.0 * x542 + x674 * x94 +
                  x678 * x700 + x682 * x77 + x705 - x77 * x78 + x779 + x79 * x95 - x82 * x95) +
             u * x414 *
                 (-MUs * x180 * x780 + x139 * x180 - x139 * x403 - x16 * x673 - x180 * x676 +
                  x180 * x683 + x180 * x90 + x180 * x91 - x184 * x674 + x185 * x667 + x185 * x91 -
                  x206 * x79 + x206 * x82 + x396 * x675 + x400 * x679 + x403 * x676 - x403 * x90 +
                  x650 * x747 - x682 * x750 + x707 + x746 + x750 * x78 + x760) +
             x245 * (AXG * x216 - AXG * x220 + AXG * x224 - AXG * x227 - AXG * x229 - AXG * x230 +
                     AXG * x235 + AXG * x238 + x181 * x228 - x181 * x236 + x186 * x223 - x216 +
                     x220 - x221 * x222 + x221 * x26 * x89 - x224 + x227 + x229 + x230 + x232 -
                     x235 - x238 - x239 + x244) -
             x26 * x550 * (x143 - x144 - x147 + x148 + x161 - x29 + x30 + x322 - x323 + x34 - x35) +
             x265 * (AXG * x247 + AXG * x249 - AXG * x254 - AXG * x255 - AXG * x257 + AXG * x259 +
                     AXG * x260 - AXG * x262 + AXG * x264 + x244 - x247 - x249 - x250 + x251 * x26 +
                     x254 + x255 + x257 - x259 - x260 + x261 + x262 - x264) -
             x28 * x308 *
                 (-AXG * x289 + AXG * x672 + AXG * x677 - AXG * x841 + AXG * x842 - AXG * x843 +
                  MUs * x844 - MUs * x845 + MXs * x844 - MXs * x845 + x118 * x252 - x118 * x258 +
                  x118 * x288 - x139 * x3 + x289 - x58 * x82 + x64 * x82 - x672 + x677 + x841 -
                  x842 - x843) +
             x284 * (-AXG * x156 + AXG * x273 + AXG * x279 + x109 * x151 * x64 - x109 * x280 -
                     x151 * x277 + x156 + x266 * x268 + x270 * x56 - x270 * x61 - x271 * x3 - x273 +
                     x274 * x33 + x275 + x276 * x64 + 5.0 * x277 + 5.0 * x278 - x279 + x280 * x50 -
                     x280 * x9 - x281 * x64 * x7 - 6.0 * x282 - 5.0 * x283) -
             x308 * (-MUs * x500 + x107 * x316 + x108 * x325 + x122 * x285 + x122 * x286 +
                     x123 * x58 - x135 * x285 - x138 * x269 - 16.0 * x144 * x288 + x144 * x324 -
                     x147 * x292 - x252 * x320 + x252 * x326 - x252 * x327 - x252 * x328 +
                     x288 * x326 - x293 * x328 + x293 * x501 + x505) -
             x308 * (AXG * x490 - x125 * x58 + x269 * x50 - x269 * x9 - 3.0 * x277 - 3.0 * x278 +
                     9.0 * x282 + 3.0 * x283 + x33 * x492 - x33 * x493 + x33 * x494 - x33 * x495 +
                     x364 - 9.0 * x370 + 3.0 * x371 + 4.0 * x373 + x489 * x84 - x489 * x89 - x490 +
                     x491 - 4.0 * x496 + x499) +
             x308 * (AXG * x331 - AXG * x332 - AXG * x335 + AXG * x339 + s * x342 - s * x349 -
                     x158 * x340 - x165 * x291 + x165 * x341 + x269 * x336 - x269 * x337 - x331 +
                     x332 + x335 + x336 * x58 - x337 * x58 - x339 + x340 * x343 - 3.0 * x344 -
                     3.0 * x346 + 7.0 * x347 + 7.0 * x348 + x350) +
             x308 * (-MUs * x306 - MXs * x306 + x132 * x304 - x139 * x288 - x139 * x293 -
                     x145 * x269 * x307 + x145 * x297 + x19 * x287 - x19 * x301 + x21 * x289 -
                     x22 * x292 + x252 * x91 + x285 * x79 + x286 * x82 - x290 * x3 - x292 * x302 +
                     x293 * x90 + x294 * x57 + x297 * x298 - x299 * x79 - x300 + x301 * x303 +
                     x305 * x79) +
             x308 * (-x10 * x148 * x84 + x100 * x316 + x101 * x325 + x122 * x299 + x123 * x64 -
                     x217 * x53 + x24 * x319 - x258 * x320 + x258 * x326 - x258 * x327 -
                     x258 * x328 - x269 * x321 + x287 * x311 - x287 * x41 + x293 * x309 -
                     x293 * x313 - x301 * x311 - x301 * x315 + x301 * x41 - x304 * x312 +
                     x310 * x35 + x314 * x315 + x317 - x323 * x324) -
             x359 * (x114 * x225 - x114 * x351 + x115 * x225 - x115 * x351 + x119 * x225 -
                     x119 * x351 + x123 * x351 - x126 * x225 + x126 * x351 + x132 * x225 -
                     x132 * x351 - x140 * x225 + x140 * x351 - x141 * x225 + x141 * x351 +
                     x150 * x351 + x300 + x505) -
             x359 * (AXG * x155 + AXG * x162 - AXG * x514 - AXG * x515 + AXG * x72 - x155 - x162 -
                     x225 * x50 - x225 * x56 + x225 * x61 + x225 * x9 + x351 * x56 - x351 * x61 +
                     x447 * x64 + x453 * x64 - x497 * x56 - 4.0 * x50 * x64 + x512 + x514 + x515 +
                     x56 * x58 - x58 * x61 - x72) +
             x359 * (-x101 * x318 - x225 * x312 + x225 * x321 + x225 * x353 - x225 * x355 +
                     x225 * x357 - x225 * x358 + x24 * x30 * x64 + x258 * x313 + x258 * x352 -
                     x258 * x356 + x312 * x351 - x312 * x64 + x317 + x318 * x354 - x321 * x351 +
                     x336 * x351 - x337 * x351 - x351 * x353 + x351 * x355 - x351 * x357 -
                     x353 * x64 + x357 * x64) +
             x359 * (AXG * x330 + AXG * x338 - AXG * x376 - AXG * x377 - AXG * x380 + x113 * x378 +
                     x113 * x379 - x127 * x378 - x158 * x378 - x176 * x291 + x176 * x341 -
                     x225 * x336 + x225 * x337 - x330 - x338 + x343 * x378 - x344 - x346 + x347 +
                     x348 + x350 + x376 + x377 + x380) +
             x359 *
                 (-x112 * x351 - x225 * x365 - x272 + x275 + x277 + x278 + x280 * x372 -
                  10.0 * x282 - x283 + x351 * x365 - x351 * x50 + x351 * x9 + x360 + x361 * x44 +
                  x362 - x363 - x364 - x366 * x47 - 3.0 * x368 + x369 + x370 - x371 - x373 + x375) +
             x402 * x434 *
                 (AXG * x888 - x111 * x183 + x111 * x93 + x183 * x887 + x433 - x439 - x517 + x518 -
                  x887 * x93 - x888) -
             x402 * (MUs * x541 + MXs * x541 + x122 * x25 * x81 - x132 * x394 + x139 * x16 +
                     x139 * x392 - x145 * x539 + x145 * x542 - x19 * x536 + x25 * x290 -
                     x294 * x94 - x298 * x539 + x302 * x86 - x303 * x540 + x392 * x501 +
                     x502 * x94 + x537 * x79 + x545) -
             x402 * (AXG * x530 - AXG * x531 + AXG * x532 + AXG * x533 + x400 * x50 + x400 * x56 -
                     x400 * x61 - x400 * x9 - 3.0 * x417 - 3.0 * x418 - x423 + 3.0 * x428 +
                     3.0 * x429 + 9.0 * x432 - x462 * x8 - x530 + x531 - x532 - x533 + x535 +
                     x8 * x85) -
             x402 * (-x120 * x195 + x120 * x393 - x120 * x395 - x203 * x336 + x203 * x337 -
                     x22 * x528 + x24 * x529 + x313 * x95 - x336 * x400 + x337 * x400 + x352 * x95 -
                     x356 * x95 - x404 - x405 + x406 + x407 + x408 - x409 - x411 + x413 - x466 +
                     x467) +
             x402 * (2.0 * MUs * x401 - x122 * x399 + x123 * x185 + 2.0 * x147 * x206 - x16 * x313 -
                     x16 * x352 + x16 * x356 + x203 * x312 + x203 * x353 - x203 * x357 +
                     x203 * x358 + x206 * x326 - x206 * x327 - x206 * x328 + 2.0 * x21 * x393 -
                     4.0 * x21 * x395 + x24 * x397 + x309 * x392 - x312 * x394 + x315 * x396 -
                     x321 * x400 - x396 * x41 + x398) -
             x414 * (x114 * x180 - x114 * x403 + x115 * x180 - x115 * x403 + x119 * x180 -
                     x119 * x403 + x123 * x403 - x126 * x180 + x126 * x403 + x132 * x180 -
                     x132 * x403 - x140 * x180 + x140 * x403 - x141 * x180 + x141 * x403 +
                     x150 * x403 + x545) +
             x414 * (-x180 * x312 + x180 * x321 + x180 * x353 - x180 * x355 + x180 * x357 +
                     x312 * x403 - x318 * x412 - x321 * x403 + x336 * x403 - x337 * x403 -
                     x353 * x403 + x355 * x403 - x357 * x403 + x358 * x403 + x398 + x404 + x405 -
                     x406 - x407 - x408 + x409 + x411 - x413) +
             x414 * (-AXG * x438 + AXG * x446 - x168 * x448 + x168 * x451 + x180 * x50 +
                     x180 * x56 - x180 * x61 - x180 * x9 - x185 * x447 - x203 * x56 - x403 * x56 +
                     x403 * x61 + x435 + x436 + x438 - x440 - x441 + x443 + x444 - x446 + x449 +
                     x450 * x50 - x452 - x454) +
             x414 * (x112 * x180 - x112 * x403 - x180 * x365 + 5.0 * x185 * x372 - x221 * x393 +
                     x365 * x403 - x367 * x425 - x403 * x50 + x403 * x9 - x415 + x416 + x417 +
                     x418 + x419 - x421 - x422 + x424 + x426 + x427 - x428 - x429 - x430 +
                     x431 * x47 - 10.0 * x432) +
             x414 * (AXG * x455 - AXG * x456 - AXG * x457 + AXG * x459 - AXG * x463 + AXG * x464 +
                     AXG * x465 - x127 * x461 - x158 * x461 + x176 * x462 - x176 * x85 -
                     x180 * x336 + x180 * x337 + x343 * x461 + x394 * x460 + x400 * x460 - x455 +
                     x456 + x457 - x459 + x463 - x464 - x465 + x466 - x467) -
             x488 * (AXG * x223 * x65 + AXG * x480 - AXG * x481 + AXG * x482 - AXG * x485 +
                     AXG * x486 - AXG * x487 - MUs * x236 * x26 + MUs * x483 + MXs * x483 - x232 +
                     x239 - x242 + x243 + x250 - x261 + x31 * x84 - x31 * x89 - x480 + x481 - x482 +
                     x485 - x486 + x487) -
             x511 * (AXG * x507 + AXG * x508 + x112 * x506 - x129 * x506 + x131 * x64 -
                     x136 * x506 + x142 * x506 - x148 * x510 - x161 * x510 + x174 * x389 +
                     x190 * x506 - x212 * x64 + x225 * x367 - x225 * x372 + x282 - x362 -
                     x367 * x506 + x368 + x372 * x506 + x491 - x496 - x507 - x508 + x509) -
             x513 * (-AXG * x498 + x101 * x136 + x112 * x354 + x125 * x225 - x136 * x354 + x272 -
                     x275 - x277 - x278 + x283 - x360 - 2.0 * x362 + x363 - x369 + x371 -
                     x374 * x64 - x375 + x497 * x9 + x499 + 2.0 * x509 + x512) +
             x520 * x7 *
                 (-AXG * x211 - AXG * x734 + AXG * x770 - AXG * x771 - x185 * x773 + x185 * x778 +
                  x203 * x281 + x203 * x773 - x203 * x778 + x211 + x315 * x400 + x322 * x768 -
                  x323 * x768 + x734 + x736 * x772 - x736 * x777 + x769 - x770 + x771 -
                  x772 * x776 - 3.0 * x774 - 3.0 * x775 + x776 * x777) -
             x520 * (-AXG * x519 - x1 * x25 * x6 + x1 * x516 + x29 * x43 - 5.0 * x30 * x43 +
                     x33 * x393 + x434 * x517 - x434 * x518 - x435 + x440 - x444 - 7.0 * x449 +
                     x454 + x519) -
             x520 * (AXG * x523 + x112 * x461 + x127 * x208 * x25 + x129 * x524 + x131 * x524 +
                     x136 * x24 * x25 + x136 * x461 - x142 * x527 - x149 * x361 - x197 * x525 +
                     x197 * x526 + x203 * x276 - x212 * x527 - 4.0 * x25 * x268 + x25 * x271 +
                     x367 * x450 - x372 * x450 - x393 * x48 - 4.0 * x419 + 2.0 * x422 - 6.0 * x427 +
                     4.0 * x430 + 6.0 * x432 + x521 - x523) +
             x547 * x7 * (AXG * x783 + 3.0 * x782 - x783 - 3.0 * x784 + x785) -
             x547 * (-x112 * x412 + x125 * x185 + x133 * x43 + x136 * x412 + x415 - x416 - x417 -
                     x418 + x421 - x426 + x428 + x429 + 2.0 * x432 - 3.0 * x436 + 3.0 * x441 -
                     x443 + x452 - x47 * x546 + x521 + x535) +
             x550 * (-AXG * x548 - AXG * x549 + x462 * x470 - x470 * x85 - x473 * x84 + x473 * x89 +
                     x479 + x548 + x549) +
             x559 * (AXG * x551 - AXG * x552 - AXG * x553 - AXG * x555 + AXG * x556 + AXG * x557 +
                     AXG * x558 + x0 * x65 + x153 - x169 + x171 - x36 + x476 - x551 + x552 + x553 +
                     x555 - x556 - x557 - x558 + x59 * x6 - x6 * x65 - x67) +
             x574 * (-AXG * x560 + AXG * x561 + AXG * x564 - AXG * x566 + AXG * x572 + x560 - x561 +
                     x562 - x564 - x565 + x566 - x567 + x568 - x569 + x570 + x571 - x572 - x573) +
             x574 * (-AXG * x219 + AXG * x234 - AXG * x575 - AXG * x576 - AXG * x577 + AXG * x578 +
                     AXG * x579 + AXG * x581 - AXG * x583 - AXG * x584 + AXG * x585 + x219 - x234 +
                     x528 * x89 + x575 + x576 + x577 - x578 - x579 - x581 + x582 + x583 + x584 -
                     x585 - x586) -
             x610 * (x101 * x143 + x101 * x322 - x121 * x3 - 4.0 * x195 * x44 + x29 * x378 +
                     x30 * x378 + x30 * x385 - x30 * x604 + x323 * x604 - x34 * x596 - x35 * x379 +
                     x352 * x64 + x356 * x58 - x356 * x64 + 2.0 * x44 * x468 - 4.0 * x633 -
                     6.0 * x640 + x819 - 6.0 * x820 + x823) +
             x610 * (AXG * x587 + MUs * x591 + MXs * x591 + x0 * x310 + x108 * x143 + x108 * x322 +
                     x143 * x379 + x144 * x379 - x144 * x604 - x147 * x596 - x288 * x41 -
                     x293 * x41 - x310 * x6 - x587 - x589 + x592 + x595 - 4.0 * x597 - 6.0 * x600 -
                     6.0 * x602 + x603 + x609) -
             x610 * (-AXG * x812 - AXG * x813 - SS * x24 * x293 + x114 * x815 + x115 * x815 +
                     x117 * x814 - x118 * x814 + x126 * x3 - x132 * x3 - x19 * x288 + x252 * x692 +
                     x305 * x685 - x342 + x349 - x613 - x614 + x620 + x621 - 2.0 * x623 +
                     6.0 * x655 - x699 * x817 - x702 * x817 + x812 + x813) -
             x627 * (x10 * x318 + x252 * x311 + x320 * x58 - x320 * x64 - x500 + x592 - x603 +
                     3.0 * x637 - 7.0 * x640 - x819 - 7.0 * x820 + x823 - 5.0 * x824 + 3.0 * x825 -
                     3.0 * x826 - 3.0 * x827 + 7.0 * x828 + 7.0 * x829 + x830) +
             x627 * (MUs * x617 + MXs * x617 - x147 * x378 + x148 * x378 + x151 * x593 -
                     x151 * x623 - 7.0 * x600 - 7.0 * x602 - x608 + x609 + x611 + x613 + x614 +
                     3.0 * x615 + 3.0 * x619 - x620 - x621 - 3.0 * x622 - 7.0 * x625 + x626) +
             x670 * (AXG * x654 - AXG * x657 - AXG * x660 + AXG * x664 - AXG * x666 + MUs * x658 +
                     MXs * x658 + x151 * x665 - x225 * x667 + x225 * x668 + x626 + x651 - x654 -
                     x655 + x657 + x660 - x662 * x89 - x664 - 5.0 * x665 + x666 + x669) -
             x670 * (x147 * x354 - x148 * x354 + x151 * x833 + x161 * x280 + x225 * x765 -
                     x269 * x34 + x319 - x328 * x64 - 3.0 * x633 - 2.0 * x635 + 4.0 * x637 +
                     3.0 * x638 + 3.0 * x639 - 8.0 * x640 - 3.0 * x641 - 7.0 * x642 + 2.0 * x645 -
                     x824 + 4.0 * x825 + x830 - 3.0 * x831 - 3.0 * x832 - 5.0 * x833 + 3.0 * x834 +
                     7.0 * x835) -
             x670 * (-AXG * x836 - AXG * x837 - AXG * x838 - MGs * x269 * x315 + x147 * x269 -
                     x148 * x269 - x161 * x58 - x291 * x40 + x341 * x40 + x589 - x595 + x597 -
                     x615 - x619 + x622 + x625 + x681 - x826 - x827 + x828 + x829 + x836 + x837 +
                     x838 + x839 - x840) +
             x709 *
                 (AXG * x684 + AXG * x686 - x105 * x685 - x117 * x691 + x118 * x691 - x126 * x25 +
                  x206 * x692 - x684 - x686 + 2.0 * x687 + x689 + x690 - x692 * x95 + x693 - x695 -
                  x696 - x697 - 6.0 * x698 + x699 * x701 + x701 * x702 + x705 + x708) +
             x709 * (x121 * x25 - x147 * x149 * x25 - x16 * x41 - x392 * x41 + x543 * x685 +
                     2.0 * x710 + 2.0 * x711 + x712 - x713 - 2.0 * x714 - x715 + 2.0 * x716 +
                     2.0 * x717 + 2.0 * x718 + x719 - 4.0 * x720 - 6.0 * x723 - 6.0 * x725 -
                     6.0 * x726 - 6.0 * x727 + x730 + x732) +
             x737 * (MUs * x734 + MXs * x734 - x122 * x736 + x148 * x461 + x151 * x728 + x539 +
                     x708 + 3.0 * x710 + 3.0 * x711 - x712 - 3.0 * x714 - 7.0 * x716 + 3.0 * x717 +
                     3.0 * x718 - 7.0 * x723 - 7.0 * x725 - 7.0 * x726 - 7.0 * x727 - x731 + x732 +
                     x733 - 2.0 * x735) -
             x737 * (AXG * x522 - x122 * x776 - x122 * x85 + x185 * x309 + x185 * x326 -
                     x185 * x356 + x203 * x29 - x203 * x309 + x203 * x356 + x25 * x358 -
                     x25 * x818 + x311 * x95 - x522 - 9.0 * x529 + x713 + x715 - x719 + 3.0 * x739 +
                     3.0 * x740 - 7.0 * x742 - 7.0 * x743 - 5.0 * x850 + x851) -
             x737 * (-AXG * x846 - AXG * x847 - AXG * x848 + MUs * x770 - MUs * x849 + MXs * x770 -
                     MXs * x849 + x115 * x803 + x151 * x687 + x538 - x687 - x689 - x690 - x693 -
                     x694 + x695 + x696 + x697 - x703 + x704 + x779 + x846 + x847 + x848) +
             x758 * (AXG * x748 - AXG * x749 - AXG * x751 + AXG * x753 - AXG * x756 + x117 * x399 +
                     x117 * x543 + x151 * x754 - x180 * x667 + x180 * x668 - x538 + x539 -
                     x685 * x752 + x687 - x698 + x746 - x748 + x749 + x751 - x753 - 5.0 * x754 +
                     x756 + x757) +
             x758 * (AXG * x762 + AXG * x763 + AXG * x764 + x148 * x400 + x161 * x203 - x40 * x462 +
                     x40 * x85 + x403 * x765 + x710 + x711 - x714 - x716 + x717 + x718 - x720 -
                     x726 - x727 + x730 - x735 + x760 - x762 - x763 - x764 + x767) -
             x758 *
                 (-MGs * x769 - MUs * x195 - x147 * x412 + x148 * x412 + x151 * x854 + x180 * x765 -
                  x185 * x328 + x397 - 3.0 * x525 + 5.0 * x526 + 4.0 * x739 + 4.0 * x740 -
                  7.0 * x741 - 8.0 * x742 - x850 + x851 + 2.0 * x852 - 3.0 * x853 - 5.0 * x854 +
                  3.0 * x855 + 3.0 * x856 + 3.0 * x857 + 7.0 * x858 - 2.0 * x859 - 3.0 * x860) +
             x863 * (AXG * x653 + AXG * x663 + AXG * x665 - AXG * x861 - AXG * x862 + x594 + x599 -
                     x600 + x601 - x602 + x611 + x649 - x653 - x663 - x665 + x669 - x839 + x840 +
                     x861 + x862) -
             x863 * (AXG * x833 + x133 * x351 + x147 * x351 - x148 * x351 - x161 * x351 -
                     x34 * x351 + x366 - x633 - x634 - x635 + x637 + x638 + x639 - x640 - x641 -
                     x642 + x643 + x644 + x645 - x820 + x825 - x831 - x832 - x833 + x834 + x835) +
             x867 * (-AXG * x854 + MGs * x195 - MGs * x393 + MGs * x395 - MGs * x468 + x147 * x180 -
                     x148 * x180 - x161 * x180 + x188 + x744 + x853 + x854 - x856 - x857 - x858 +
                     x859 + x860) +
             x867 * (AXG * x747 + AXG * x754 - AXG * x864 - AXG * x865 + AXG * x866 - MGs * x201 +
                     MGs * x469 + x161 * x403 + x722 - x723 + x724 - x725 + x729 + x733 - x754 +
                     x757 + x767 - x852 - x855 + x864 + x865 - x866)) *
            Denom(768.0 * Nc * s * (MUs - u) * pow(Pi, 2)));
    //}
  }

  return ret.real();
}
