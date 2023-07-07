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
#define D00(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd00, a, b, c, d, e, f, g, h, i, j)
#define D11(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd11, a, b, c, d, e, f, g, h, i, j)
#define D12(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd12, a, b, c, d, e, f, g, h, i, j)
#define D22(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd22, a, b, c, d, e, f, g, h, i, j)
#define D13(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd13, a, b, c, d, e, f, g, h, i, j)
#define D23(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd23, a, b, c, d, e, f, g, h, i, j)
#define D33(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd33, a, b, c, d, e, f, g, h, i, j)

ComplexType ME_us_box_QGqq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2,
                           double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){
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

    auto Denom = [](auto a) { return 1. / a; };

    //*

#define syFC1 C0(MUs, 0, MUs + MXs - s - t, MGs, Mqis, Mqis)
#define syFC2 D0(MUs, 0, MXs, 0, MUs + MXs - s - t, t, MGs, Mqis, Mqis, MQis)
#define syFC3 C1(MUs, MUs + MXs - s - t, 0, Mqis, MGs, Mqis)
#define syFC4 D1(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC5 C2(0, MUs + MXs - s - t, MXs, MQis, MGs, Mqis)
#define syFC6 C2(MUs, MUs + MXs - s - t, 0, Mqis, MGs, Mqis)
#define syFC7 D2(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC8 D3(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC9 D00(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC10 D11(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC11 D12(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC12 D13(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC13 D22(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC14 D23(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)
#define syFC15 D33(MUs, MUs + MXs - s - t, MXs, t, 0, 0, Mqis, MGs, Mqis, MQis)

    _EPS0_(
        ret, auto x0 = L * jRGp; auto x1 = iLGp * x0; auto x2 = kR * x1;
        auto x3 = pow(s, 3); auto x4 = UU * u; auto x5 = x3 * x4;
        auto x6 = x2 * x5; auto x7 = 3.0 * x6; auto x8 = pow(s, 2);
        auto x9 = kR * x8; auto x10 = pow(u, 2); auto x11 = UU * x1;
        auto x12 = x10 * x11; auto x13 = 7.0 * x12; auto x14 = x13 * x9;
        auto x15 = MXs * x3; auto x16 = x11 * x15; auto x17 = kR * x16;
        auto x18 = 3.0 * x17; auto x19 = AXG * x18; auto x20 = AXG * x7;
        auto x21 = R * jLGp; auto x22 = MG * iLGp; auto x23 = x21 * x22;
        auto x24 = UU * x23; auto x25 = MX * x3; auto x26 = kR * x25;
        auto x27 = x24 * x26; auto x28 = 3.0 * x27; auto x29 = UU * iRGp;
        auto x30 = x21 * x29; auto x31 = Mqi * x30; auto x32 = x26 * x31;
        auto x33 = 3.0 * x32; auto x34 = AXG * x14; auto x35 = 2.0 * SS;
        auto x36 = u * x35; auto x37 = MX * x23; auto x38 = x36 * x37;
        auto x39 = x38 * x9; auto x40 = x4 * x8; auto x41 = x37 * x40;
        auto x42 = 7.0 * kR; auto x43 = AXG * x40; auto x44 = x10 * x35;
        auto x45 = x1 * x44; auto x46 = x45 * x9; auto x47 = MUs * SS;
        auto x48 = pow(MX, 5); auto x49 = Mqi * x48; auto x50 = x47 * x49;
        auto x51 = 4.0 * x1; auto x52 = kL * x51; auto x53 = x50 * x52;
        auto x54 = x2 * x40; auto x55 = 3.0 * x54; auto x56 = AXG * x46;
        auto x57 = MXs * x8; auto x58 = UU * x57; auto x59 = 3.0 * x58;
        auto x60 = x2 * x59; auto x61 = AXG * x59; auto x62 = x2 * x61;
        auto x63 = t * x62; auto x64 = x2 * x43; auto x65 = 3.0 * x64;
        auto x66 = t * x65; auto x67 = UU * t; auto x68 = x37 * x67;
        auto x69 = 3.0 * x9; auto x70 = x68 * x69; auto x71 = MX * t;
        auto x72 = x31 * x71; auto x73 = AXG * x38; auto x74 = AXG * x72;
        auto x75 = AXG * x70 + t * x55 - t * x60 + x46 + x53 - x56 + x63 - x66 -
                   x69 * x72 + x69 * x74 - x70 + x73 * x9;
        auto x76 = 2.0 * x6; auto x77 = 4.0 * x12; auto x78 = x77 * x9;
        auto x79 = 2.0 * x11; auto x80 = kR * x15 * x79; auto x81 = AXG * x76;
        auto x82 = 2.0 * x21; auto x83 = x29 * x82; auto x84 = Mqi * x25;
        auto x85 = kR * x84; auto x86 = x83 * x85; auto x87 = 4.0 * AXG;
        auto x88 = x10 * x87; auto x89 = x11 * x88; auto x90 = x89 * x9;
        auto x91 = AXG * x83; auto x92 = x85 * x91; auto x93 = MUs * x35;
        auto x94 = x1 * x93; auto x95 = AXG * x49; auto x96 = x94 * x95;
        auto x97 = kL * x96; auto x98 = 4.0 * x37; auto x99 = kR * x40;
        auto x100 = x37 * x87; auto x101 = AXG * x80; auto x102 = x22 * x82;
        auto x103 = UU * x102; auto x104 = x103 * x26;
        auto x105 = AXG * x104 + x101 - x104; auto x106 = s * x10;
        auto x107 = x106 * x35; auto x108 = UU * s; auto x109 = x10 * x108;
        auto x110 = 4.0 * x109; auto x111 = s * t; auto x112 = x111 * x4;
        auto x113 = 3.0 * x112; auto x114 = AXG * x107; auto x115 = x109 * x87;
        auto x116 = SS * u; auto x117 = MXs * s; auto x118 = 4.0 * x117;
        auto x119 = x116 * x118; auto x120 = MXs * x111; auto x121 = UU * x120;
        auto x122 = 3.0 * x121; auto x123 = x117 * x4; auto x124 = 7.0 * x123;
        auto x125 = AXG * x117; auto x126 = x125 * x36; auto x127 = AXG * x122;
        auto x128 = AXG * x123; auto x129 = 11.0 * x128; auto x130 = AXG * x113;
        auto x131 = pow(MX, 6); auto x132 = UU * x131; auto x133 = x132 * x87;
        auto x134 = iRGp * x21; auto x135 = kL * x134; auto x136 = pow(MX, 4);
        auto x137 = 4.0 * x136; auto x138 = x116 * x137; auto x139 = MUs * x44;
        auto x140 = x137 * x47; auto x141 = x135 * x140; auto x142 = MXs * x10;
        auto x143 = 6.0 * SS; auto x144 = x142 * x143; auto x145 = AXG * x136;
        auto x146 = x145 * x93; auto x147 = AXG * x139; auto x148 = x142 * x87;
        auto x149 = SS * x148; auto x150 = MXs * x47; auto x151 = u * x150;
        auto x152 = 6.0 * x151; auto x153 = x145 * x36; auto x154 = 8.0 * x4;
        auto x155 = x145 * x154; auto x156 = x151 * x87;
        auto x157 = x133 * x135 + x135 * x138 - x135 * x139 - x135 * x144 +
                    x135 * x146 + x135 * x147 + x135 * x149 + x135 * x152 -
                    x135 * x153 - x135 * x155 - x135 * x156 - x141;
        auto x158 = pow(u, 3); auto x159 = x158 * x35; auto x160 = x134 * x159;
        auto x161 = x108 * x136; auto x162 = x135 * x161; auto x163 = x2 * x58;
        auto x164 = t * x44; auto x165 = AXG * x159; auto x166 = AXG * x58;
        auto x167 = x166 * x2; auto x168 = x112 * x2; auto x169 = x121 * x2;
        auto x170 = MXs * t; auto x171 = x170 * x93; auto x172 = kL * x29;
        auto x173 = x172 * x21; auto x174 = MUs * x36; auto x175 = t * x174;
        auto x176 = x170 * x36; auto x177 = AXG * x137; auto x178 = t * x173;
        auto x179 = AXG * x121; auto x180 = x179 * x2; auto x181 = AXG * x112;
        auto x182 = x181 * x2; auto x183 = AXG * x170; auto x184 = x154 * x183;
        auto x185 = iRGp * x82; auto x186 = kL * x185;
        auto x187 = x108 * x131 * x186; auto x188 = x161 * x186;
        auto x189 = MQis * x188; auto x190 = pow(MU, 4); auto x191 = x190 * x35;
        auto x192 = x117 * x135; auto x193 = MG * x0; auto x194 = iRGp * x193;
        auto x195 = x194 * x48; auto x196 = 2.0 * kL; auto x197 = x108 * x196;
        auto x198 = x1 * x108; auto x199 = x198 * x49; auto x200 = s * x135;
        auto x201 = x200 * x93; auto x202 = MUs * x161; auto x203 = x186 * x202;
        auto x204 = 4.0 * x30; auto x205 = x204 * x85; auto x206 = 2.0 * AXG;
        auto x207 = kL * x206; auto x208 = x199 * x207; auto x209 = s * x145;
        auto x210 = pow(MX, 3); auto x211 = x194 * x210; auto x212 = kL * s;
        auto x213 = Mqi * x102; auto x214 = kL * x213; auto x215 = x32 * x87;
        auto x216 = AXG * x161; auto x217 = Mqi * x210; auto x218 = x1 * x217;
        auto x219 = x218 * x93; auto x220 = AXG * s; auto x221 = x125 * x93;
        auto x222 = MQis * x134; auto x223 = kL * x222; auto x224 = Mqi * x23;
        auto x225 = kL * x221; auto x226 = x224 * x93; auto x227 = kL * x117;
        auto x228 = x117 * x223 * x93 + x226 * x227; auto x229 = 2.0 * x163;
        auto x230 = kR * x51; auto x231 = x108 * x51; auto x232 = kL * x217;
        auto x233 = u * x111; auto x234 = SS * x233; auto x235 = x210 * x23;
        auto x236 = 4.0 * x235; auto x237 = kR * x108; auto x238 = 4.0 * x134;
        auto x239 = x237 * x238; auto x240 = MX * x103; auto x241 = x240 * x9;
        auto x242 = MX * Mqi; auto x243 = x242 * x9; auto x244 = x161 * x2;
        auto x245 = x111 * x51; auto x246 = x123 * x2; auto x247 = kR * s;
        auto x248 = x116 * x37; auto x249 = 4.0 * x248; auto x250 = x116 * x238;
        auto x251 = x242 * x250; auto x252 = x198 * x87;
        auto x253 = x108 * x134; auto x254 = kR * x253 * x87;
        auto x255 = s * x47; auto x256 = x255 * x52; auto x257 = x23 * x47;
        auto x258 = MX * x257; auto x259 = 4.0 * x258; auto x260 = u * x237;
        auto x261 = x238 * x47; auto x262 = x242 * x261; auto x263 = u * x242;
        auto x264 = MX * x24; auto x265 = x264 * x9; auto x266 = x2 * x93;
        auto x267 = x233 * x35; auto x268 = x111 * x93; auto x269 = x172 * x193;
        auto x270 = 2.0 * MX * x269; auto x271 = MX * x194;
        auto x272 = x271 * x36; auto x273 = kR * x240; auto x274 = x271 * x93;
        auto x275 = x212 * x274; auto x276 = x37 * x93; auto x277 = MX * x102;
        auto x278 = x260 * x277; auto x279 = u * x108; auto x280 = 5.0 * x25;
        auto x281 = x24 * x280; auto x282 = kR * x281; auto x283 = 5.0 * x17;
        auto x284 = 5.0 * t; auto x285 = x163 * x284; auto x286 = MX * x31;
        auto x287 = x284 * x286 * x9; auto x288 = -AXG * x287;
        auto x289 = x285 + x287 + x288 + x56 - x76 + x81;
        auto x290 = x280 * x31; auto x291 = kR * x290;
        auto x292 = -AXG * x291 + x283 + x291; auto x293 = x10 * x79;
        auto x294 = x293 * x9; auto x295 = s * x141; auto x296 = 2.0 * t;
        auto x297 = x68 * x9; auto x298 = 5.0 * x297; auto x299 = x1 * x47;
        auto x300 = x299 * x87; auto x301 = x203 - x296 * x54;
        auto x302 = AXG * x294 - AXG * x295 - AXG * x298 - s * x232 * x300 -
                    x167 * x284 + x217 * x256 - x294 + x295 + x296 * x64 +
                    x298 + x301 - x46;
        auto x303 = 4.0 * x17 + x205 - x215 + x289; auto x304 = pow(MXs, 2);
        auto x305 = pow(MUs, 2); auto x306 = x305 * x35; auto x307 = t * x83;
        auto x308 = x243 * x307; auto x309 = 3.0 * x40; auto x310 = kR * x37;
        auto x311 = x242 * x94; auto x312 = 3.0 * x43; auto x313 = x1 * x242;
        auto x314 = AXG * x27; auto x315 = -x27 + x314; auto x316 = t * x1;
        auto x317 = x316 * x93; auto x318 = x316 * x36 * x9 - x317 * x9;
        auto x319 = x102 * x67; auto x320 = MX * x319 * x9;
        auto x321 = -AXG * x320 + x140 * x223 - x146 * x223 - x284 * x54 +
                    x318 + x320 - x53 - x63 + x66 + x97;
        auto x322 = x133 * x222; auto x323 = x136 * x8;
        auto x324 = x210 * x269 * x8; auto x325 = x11 * x217;
        auto x326 = x325 * x8; auto x327 = kL * x326; auto x328 = x135 * x57;
        auto x329 = AXG * x57; auto x330 = x135 * x329; auto x331 = MX * x8;
        auto x332 = x193 * x331; auto x333 = iRGp * x332;
        auto x334 = x333 * x47; auto x335 = Mqi * x331; auto x336 = kL * x335;
        auto x337 = x299 * x336; auto x338 = x1 * x116 * x336;
        auto x339 = AXG * x332; auto x340 = iRGp * x339;
        auto x341 = 2.0 * syFC13; auto x342 = SS * x106; auto x343 = AXG * x342;
        auto x344 = u * x299; auto x345 = x247 * x344; auto x346 = x116 * x117;
        auto x347 = u * x47; auto x348 = x200 * x347; auto x349 = kR * x299;
        auto x350 = x125 * x135; auto x351 = 2.0 * syFC6;
        auto x352 = x29 * x339; auto x353 = t * x352; auto x354 = 3.0 * x353;
        auto x355 = 7.0 * x30; auto x356 = x10 * x355; auto x357 = x356 * x8;
        auto x358 = x15 * x30; auto x359 = 3.0 * x358; auto x360 = x193 * x29;
        auto x361 = x25 * x360; auto x362 = 3.0 * x361; auto x363 = x11 * x84;
        auto x364 = 3.0 * x363; auto x365 = x271 * x40; auto x366 = x271 * x43;
        auto x367 = 3.0 * x134; auto x368 = x367 * x5; auto x369 = t * x134;
        auto x370 = t * x367 * x43; auto x371 = x369 * x61 - x370;
        auto x372 = AXG * x359 - AXG * x368 + x368 + x371;
        auto x373 = x134 * x44; auto x374 = x373 * x8; auto x375 = t * x40;
        auto x376 = x29 * x332; auto x377 = t * x376; auto x378 = Mqi * x11;
        auto x379 = x331 * x378; auto x380 = t * x379; auto x381 = 3.0 * x380;
        auto x382 = -AXG * x374 + AXG * x381 + x367 * x375 - x369 * x59 + x374 -
                    3.0 * x377 - x381;
        auto x383 = kL * syFC10; auto x384 = x210 * x360;
        auto x385 = x384 * x88; auto x386 = t * x93; auto x387 = -x218 * x386;
        auto x388 = x159 * x194; auto x389 = x134 * x136;
        auto x390 = x389 * x93; auto x391 = x134 * x191; auto x392 = x218 * x36;
        auto x393 = x360 * x48; auto x394 = t * x393; auto x395 = t * x11;
        auto x396 = x49 * x87; auto x397 = x177 * x30; auto x398 = MUs * x397;
        auto x399 = 6.0 * x134; auto x400 = x170 * x347; auto x401 = x1 * x175;
        auto x402 = AXG * t; auto x403 = x154 * x218; auto x404 = MUs * x134;
        auto x405 = pow(u, 4); auto x406 = x35 * x405; auto x407 = x134 * x406;
        auto x408 = x158 * x47; auto x409 = SS * x158; auto x410 = 8.0 * x409;
        auto x411 = x134 * x410; auto x412 = x134 * x87; auto x413 = x158 * x30;
        auto x414 = x413 * x87; auto x415 = MXs * x134; auto x416 = AXG * x415;
        auto x417 = x143 * x158; auto x418 = x1 * x159; auto x419 = x1 * x165;
        auto x420 = -AXG * x407 - MXs * x411 + MXs * x414 - x238 * x408 -
                    x242 * x418 + x242 * x419 + x407 + x408 * x412 +
                    x416 * x417;
        auto x421 = -MX * x388 - t * x390 + t * x392 + t * x398 + x133 * x369 +
                    x165 * x271 - x170 * x391 - x184 * x404 - x211 * x386 +
                    x242 * x401 + x387 + x394 * x87 + x395 * x396 +
                    x399 * x400 - x402 * x403 + x420;
        auto x422 = 7.0 * x216; auto x423 = t * x160; auto x424 = -t * x414;
        auto x425 = x10 * x47; auto x426 = x238 * x425; auto x427 = SS * x142;
        auto x428 = t * x427; auto x429 = -x238 * x428;
        auto x430 = x134 * x140 * x220; auto x431 = x30 * x88;
        auto x432 = MUs * x431; auto x433 = t * x432; auto x434 = x142 * x30;
        auto x435 = x402 * x434; auto x436 = 12.0 * x435;
        auto x437 = -x164 * x313; auto x438 = x242 * x395;
        auto x439 = x438 * x88; auto x440 = -t * x426 + x423 + x424 + x429 +
                                            x430 + x433 + x436 + x437 + x439;
        auto x441 = -3.0 * x199 + x404 * x422 + x440; auto x442 = x36 * x389;
        auto x443 = x145 * x4; auto x444 = x369 * x443;
        auto x445 = t * x442 - 12.0 * x444; auto x446 = x131 * x253;
        auto x447 = 3.0 * x161; auto x448 = AXG * x199;
        auto x449 = -x404 * x447 - 3.0 * x446 + 11.0 * x448;
        auto x450 = x190 * x36; auto x451 = 8.0 * x134; auto x452 = x255 * x451;
        auto x453 = x211 * x36; auto x454 = s * x217; auto x455 = x299 * x454;
        auto x456 = x360 * x71; auto x457 = x154 * x211;
        auto x458 = t * x453 - x136 * x452 - x164 * x271 + x175 * x271 +
                    x300 * x454 + x369 * x450 - x402 * x457 - 6.0 * x455 +
                    x456 * x88;
        auto x459 = x47 * x57; auto x460 = x238 * x459; auto x461 = x404 * x61;
        auto x462 = x111 * x373; auto x463 = x136 * x30; auto x464 = 3.0 * x463;
        auto x465 = x10 * x30; auto x466 = 3.0 * x111; auto x467 = x335 * x94;
        auto x468 = 3.0 * x384; auto x469 = x134 * x93; auto x470 = x329 * x469;
        auto x471 = x134 * x36; auto x472 = -x120 * x471; auto x473 = x120 * x4;
        auto x474 = x145 * x355; auto x475 = AXG * x356;
        auto x476 = x111 * x475; auto x477 = AXG * x120; auto x478 = x134 * x4;
        auto x479 = x477 * x478; auto x480 = 7.0 * x384; auto x481 = AXG * x480;
        auto x482 = AXG * x325; auto x483 = 7.0 * x482; auto x484 = 7.0 * x181;
        auto x485 = x111 * x464 + x111 * x468 - x111 * x474 - x111 * x481 -
                    x111 * x483 - x113 * x271 - x113 * x313 - x113 * x404 +
                    x130 * x404 - x233 * x469 + x271 * x484 + x313 * x484 +
                    x325 * x466 - x399 * x473 + x462 + x465 * x466 + x467 -
                    x470 + x472 - x476 + 14.0 * x479;
        auto x486 = x136 * x355; auto x487 = 11.0 * x8; auto x488 = x145 * x30;
        auto x489 = x10 * x83; auto x490 = x111 * x489; auto x491 = 4.0 * x334;
        auto x492 = x47 * x51; auto x493 = x335 * x492; auto x494 = AXG * x489;
        auto x495 = MUs * x112; auto x496 = x185 * x495;
        auto x497 = -x261 * x329; auto x498 = x238 * x473;
        auto x499 = x111 * x488; auto x500 = x181 * x185;
        auto x501 = x154 * x477; auto x502 = 2.0 * x271;
        auto x503 = x112 * x502; auto x504 = 2.0 * x1; auto x505 = x112 * x242;
        auto x506 = x504 * x505; auto x507 = x300 * x335;
        auto x508 = 6.0 * x384; auto x509 = x111 * x508; auto x510 = 6.0 * x271;
        auto x511 = 6.0 * x1; auto x512 = x181 * x242; auto x513 = x185 * x58;
        auto x514 = MUs * x513; auto x515 = x166 * x185;
        auto x516 = -MUs * x515 + x460 + x514;
        auto x517 = -AXG * x509 + MUs * x500 + 2.0 * x111 * x384 - x111 * x494 +
                    x134 * x501 + x181 * x510 + x490 + x491 + x493 - x496 +
                    x497 - x498 - 6.0 * x499 - x503 - x506 - x507 +
                    x511 * x512 + x516;
        auto x518 = kL * syFC11; auto x519 = x238 * x409;
        auto x520 = x116 * x190; auto x521 = x238 * x520; auto x522 = 8.0 * AXG;
        auto x523 = x413 * x522; auto x524 = t * x451; auto x525 = x255 * x389;
        auto x526 = 4.0 * x194; auto x527 = x210 * x526;
        auto x528 = x116 * x527; auto x529 = t * x522; auto x530 = 12.0 * x134;
        auto x531 = SS * x10; auto x532 = x526 * x71; auto x533 = Mqi * x531;
        auto x534 = MX * x533; auto x535 = t * x51; auto x536 = x10 * x522;
        auto x537 = 8.0 * x299; auto x538 = x217 * x220; auto x539 = x211 * x4;
        auto x540 = 16.0 * AXG; auto x541 = t * x540; auto x542 = x218 * x4;
        auto x543 = kL * syFC12; auto x544 = SS * x405; auto x545 = x544 * x87;
        auto x546 = x408 * x451; auto x547 = 16.0 * x409;
        auto x548 = AXG * x132; auto x549 = 12.0 * x409; auto x550 = MX * x526;
        auto x551 = x409 * x51; auto x552 = SS * x190; auto x553 = x170 * x552;
        auto x554 = x194 * x87; auto x555 = MX * x554; auto x556 = x1 * x87;
        auto x557 = x409 * x556; auto x558 = x116 * x217; auto x559 = t * x488;
        auto x560 = 8.0 * MUs; auto x561 = x47 * x527; auto x562 = x217 * x47;
        auto x563 = x522 * x531; auto x564 = x347 * x51; auto x565 = t * x564;
        auto x566 = 16.0 * MUs; auto x567 = x531 * x8; auto x568 = 5.0 * x134;
        auto x569 = x5 * x568; auto x570 = x465 * x487; auto x571 = 6.0 * x369;
        auto x572 = 5.0 * x358; auto x573 = x8 * x88; auto x574 = SS * x134;
        auto x575 = x369 * x58; auto x576 = 5.0 * x361; auto x577 = x116 * x333;
        auto x578 = 6.0 * x11; auto x579 = t * x335; auto x580 = x578 * x579;
        auto x581 = 11.0 * x43; auto x582 = 5.0 * x363;
        auto x583 = AXG * x582 - x572 - x582; auto x584 = 6.0 * x158;
        auto x585 = x253 * x584; auto x586 = x111 * x238;
        auto x587 = 6.0 * x463; auto x588 = 6.0 * x465; auto x589 = x217 * x578;
        auto x590 = MUs * x399; auto x591 = -x233 * x261;
        auto x592 = AXG * x111; auto x593 = AXG * x384; auto x594 = 14.0 * x111;
        auto x595 = 14.0 * x181; auto x596 = x463 * x8; auto x597 = 5.0 * x325;
        auto x598 = 9.0 * x325; auto x599 = kL * syFC14; auto x600 = x185 * x3;
        auto x601 = x4 * x600; auto x602 = x10 * x204; auto x603 = x15 * x83;
        auto x604 = AXG * x601; auto x605 = x79 * x84; auto x606 = AXG * x605;
        auto x607 = AXG * x603; auto x608 = 2.0 * x361;
        auto x609 = AXG * x608 + x354 + x607 - x608; auto x610 = kL * syFC15;
        auto x611 = x134 * x552; auto x612 = AXG * x446; auto x613 = 5.0 * x216;
        auto x614 = -x199 + x404 * x613 + 9.0 * x612; auto x615 = AXG * x467;
        auto x616 = x36 * x8; auto x617 = x369 * x616; auto x618 = -x617;
        auto x619 = x375 * x568; auto x620 = 2.0 * x376; auto x621 = t * x620;
        auto x622 = x284 * x379; auto x623 = x296 * x352;
        auto x624 = AXG * x622; auto x625 = AXG * x361 - x361;
        auto x626 = kL * syFC4; auto x627 = 10.0 * x531;
        auto x628 = x389 * x627; auto x629 = x190 * x373;
        auto x630 = AXG * x629; auto x631 = x10 * x143; auto x632 = x134 * x145;
        auto x633 = x631 * x632; auto x634 = x145 * x465;
        auto x635 = 12.0 * x634; auto x636 = UU * x148;
        auto x637 = MQis * x160 - x139 * x224 - x144 * x224 + x147 * x224 +
                    x149 * x224 + x159 * x224 - x165 * x222 - x165 * x224 +
                    x171 * x222 + x224 * x636 + x420 + x629 - x630 - x633 -
                    x635;
        auto x638 = x118 * x47; auto x639 = x67 * x88;
        auto x640 = x164 * x224 - x221 * x222 + x222 * x638 - x224 * x639 +
                    6.0 * x299 * x538 - x454 * x537 - 6.0 * x525;
        auto x641 = x131 * x238; auto x642 = x47 * x641;
        auto x643 = x137 * x611; auto x644 = pow(MX, 8); auto x645 = x30 * x644;
        auto x646 = x645 * x87; auto x647 = x138 * x222; auto x648 = AXG * x131;
        auto x649 = x469 * x648; auto x650 = x145 * x391;
        auto x651 = pow(MX, 7); auto x652 = x651 * x87; auto x653 = x378 * x652;
        auto x654 = x133 * x134; auto x655 = MUs * x654; auto x656 = x116 * x49;
        auto x657 = x51 * x656; auto x658 = u * x389; auto x659 = x47 * x658;
        auto x660 = 14.0 * x659; auto x661 = x1 * x36; auto x662 = x661 * x95;
        auto x663 = x133 * x224; auto x664 = x154 * x95; auto x665 = x1 * x664;
        auto x666 = x347 * x451; auto x667 = x145 * x666;
        auto x668 = x155 * x404; auto x669 = x152 * x222;
        auto x670 = x153 * x222; auto x671 = x140 * x224;
        auto x672 = x155 * x222; auto x673 = 6.0 * x344;
        auto x674 = x217 * x673; auto x675 = x146 * x224;
        auto x676 = x344 * x87; auto x677 = x217 * x676;
        auto x678 = x156 * x222; auto x679 = 4.0 * x358; auto x680 = t * x185;
        auto x681 = x40 * x680; auto x682 = 4.0 * x363; auto x683 = x363 * x87;
        auto x684 = kL * syFC7; auto x685 = x471 * x648;
        auto x686 = x218 * x631; auto x687 = x4 * x648; auto x688 = x530 * x687;
        auto x689 = x134 * x47; auto x690 = 14.0 * x142;
        auto x691 = x689 * x690; auto x692 = x139 * x222;
        auto x693 = x144 * x222; auto x694 = u * x143; auto x695 = x190 * x694;
        auto x696 = x415 * x695; auto x697 = x147 * x222;
        auto x698 = x139 * x313; auto x699 = x149 * x222;
        auto x700 = x148 * x30; auto x701 = MQis * x700;
        auto x702 = x415 * x520; auto x703 = x702 * x87;
        auto x704 = x138 * x224; auto x705 = SS * x217; auto x706 = x705 * x88;
        auto x707 = x1 * x706; auto x708 = x325 * x88; auto x709 = MUs * x700;
        auto x710 = 10.0 * AXG; auto x711 = x142 * x710;
        auto x712 = x689 * x711; auto x713 = x152 * x224;
        auto x714 = x147 * x313; auto x715 = x153 * x224;
        auto x716 = x155 * x224; auto x717 = x156 * x224;
        auto x718 = kL * syFC8; auto x719 = x335 * x79; auto x720 = 2.0 * x352;
        auto x721 = x116 * x550; auto x722 = x116 * x242;
        auto x723 = x51 * x722; auto x724 = x211 * x87; auto x725 = MUs * x1;
        auto x726 = x447 * x725; auto x727 = x253 * x49; auto x728 = t * x418;
        auto x729 = x11 * x158; auto x730 = x729 * x87; auto x731 = t * x730;
        auto x732 = x425 * x51; auto x733 = t * x732; auto x734 = x428 * x51;
        auto x735 = x136 * x299; auto x736 = s * x735; auto x737 = x235 * x36;
        auto x738 = x1 * x140; auto x739 = AXG * x738; auto x740 = s * x739;
        auto x741 = MUs * x395; auto x742 = x741 * x88; auto x743 = x11 * x142;
        auto x744 = x402 * x743; auto x745 = 12.0 * x744;
        auto x746 = x134 * x242; auto x747 = x164 * x746;
        auto x748 = x217 * x399; auto x749 = x255 * x748; auto x750 = x72 * x88;
        auto x751 = x154 * x235; auto x752 = x1 * x450; auto x753 = x1 * x136;
        auto x754 = x36 * x753; auto x755 = x316 * x443;
        auto x756 = x217 * x369;
        auto x757 = -AXG * x154 * x756 + t * x752 + t * x754 - 12.0 * x755;
        auto x758 = t * x737 - x164 * x37 + x175 * x37 + x217 * x255 * x412 +
                    x37 * x639 - x402 * x751 + x728 - x731 - x733 - x734 -
                    8.0 * x736 + x740 + x742 + x745 - x747 - x749 + x750 + x757;
        auto x759 = kR * syFC10; auto x760 = x408 * x51; auto x761 = x1 * x410;
        auto x762 = MXs * x761; auto x763 = x158 * x87; auto x764 = x299 * x763;
        auto x765 = MXs * x730; auto x766 = MXs * x1; auto x767 = AXG * x766;
        auto x768 = x417 * x767; auto x769 = x160 * x242;
        auto x770 = x134 * x165; auto x771 = x242 * x770; auto x772 = x23 * x67;
        auto x773 = x48 * x772; auto x774 = SS * x88; auto x775 = x210 * x24;
        auto x776 = x1 * x406; auto x777 = x1 * x191; auto x778 = x31 * x48;
        auto x779 = t * x87; auto x780 = x369 * x93; auto x781 = x134 * x175;
        auto x782 = -AXG * x776 + x133 * x316 + x170 * x673 - x170 * x777 +
                    x177 * x741 - x184 * x725 - x217 * x780 + x242 * x781 +
                    x36 * x756 - x386 * x753 + x776 + x778 * x779;
        auto x783 = -x159 * x37 + x165 * x37 - x235 * x386 - x235 * x774 -
                    x760 - x762 + x764 + x765 + x768 - x769 + x771 +
                    x773 * x87 - x775 * x88 + x782;
        auto x784 = x459 * x51; auto x785 = x11 * x136; auto x786 = 3.0 * x785;
        auto x787 = x111 * x786; auto x788 = x12 * x466;
        auto x789 = x122 * x725; auto x790 = 3.0 * x775;
        auto x791 = x111 * x790; auto x792 = x210 * x31; auto x793 = 3.0 * x792;
        auto x794 = x111 * x793; auto x795 = x233 * x94;
        auto x796 = x113 * x725; auto x797 = x473 * x511;
        auto x798 = x11 * x145; auto x799 = 7.0 * x798; auto x800 = x111 * x799;
        auto x801 = x130 * x725; auto x802 = x1 * x4; auto x803 = x477 * x802;
        auto x804 = 14.0 * x803; auto x805 = x127 * x725;
        auto x806 = x113 * x37; auto x807 = 7.0 * x775; auto x808 = AXG * x807;
        auto x809 = x111 * x808; auto x810 = x217 * x355;
        auto x811 = AXG * x810; auto x812 = x111 * x811; auto x813 = x37 * x484;
        auto x814 = x484 * x746; auto x815 = x329 * x94; auto x816 = -x815;
        auto x817 = x13 * x592; auto x818 = x113 * x746;
        auto x819 = x816 - x817 - x818; auto x820 = x111 * x45;
        auto x821 = x120 * x94; auto x822 = x120 * x661;
        auto x823 = x820 + x821 - x822; auto x824 = 2.0 * x725;
        auto x825 = x166 * x824; auto x826 = -x825; auto x827 = 2.0 * x785;
        auto x828 = x121 * x824; auto x829 = x103 * x210;
        auto x830 = x111 * x293; auto x831 = x112 * x824;
        auto x832 = x329 * x492; auto x833 = x4 * x51; auto x834 = x120 * x833;
        auto x835 = x111 * x798; auto x836 = 6.0 * x775; auto x837 = 6.0 * x792;
        auto x838 = 6.0 * x37; auto x839 = x335 * x412;
        auto x840 = x261 * x335 - x47 * x839; auto x841 = 2.0 * x58;
        auto x842 = x725 * x841; auto x843 = x112 * x185;
        auto x844 = -x112 * x277 - x242 * x843 + x830 + x842;
        auto x845 = -AXG * x830 + x1 * x501 + x111 * x217 * x83 + x111 * x827 +
                    x111 * x829 - x179 * x824 + x181 * x824 + x181 * x838 +
                    x399 * x512 - x592 * x836 - x592 * x837 + x784 + x828 -
                    x831 - x832 - x834 - 6.0 * x835 + x840 + x844;
        auto x846 = kR * syFC11; auto x847 = x51 * x520; auto x848 = x1 * x138;
        auto x849 = x238 * x49; auto x850 = x522 * x729; auto x851 = 8.0 * x1;
        auto x852 = x116 * x236; auto x853 = MUs * x12; auto x854 = x531 * x98;
        auto x855 = x238 * x534; auto x856 = x217 * x530;
        auto x857 = x347 * x98; auto x858 = x217 * x452; auto x859 = x235 * x4;
        auto x860 = x217 * x4; auto x861 = kR * syFC12; auto x862 = x158 * x537;
        auto x863 = x1 * x132; auto x864 = AXG * x863; auto x865 = 8.0 * x864;
        auto x866 = x409 * x412; auto x867 = x217 * x250;
        auto x868 = 12.0 * x344; auto x869 = x210 * x257;
        auto x870 = 4.0 * x869; auto x871 = x217 * x261;
        auto x872 = x238 * x347; auto x873 = x242 * x872; auto x874 = MUs * x4;
        auto x875 = x1 * x874; auto x876 = 16.0 * x875; auto x877 = 5.0 * x1;
        auto x878 = x5 * x877; auto x879 = x12 * x487; auto x880 = 5.0 * x16;
        auto x881 = SS * x1; auto x882 = 6.0 * x316; auto x883 = 4.0 * x331;
        auto x884 = 6.0 * x331; auto x885 = x772 * x884; auto x886 = x31 * x884;
        auto x887 = x116 * x87; auto x888 = 11.0 * x134; auto x889 = x242 * x43;
        auto x890 = x792 * x8; auto x891 = x826 + 5.0 * x890;
        auto x892 = kR * syFC14; auto x893 = kR * syFC15;
        auto x894 = x613 * x725; auto x895 = 9.0 * AXG; auto x896 = 10.0 * x8;
        auto x897 = MUs * x58; auto x898 = MQis * x1; auto x899 = Mqi * x194;
        auto x900 = AXG * x792; auto x901 = 14.0 * x900; auto x902 = x57 * x94;
        auto x903 = x111 * x12; auto x904 = x23 * x93; auto x905 = x331 * x904;
        auto x906 = -AXG * x905 + x181 * x277 + x840;
        auto x907 = -MUs * x103 * x331 + x113 * x899 - x130 * x899 +
                    9.0 * x181 * x746 + x819 + x902 + 5.0 * x903 + x906;
        auto x908 = kR * syFC4; auto x909 = x116 * x131; auto x910 = x51 * x909;
        auto x911 = x648 * x661; auto x912 = 12.0 * x1; auto x913 = x687 * x912;
        auto x914 = x299 * x690; auto x915 = x250 * x49;
        auto x916 = x695 * x766; auto x917 = x471 * x95;
        auto x918 = x139 * x746; auto x919 = x11 * x148;
        auto x920 = x116 * x766; auto x921 = x190 * x920;
        auto x922 = x87 * x921; auto x923 = u * x537; auto x924 = x145 * x923;
        auto x925 = x155 * x725; auto x926 = x134 * x664;
        auto x927 = MUs * x919; auto x928 = x299 * x711;
        auto x929 = MQis * x919 + x138 * x899 - x139 * x898 - x144 * x898 +
                    x147 * x898 + x149 * x898 + x152 * x899 - x153 * x899 -
                    x155 * x899 - x156 * x899 - x910 + x911 + x913 + x914 -
                    x915 - x916 + x917 + x918 + x922 + x924 + x925 + x926 -
                    x927 - x928;
        auto x930 = MQis * x395; auto x931 = Mqi * x360;
        auto x932 = Mqi * x388 - t * x177 * x931 - x165 * x899 + x171 * x898 +
                    x171 * x899 - x175 * x898 - x176 * x898 - x177 * x930 +
                    x184 * x898 + x757 + x782;
        auto x933 = 6.0 * x785; auto x934 = x790 * x8; auto x935 = x8 * x837;
        auto x936 = -x820; auto x937 = x177 * x8; auto x938 = 2.0 * x899;
        auto x939 = x112 * x938; auto x940 = x242 * x500;
        auto x941 = x841 * x898; auto x942 = 2.0 * x898; auto x943 = 2.0 * x166;
        auto x944 = x166 * x942 - x841 * x899 + x899 * x943 - x941;
        auto x945 = kR * syFC7; auto x946 = 7.0 * x785; auto x947 = kR * syFC8;
        auto x948 = x211 * x347; auto x949 = 6.0 * x948; auto x950 = x48 * x526;
        auto x951 = AXG * x195;
        auto x952 = -x686 - x691 - x698 + x707 + x708 + x709 + x712 + x714;
        auto x953 = x116 * x950 - x139 * x271 + x147 * x271 - x154 * x951 -
                    x211 * x631 + x211 * x774 - x347 * x724 - x36 * x951 -
                    x628 - x629 + x630 + x633 + x635 + x952;
        auto x954 = x116 * x641; auto x955 = -x685 - x688 + x696 - x703 + x954;
        auto x956 = -x642 - x643 + x646 + x649 + x650 + x653 + x655 + x657 +
                    x660 - x662 - x665 - x667 - x668 + x674 - x677;
        auto x957 = x360 * x652 - x47 * x950 + x93 * x951 + x955 + x956;
        auto x958 = x10 * x552; auto x959 = 20.0 * x531; auto x960 = x116 * x48;
        auto x961 = 8.0 * x194; auto x962 = x145 * x531;
        auto x963 = 12.0 * x531; auto x964 = 28.0 * x142;
        auto x965 = MUs * x522; auto x966 = 20.0 * x134; auto x967 = AXG * x142;
        auto x968 = x116 * x451; auto x969 = MXs * x968; auto x970 = x4 * x540;
        auto x971 = x299 * x88; auto x972 = x131 * x47; auto x973 = 8.0 * x611;
        auto x974 = x522 * x651; auto x975 = x47 * x48; auto x976 = 16.0 * x95;
        auto x977 = x145 * x478; auto x978 = AXG * x217; auto x979 = x4 * x57;
        auto x980 = x335 * x661; auto x981 = x242 * x40; auto x982 = kL * syFC2;
        auto x983 = x304 * x83; auto x984 = x469 * x57; auto x985 = x179 * x242;
        auto x986 = MQis * x513 - x112 * x213 + x213 * x58 + x462 - x490 +
                    x496 + x498 - x504 * x512 - x514 + x591;
        auto x987 = x384 * x8; auto x988 = 14.0 * AXG;
        auto x989 = x470 - x493 + x507; auto x990 = x589 * x8;
        auto x991 = x185 * x4; auto x992 = x23 * x960; auto x993 = x143 * x235;
        auto x994 = x235 * x347; auto x995 = AXG * x48; auto x996 = x23 * x36;
        auto x997 = x190 * x45; auto x998 = x1 * x631; auto x999 = x10 * x798;
        auto x1000 = x134 * x217;
        auto x1001 = AXG * x997 - x1000 * x631 + x134 * x706 + x145 * x998 +
                     x147 * x746 - x627 * x753 + x792 * x88 - x997 +
                     12.0 * x999;
        auto x1002 = -x10 * x993 + x1001 - x139 * x37 + x147 * x37 -
                     x154 * x23 * x995 - x87 * x994 - x914 - x918 - x922 +
                     x927 + x928 + 4.0 * x992 + 6.0 * x994 - x995 * x996;
        auto x1003 = x257 * x48; auto x1004 = x11 * x644;
        auto x1005 = x1 * x137; auto x1006 = u * x735; auto x1007 = x347 * x412;
        auto x1008 = x1004 * x87 - x1005 * x552 + 14.0 * x1006 - x1007 * x217 +
                     x133 * x725 + x145 * x777 + x31 * x652 + x347 * x748 +
                     x469 * x95 - x47 * x849 - x51 * x972 + x648 * x94;
        auto x1009 = -4.0 * x1003 + x1008 + x24 * x652 + x904 * x995 + x910 -
                     x911 - x913 + x915 + x916 - x917 - x924 - x925 - x926;
        auto x1010 = x1 * x552; auto x1011 = AXG * x451;
        auto x1012 = x1011 * x217; auto x1013 = x242 * x689;
        auto x1014 = 8.0 * x753; auto x1015 = x145 * x299;
        auto x1016 = x198 * x584; auto x1017 = x233 * x492;
        auto x1018 = 6.0 * x725; auto x1019 = AXG * x775;
        auto x1020 = x304 * x79; auto x1021 = 2.0 * x109;
        auto x1022 = AXG * x1021; auto x1023 = x899 * x93;
        auto x1024 = kR * syFC2;
        auto x1025 = x1008 - x133 * x898 - x133 * x899 - x138 * x898 +
                     x140 * x898 + x140 * x899 - x146 * x898 - x146 * x899 -
                     x152 * x898 + x153 * x898 + x155 * x898 + x156 * x898;
        auto x1026 = -MQis * x418 + x1001 + x139 * x899 + x144 * x899 -
                     x147 * x899 - x148 * x931 - x149 * x899 + x165 * x898 +
                     x760 + x762 - x764 - x765 - x768 + x769 - x771;
        auto x1027 = x88 * x931;
        auto x1028 = -AXG * x749 + t * x1027 - x164 * x898 - x164 * x899 +
                     x175 * x899 + x176 * x899 - x184 * x899 + x221 * x898 -
                     x638 * x898 - x728 + x731 + x733 + x734 + 6.0 * x736 -
                     x740 - x742 - x745 + x747 - x750 + x858 + x88 * x930;
        auto x1029 = x331 * x996; auto x1030 = 3.0 * x331 * x772;
        auto x1031 = 2.0 * x4; auto x1032 = x1 * x1031; auto x1033 = x15 * x568;
        auto x1034 = x194 * x280; auto x1035 = x84 * x877;
        auto x1036 = u * x600; auto x1037 = u * x680 * x8;
        auto x1038 = x579 * x877;
        auto x1039 = -AXG * x1036 - AXG * x1037 + AXG * x1038 +
                     t * x329 * x568 + x1036 + x1037 - x1038 - x284 * x333 +
                     x284 * x340;
        auto x1040 = 4.0 * x10; auto x1041 = MX * x360;
        auto x1042 = x1040 * x1041; auto x1043 = x1041 * x88;
        auto x1044 = AXG * x539; auto x1045 = 19.0 * AXG; auto x1046 = -x866;
        auto x1047 = x399 * x425; auto x1048 = -x1047;
        auto x1049 = x1011 * x427; auto x1050 = AXG * x1047;
        auto x1051 = -x531 * x550;
        auto x1052 = x1046 + x1048 + x1049 + x1050 + x1051 + x519;
        auto x1053 =
            -MUs * x602 + x1045 * x434 + x1052 - x242 * x77 + x242 * x89 + x432;
        auto x1054 = x211 * x694; auto x1055 = x116 * x211;
        auto x1056 = -x1055 * x87; auto x1057 = x1054 + x1056;
        auto x1058 = x158 * x204 - x414; auto x1059 = x134 * x171;
        auto x1060 = x242 * x88;
        auto x1061 = x1060 * x881 + x271 * x774 - x51 * x534;
        auto x1062 = x1059 + x1061 - x427 * x530; auto x1063 = s * x383;
        auto x1064 = x389 * x4; auto x1065 = MXs * x552;
        auto x1066 = 7.0 * x542; auto x1067 = x415 * x874;
        auto x1068 = AXG * x393; auto x1069 = 26.0 * x977;
        auto x1070 = AXG * x542; auto x1071 = x191 * x415;
        auto x1072 = AXG * x1071; auto x1073 = x134 * x151;
        auto x1074 = x1073 * x710; auto x1075 = -x1074;
        auto x1076 = x416 * x874; auto x1077 = x1072 + x1075 - 11.0 * x1076;
        auto x1078 = x218 * x694; auto x1079 = x211 * x47;
        auto x1080 = 6.0 * x1079; auto x1081 = x134 * x450;
        auto x1082 = AXG * x1081; auto x1083 = -x1082; auto x1084 = x134 * x138;
        auto x1085 = AXG * x1084; auto x1086 = -x1085;
        auto x1087 = x1081 + x1083 + x1086; auto x1088 = x347 * x550;
        auto x1089 = x1088 - x347 * x555;
        auto x1090 = x1079 * x87 + x242 * x564 - x242 * x676 - x556 * x558;
        auto x1091 =
            16.0 * x1073 + x1078 - x1080 + x1087 + x1089 + x1090 + x136 * x968;
        auto x1092 = x142 * x83; auto x1093 = x134 * x144;
        auto x1094 = x1031 * x211; auto x1095 = t * x79; auto x1096 = MUs * x83;
        auto x1097 = x1096 * x170; auto x1098 = 6.0 * x1044;
        auto x1099 = AXG * x1092 + AXG * x1093 + x1051 + x1056 + x1061 + x1081 +
                     x1083 + x1089 - x1092 - x1093 + x1094 + x1095 * x217 -
                     x1096 * x183 + x1097 - x1098 + x136 * x307 + x160 -
                     x402 * x589 - x426 + x528 + x689 * x88 - x770;
        auto x1100 = s * x518; auto x1101 = 2.0 * x393;
        auto x1102 = x185 * x874; auto x1103 = MXs * x1102;
        auto x1104 = x137 * x478; auto x1105 = x1072 + x1104;
        auto x1106 = -x1071 - x561; auto x1107 = x1031 * x218;
        auto x1108 = x1103 + x1107;
        auto x1109 = -AXG * x1103 - x1011 * x151 + 6.0 * x1068 - 6.0 * x1070 +
                     x1084 + x1086 + x1090 + x1096 * x145 - x1101 + x1105 +
                     x1106 + x1108 - x132 * x185 - x134 * x155 + x151 * x451 +
                     x399 * x548 - x49 * x79 + x51 * x558 + x578 * x95;
        auto x1110 = 16.0 * x116; auto x1111 = 14.0 * x1064;
        auto x1112 = -4.0 * x393; auto x1113 = MXs * x87;
        auto x1114 = 12.0 * x116; auto x1115 = x263 * x537;
        auto x1116 = x116 * x522; auto x1117 = s * x543; auto x1118 = s * x599;
        auto x1119 = x158 * x83; auto x1120 = AXG * x1119;
        auto x1121 = 3.0 * x4; auto x1122 = x170 * x30;
        auto x1123 = MUs * x1122; auto x1124 = x183 * x30;
        auto x1125 = 3.0 * x1124; auto x1126 = MUs * x1125; auto x1127 = -x1126;
        auto x1128 = MUs * x489; auto x1129 = 13.0 * x434;
        auto x1130 = x242 * x293;
        auto x1131 =
            AXG * x1128 + AXG * x1129 + AXG * x1130 + x1119 - x1128 - x1130;
        auto x1132 = x10 * x1041; auto x1133 = 2.0 * x1132;
        auto x1134 = x1132 * x206; auto x1135 = -x1133 + x1134;
        auto x1136 = s * x610; auto x1137 = x367 * x874;
        auto x1138 = 7.0 * x1076; auto x1139 = x1121 * x218 - 20.0 * x977;
        auto x1140 = Mqi * x1040; auto x1141 = MQis * x431;
        auto x1142 = Mqi * x24; auto x1143 = x134 * x427;
        auto x1144 = x224 * x44; auto x1145 = x242 * x998;
        auto x1146 = x271 * x44;
        auto x1147 = -AXG * x1144 + AXG * x1145 + AXG * x1146 - 10.0 * x1143 +
                     x1144 - x1145;
        auto x1148 = s * x626; auto x1149 = 4.0 * x384;
        auto x1150 = MUs * x1149; auto x1151 = x558 * x851;
        auto x1152 = Mqi * x136; auto x1153 = x1152 * x24;
        auto x1154 = x151 * x530; auto x1155 = x550 * x552;
        auto x1156 = MXs * x250; auto x1157 = MQis * x1156;
        auto x1158 = x222 * x4; auto x1159 = 7.0 * MXs;
        auto x1160 = x1158 * x1159; auto x1161 = AXG * x93;
        auto x1162 = x1161 * x211; auto x1163 = 4.0 * x224;
        auto x1164 = x1163 * x150; auto x1165 = 6.0 * x242;
        auto x1166 = x1165 * x344; auto x1167 = AXG * MXs;
        auto x1168 = 11.0 * x1167; auto x1169 = MUs * x87;
        auto x1170 = x1169 * x384; auto x1171 = x1142 * x145;
        auto x1172 = x1167 * x226; auto x1173 = AXG * x1166;
        auto x1174 = x174 * x222; auto x1175 = AXG * x1174;
        auto x1176 = MXs * x36; auto x1177 = x1176 * x222;
        auto x1178 = AXG * x1177; auto x1179 = -x1174 + x1175 + x1178;
        auto x1180 = t * x463; auto x1181 = 5.0 * x1180;
        auto x1182 = x284 * x465; auto x1183 = x191 * x369;
        auto x1184 = t * x325; auto x1185 = 3.0 * x1184;
        auto x1186 = x1121 * x222; auto x1187 = t * x1186;
        auto x1188 = 5.0 * x30; auto x1189 = x1188 * x170;
        auto x1190 = MUs * x1189; auto x1191 = MQis * x1122;
        auto x1192 = 3.0 * x1191; auto x1193 = x568 * x874;
        auto x1194 = t * x1193; auto x1195 = 10.0 * x478;
        auto x1196 = x1195 * x170; auto x1197 = t * x474;
        auto x1198 = MQis * x1125; auto x1199 = AXG * x1137;
        auto x1200 = t * x1199; auto x1201 = x1121 * x224;
        auto x1202 = t * x1201; auto x1203 = 14.0 * x183;
        auto x1204 = x1203 * x478; auto x1205 = t * x384;
        auto x1206 = x1205 * x206; auto x1207 = AXG * x1187;
        auto x1208 = MX * x4; auto x1209 = Mqi * x1208;
        auto x1210 = x1209 * x316; auto x1211 = 3.0 * x1210;
        auto x1212 = x1142 * x170; auto x1213 = 3.0 * x1212;
        auto x1214 = x402 * x598; auto x1215 = x1208 * x194;
        auto x1216 = x1215 * x206; auto x1217 = t * x1216;
        auto x1218 = x1142 * x183; auto x1219 = 3.0 * x1218;
        auto x1220 = x1210 * x895; auto x1221 = x1201 * x402;
        auto x1222 = AXG * x1101; auto x1223 = MQis * x991;
        auto x1224 = MXs * x1223; auto x1225 = x1 * x174;
        auto x1226 = x1225 * x242; auto x1227 = x102 * x4;
        auto x1228 = Mqi * x1227; auto x1229 = MXs * x1228;
        auto x1230 = AXG * x1225; auto x1231 = x1230 * x242;
        auto x1232 = x174 * x224; auto x1233 = AXG * x1232;
        auto x1234 = x1176 * x224; auto x1235 = AXG * x1234;
        auto x1236 = -x1232 + x1233 + x1235; auto x1237 = -x1177 - x1234;
        auto x1238 = s * x684; auto x1239 = x10 * x103;
        auto x1240 = Mqi * x1239; auto x1241 = MQis * x489;
        auto x1242 = s * x718; auto x1243 = MQis * x373;
        auto x1244 = x143 * x658; auto x1245 = AXG * x1243;
        auto x1246 = x271 * x347; auto x1247 = AXG * x453;
        auto x1248 = x550 * x874; auto x1249 = MXs * x116 * x1163;
        auto x1250 = x555 * x874; auto x1251 = x224 * x4;
        auto x1252 = x1159 * x1251; auto x1253 = AXG * x174 * x271;
        auto x1254 = x134 * x174; auto x1255 = x4 * x415;
        auto x1256 = 17.0 * x1255; auto x1257 = AXG * x1256;
        auto x1258 = x1 * x1209; auto x1259 = 9.0 * x1258;
        auto x1260 = AXG * x1215; auto x1261 = AXG * x373;
        auto x1262 = x36 * x415; auto x1263 = AXG * x1254;
        auto x1264 = AXG * x1262; auto x1265 = AXG * x272;
        auto x1266 = x556 * x722;
        auto x1267 = x1264 + x1265 + x1266 - x723 + x780;
        auto x1268 = -x1189 - x1261 - x1262 + x1263 + x1267 + x373;
        auto x1269 = AXG * x1223; auto x1270 = 3.0 * x360;
        auto x1271 = Mqi * x103 * x1167; auto x1272 = AXG * x274;
        auto x1273 = -AXG * x1228 + x1228; auto x1274 = 11.0 * x1255;
        auto x1275 = AXG * x1102; auto x1276 = AXG * x1274;
        auto x1277 = x1209 * x511; auto x1278 = AXG * x1277;
        auto x1279 = x1208 * x526; auto x1280 = -x1215 * x87 + x1279;
        auto x1281 = x718 * x8; auto x1282 = 4.0 * x853; auto x1283 = MUs * x89;
        auto x1284 = x1045 * x743; auto x1285 = x1040 * x264;
        auto x1286 = x242 * x602; auto x1287 = x264 * x88;
        auto x1288 = x286 * x88; auto x1289 = 15.0 * AXG;
        auto x1290 = x10 * x299; auto x1291 = 6.0 * x1290;
        auto x1292 = x1 * x427; auto x1293 = x1292 * x522;
        auto x1294 = AXG * x1291; auto x1295 = x37 * x774;
        auto x1296 = x1060 * x574; auto x1297 = x235 * x887;
        auto x1298 = x100 * x347; auto x1299 = x551 - x557;
        auto x1300 = u * x993 - x1291 - 12.0 * x1292 + x1293 + x1294 + x1295 +
                     x1296 - x1297 - x1298 + x1299 - x854 - x855;
        auto x1301 = s * x759; auto x1302 = 3.0 * x863; auto x1303 = x4 * x753;
        auto x1304 = x24 * x48; auto x1305 = x766 * x874;
        auto x1306 = 7.0 * x478; auto x1307 = AXG * x1304;
        auto x1308 = AXG * x778; auto x1309 = 11.0 * x1308;
        auto x1310 = x145 * x802; auto x1311 = x767 * x874;
        auto x1312 = x217 * x478; auto x1313 = x191 * x766;
        auto x1314 = AXG * x1313; auto x1315 = x1 * x151;
        auto x1316 = x1315 * x710; auto x1317 = -x1316;
        auto x1318 = x1314 + x1317 + 11.0 * x864; auto x1319 = AXG * x752;
        auto x1320 = AXG * x848; auto x1321 = 6.0 * x869;
        auto x1322 = -x1007 * x242 - x412 * x558 + x869 * x87 + x873;
        auto x1323 = x1000 * x694 + x1014 * x116 - x1065 * x51 + 16.0 * x1315 -
                     x1319 - x1320 - x1321 + x1322 + x752;
        auto x1324 = x103 * x48; auto x1325 = MUs * x79;
        auto x1326 = 2.0 * x1305; auto x1327 = x217 * x991;
        auto x1328 = -x1313 - x870;
        auto x1329 =
            s * (-MUs * x827 + 6.0 * x1307 + 6.0 * x1308 - 2.0 * x1311 + x1314 -
                 x1315 * x522 + x1322 - x1324 + x1325 * x145 + x1326 + x1327 +
                 x1328 + x151 * x851 + x412 * x562 - x49 * x83 - x738 + x739 -
                 2.0 * x863 + 6.0 * x864 + x867 - x871);
        auto x1330 = 6.0 * x853; auto x1331 = 12.0 * x1290;
        auto x1332 = 32.0 * AXG; auto x1333 = x37 * x531;
        auto x1334 = 6.0 * x10; auto x1335 = x1334 * x264;
        auto x1336 = x1334 * x286; auto x1337 = x347 * x37;
        auto x1338 = 26.0 * AXG; auto x1339 = s * x861;
        auto x1340 = -4.0 * x1304; auto x1341 = x552 * x766;
        auto x1342 = 2.0 * x729; auto x1343 = x10 * x1325;
        auto x1344 = AXG * x1343; auto x1345 = 13.0 * x743;
        auto x1346 = AXG * x1345; auto x1347 = MX * x1239;
        auto x1348 = x242 * x489; auto x1349 = AXG * x1347;
        auto x1350 = x242 * x494; auto x1351 = AXG * x859;
        auto x1352 = s * x893; auto x1353 = x1005 * x4; auto x1354 = 9.0 * x864;
        auto x1355 = AXG * x860; auto x1356 = -x1304 - x863;
        auto x1357 = MXs * x45; auto x1358 = x1 * x304; auto x1359 = x142 * x79;
        auto x1360 = x1176 * x899; auto x1361 = x102 * x1208;
        auto x1362 = MQis * x1361; auto x1363 = x1254 * x242;
        auto x1364 = MXs * Mqi; auto x1365 = x1031 * x1364 * x194;
        auto x1366 = x277 * x874; auto x1367 = Mqis * x1361;
        auto x1368 = x1208 * x23; auto x1369 = 4.0 * x1368;
        auto x1370 = x1263 * x242; auto x1371 = x1209 * x185;
        auto x1372 = s * x1024; auto x1373 = 4.0 * x775;
        auto x1374 = MUs * x1373; auto x1375 = x1152 * x360;
        auto x1376 = x151 * x912; auto x1377 = x174 * x898;
        auto x1378 = MXs * x51; auto x1379 = x116 * x1378;
        auto x1380 = MQis * x1379; auto x1381 = x552 * x98;
        auto x1382 = x4 * x766; auto x1383 = 7.0 * x1382;
        auto x1384 = MQis * x1383; auto x1385 = AXG * x1377;
        auto x1386 = x36 * x766; auto x1387 = AXG * x1386;
        auto x1388 = MQis * x1387; auto x1389 = x1161 * x235;
        auto x1390 = Mqi * x150 * x526; auto x1391 = 11.0 * x1382;
        auto x1392 = AXG * x1391; auto x1393 = x1169 * x775;
        auto x1394 = x145 * x931; auto x1395 = x1023 * x1167;
        auto x1396 = s * x908; auto x1397 = 5.0 * x785; auto x1398 = t * x777;
        auto x1399 = x1121 * x898; auto x1400 = t * x1399;
        auto x1401 = x11 * x170; auto x1402 = 3.0 * MQis;
        auto x1403 = x874 * x877; auto x1404 = AXG * x1401;
        auto x1405 = 3.0 * x875; auto x1406 = AXG * x1405;
        auto x1407 = x210 * x319; auto x1408 = AXG * x1407;
        auto x1409 = Mqi * x1270; auto x1410 = 9.0 * x792;
        auto x1411 = AXG * x1410;
        auto x1412 = -AXG * x1400 + AXG * x37 * x44 + 5.0 * MUs * x1401 -
                     3.0 * MUs * x1404 + t * x1397 - t * x1403 + t * x1406 -
                     t * x1411 + t * x793 - t * x799 + x1203 * x802 + x1299 -
                     x1398 + x1400 - x1401 * x1402 + x1402 * x1404 - x1408 -
                     x1409 * x170 + x1409 * x183 - 10.0 * x170 * x802 + x401;
        auto x1413 = s * x947; auto x1414 = MQis * x264;
        auto x1415 = MUs * x264; auto x1416 = Mqis * x264;
        auto x1417 = x242 * x469; auto x1418 = MXs * x83;
        auto x1419 = 2.0 * x360; auto x1420 = x242 * x471;
        auto x1421 = 3.0 * x264; auto x1422 = AXG * x1420;
        auto x1423 = AXG * x1417; auto x1424 = x1031 * x898;
        auto x1425 = 2.0 * x875; auto x1426 = AXG * x1425;
        auto x1427 = x1031 * x899; auto x1428 = AXG * x1424;
        auto x1429 = x1209 * x399; auto x1430 = AXG * MUs * x240;
        auto x1431 = x412 * x722; auto x1432 = AXG * x1429;
        auto x1433 = -x1368 * x87 + x1369; auto x1434 = MUs * x588;
        auto x1435 = x425 * x530; auto x1436 = MX * x961;
        auto x1437 = 6.0 * x1132; auto x1438 = x1165 * x12;
        auto x1439 = AXG * x1143; auto x1440 = x1 * x534;
        auto x1441 = -x1436 * x347; auto x1442 = 2.0 * x1041;
        auto x1443 = x1041 * x206; auto x1444 = MXs * x1442;
        auto x1445 = MUs * MXs; auto x1446 = -x1226; auto x1447 = MXs * x1443;
        auto x1448 = s * x982; auto x1449 = x242 * x45;
        auto x1450 = MQis * x170; auto x1451 = 2.0 * x1215;
        auto x1452 = AXG * x1078; auto x1453 = -x1081 + x1232 - x1243;
        auto x1454 = 7.0 * x1258; auto x1455 = x242 * x661;
        auto x1456 = x1254 + x1272; auto x1457 = x1156 - x1263;
        auto x1458 =
            -AXG * x1455 - x1264 - x1265 + x1455 + x1456 + x1457 + x272 - x274;
        auto x1459 = 10.0 * x1255; auto x1460 = x1102 - x1275;
        auto x1461 = -x1266 - x415 * x887 + x47 * x555 + x723;
        auto x1462 = -x116 * x555 + x1189 + x1254 + x1261 + x1457 + x1461 -
                     x373 - x489 + x494 + x721;
        auto x1463 = x1188 * x1445; auto x1464 = 9.0 * x384;
        auto x1465 = 28.0 * AXG; auto x1466 = x242 * x47;
        auto x1467 = -x1466 * x51 + x242 * x300; auto x1468 = x154 * x415;
        auto x1469 = AXG * x1468;
        auto x1470 = x1209 * x51 - x1209 * x556 + x1280 + x1460;
        auto x1471 = x1167 * x154;
        auto x1472 = AXG * x403 + MQis * x397 + MUs * x1469 - MXs * x226 +
                     x1071 - x11 * x396 + x1141 + x1142 * x177 + x1174 + x1177 +
                     x1234 - x134 * x152 + x1446 + x1453 - x1471 * x222 -
                     x1471 * x224 + x390 - x392 - x398 - x442 + x443 * x530 -
                     x654;
        auto x1473 = MQis * x1418; auto x1474 = x415 * x93;
        auto x1475 = x1 * x144; auto x1476 = 6.0 * x1351;
        auto x1477 = x1227 * x210; auto x1478 = -x1353 - x1477;
        auto x1479 = s * (-AXG * x1359 - AXG * x1475 + x1 * x155 - x1295 -
                          x1296 + x1297 + x1298 + x1319 + x1320 + x1355 * x399 +
                          x1359 + x1475 + x1476 + x1478 - x418 + x419 + x732 -
                          x752 - x848 - x852 + x854 + x855 - x857 - x971);
        auto x1480 = AXG * x1477; auto x1481 = MQis * x45; auto x1482 = -x1481;
        auto x1483 = x631 * x746; auto x1484 = x44 * x899;
        auto x1485 = AXG * x1481 + AXG * x1484 - x1484;
        auto x1486 = -AXG * x1483 + AXG * x737 + x1291 + 10.0 * x1292 - x1293 -
                     x1294 + x1482 + x1483 + x1485;
        auto x1487 = MQis * x1031; auto x1488 = x217 * x471;
        auto x1489 = x174 * x899; auto x1490 = s * x945;
        auto x1491 = x217 * x307; auto x1492 = Mqi * x1419;
        auto x1493 = MQis * x293; auto x1494 = -x1375;
        auto x1495 = MQis * x1382; auto x1496 = x11 * x1445;
        auto x1497 = AXG * x1496; auto x1498 = x1306 * x242;
        auto x1499 = -x1387; auto x1500 = 17.0 * x1382;
        auto x1501 = -AXG * x1500 + x1405 - x1406 + x1499;
        auto x1502 = AXG * x276; auto x1503 = -x1230;
        auto x1504 = x1225 + x1379 + x1503;
        auto x1505 = -x1417 + x1420 - x1422 + x1423 + x1502 + x1504 - x276;
        auto x1506 = 10.0 * x1382; auto x1507 = 6.0 * x1368;
        auto x1508 = x1425 - x1426; auto x1509 = -x1431;
        auto x1510 = x1509 + x249 + x251;
        auto x1511 = x258 * x87 - x259 - x87 * x920;
        auto x1512 = x1504 + x1508 + x1510 + x1511 - x248 * x87;
        auto x1513 = 5.0 * x1496; auto x1514 = 9.0 * x775;
        auto x1515 = 17.0 * x900; auto x1516 = x154 * x766;
        auto x1517 = x1209 * x238 - x1209 * x412;
        auto x1518 = 9.0 * x242 * x478; auto x1519 = x1121 * x899;
        auto x1520 = MQis * x785; auto x1521 = MUs * x785;
        auto x1522 = MQis * x344; auto x1523 = x217 * x689;
        auto x1524 = MQis * x920; auto x1525 = x1 * x150;
        auto x1526 = MQis * x1525; auto x1527 = x150 * x899;
        auto x1528 = -x1524 + x1526 + x1527; auto x1529 = x222 * x531;
        auto x1530 = x23 * x533; auto x1531 = SS * x305;
        auto x1532 = x264 * x304; auto x1533 = x150 * x746;
        auto x1534 = x134 * x722; auto x1535 = x1159 * x931;
        auto x1536 = x1466 * x399;
        auto x1537 = -AXG * x1000 * x143 + AXG * x1023 - AXG * x1536 -
                     AXG * x777 - SS * x1364 * x526 - x1023 +
                     x1167 * x35 * x899 + x1415 * x87 - 4.0 * x1415 +
                     x143 * x753 - x1502 + x1536 - x177 * x881 + x451 * x705 +
                     x777;
        auto x1538 = -x159; auto x1539 = -x164; auto x1540 = s * x174;
        auto x1541 = x177 * x67; auto x1542 = -x171; auto x1543 = -x636;
        auto x1544 = -x184;
        auto x1545 = -AXG * x1540 + x1538 + x1539 + x1540 + x1541 + x1542 +
                     x1543 + x1544 + x165 + x175 + x176 + x221 + x422 - x447 -
                     x638 + x639;
        auto x1546 = -x133 - x138 + x139 + x140 + x144 - x146 - x147 - x149 -
                     x152 + x153 + x155 + x156;
        auto x1547 = 3.0 * x117; auto x1548 = syFC3 * x135;
        auto x1549 = 2.0 * x106; auto x1550 = AXG - 1;
        ,
        TR * TR * gs * gs * (Nc - 1) * (Nc + 1) *
            (UU * x518 *
                 (AXG * x1033 + AXG * x1034 + AXG * x1035 - x1033 - x1034 -
                  x1035 + x1039) +
             UU * x599 *
                 (x1039 - x15 * x238 + x15 * x412 - x25 * x526 + x25 * x554 -
                  x51 * x84 + x556 * x84) +
             kL * syFC9 *
                 (AXG * x719 + s * x721 + s * x723 + x108 * x527 - x108 * x724 -
                  x112 * x238 + x121 * x238 + x123 * x412 - x231 * x263 +
                  x234 * x238 + x252 * x263 - x255 * x550 - x279 * x550 +
                  x279 * x555 - x47 * x586 + x513 + x515 - x620 - x719 + x720) +
             kL * x341 *
                 (AXG * x358 + AXG * x363 + AXG * x380 + x116 * x340 +
                  x166 * x369 + x353 - x358 - x363 - x377 - x380 - x575 - x577 +
                  x625) -
             syFC1 * x172 * (-x233 * x82 - x332 + x339 + x57 * x82) -
             syFC1 *
                 (AXG * x265 + AXG * x278 - u * x197 * x271 - x111 * x266 -
                  x117 * x2 * x36 + x117 * x266 + x117 * x270 + x117 * x273 +
                  x121 * x186 + x123 * x186 - x125 * x270 - x125 * x273 +
                  x135 * x267 - x135 * x268 - 2.0 * x168 + 2.0 * x169 - x188 -
                  x192 * x36 + x192 * x93 + x2 * x267 + x207 * x271 * x279 +
                  x212 * x272 + x229 - 2.0 * x244 + 2.0 * x246 - x247 * x276 +
                  x247 * x38 - x265 - x275 - x278) +
             syFC10 *
                 (AXG * x28 + AXG * x33 + x14 - x18 + x19 - x20 - x28 - x33 -
                  x34 + x37 * x42 * x43 - x39 - x41 * x42 + x7 + x75) -
             syFC11 * (-AXG * x282 - AXG * x283 + x282 + x289 + x292 + x302) -
             syFC14 * (-x17 * x87 - x27 * x87 + 4.0 * x27 + x302 + x303) +
             syFC15 * (x100 * x99 + x105 + x75 + x76 + x78 - x80 - x81 - x86 -
                       x90 + x92 - x97 - x98 * x99) -
             syFC2 * (-AXG * x308 - MQis * x275 + t * x229 +
                      x117 * x136 * x172 * x82 - x189 - x192 * x306 -
                      x201 * x304 - x208 - x225 * x313 + x227 * x311 + x228 +
                      x301 + x308 - x309 * x310 + x310 * x312 + x315 + x318 +
                      x39 + x80 + x86 - x92) -
             syFC3 * x2 * (x1545 + x1546) +
             syFC3 *
                 (x107 * x2 + x110 * x2 + x113 * x2 - x114 * x2 - x115 * x2 -
                  x119 * x2 - x122 * x2 - x124 * x2 + x126 * x2 + x127 * x2 +
                  x129 * x2 - x130 * x2 + x157 + x55 - x60 + x62 - x65) -
             syFC4 * (-x14 - x19 + x20 + x27 + x285 + x287 + x288 + x292 -
                      x314 + x321 + x34 + x56 - x7) -
             syFC5 * x173 *
                 (AXG * x1549 + AXG * x233 - 7.0 * u * x125 + u * x1547 + x120 -
                  x1549 + 5.0 * x209 - x233 - x329 - x477 + x57) -
             syFC5 * x2 *
                 (-AXG * x124 - x1021 + x1022 + x1538 + x1539 + x1541 + x1542 +
                  x1543 + x1544 + x1546 + x1547 * x4 - x161 + x165 + x175 +
                  x176 + x613 + x639) +
             syFC5 * (kL * x160 + x135 * x164 - x135 * x165 + x135 * x171 -
                      x135 * x175 - x135 * x176 + x135 * x184 + x148 * x173 +
                      x157 + x162 - x163 + x167 + x168 - x169 - x177 * x178 -
                      x178 * x88 + x180 - x182) +
             syFC7 * (AXG * x187 - MQis * x186 * x216 + kL * x219 * x220 +
                      x105 + x135 * x209 * x93 - x136 * x201 + x161 * x214 -
                      x187 + x189 + x191 * x192 - x195 * x197 - x196 * x199 -
                      x203 - x205 + x208 + x211 * x212 * x93 - x212 * x219 -
                      x214 * x216 + x215 - x221 * x223 - x224 * x225 + x228) -
             syFC8 * (-kL * x322 - x101 + x303 + x321 - x78 + x90) +
             syFC9 *
                 (AXG * x241 - kR * x245 * x47 + x100 * x260 - x112 * x230 +
                  x121 * x230 - x162 * x87 + 2.0 * x167 + x217 * x239 -
                  x217 * x254 + x229 + x230 * x234 + x231 * x232 - x232 * x252 -
                  x235 * x237 * x87 + x236 * x237 - x239 * x263 - x241 -
                  x242 * x256 - x243 * x83 + x243 * x91 - x244 * x87 +
                  x246 * x87 + x247 * x249 + x247 * x251 - x247 * x259 -
                  x247 * x262 + x254 * x263 - x260 * x98) -
             t * x626 * (x1472 + x219) - t * x718 * (x1472 + x426) +
             u * x1396 *
                 (SS * x236 + x1168 * x931 + 9.0 * x1496 - 11.0 * x1497 -
                  x1515 - x1535 + x1537 + 8.0 * x258 + 14.0 * x785 -
                  26.0 * x798 + x810) +
             u * x1413 *
                 (AXG * x1535 - x1270 * x1364 - 7.0 * x1497 + x1513 + x1537 +
                  10.0 * x258 + 8.0 * x785 + x793 - 20.0 * x798 - 13.0 * x900 +
                  x993) -
             x1002 * x759 - x1002 * x893 - x1009 * x759 - x1009 * x893 +
             x1024 * x8 *
                 (AXG * x1371 - AXG * x1414 + AXG * x1415 + AXG * x1416 +
                  MXs * x1325 - MXs * x1421 + x1020 - x1113 * x286 +
                  x1167 * x1421 - x1364 * x1419 - x1371 - x1378 * x4 + x1386 +
                  x1414 - x1415 - x1416 + x1417 + x1418 * x242 - x1420 + x1422 -
                  x1423 + x276 + x827) -
             x1024 *
                 (x1 * x111 * x306 - x1017 - x1020 * x111 - x1021 * x37 +
                  x1022 * x37 - x1023 * x111 + x107 * x37 + x107 * x746 -
                  x112 * x942 - x114 * x746 + x121 * x938 + x121 * x942 +
                  x185 * x985 + x267 * x898 + x267 * x899 - x268 * x898 + x823 -
                  x828 - x830 + x831 + x834 + x902 - x939 - x940 + x941) -
             x1025 * x908 - x1025 * x947 - x1026 * x908 - x1026 * x947 +
             x1063 * (-x1042 + x1043 - 15.0 * x1044 + x1053 + x1057 + x1058 +
                      x1062 - 11.0 * x434 + 7.0 * x539) +
             x1063 * (10.0 * x1064 - x1065 * x238 + x1066 + 7.0 * x1067 +
                      11.0 * x1068 - x1069 - 15.0 * x1070 + x1077 + x1091 -
                      3.0 * x393 + x548 * x888) +
             x1099 * x1100 + x1099 * x1118 + x1100 * x1109 + x1109 * x1118 -
             x1117 * (-AXG * x1434 - AXG * x1435 - AXG * x1437 - AXG * x1438 +
                      AXG * x411 + 26.0 * x1044 + x1055 * x522 - 12.0 * x1055 +
                      24.0 * x1143 + x1246 * x522 - x1332 * x434 + x1434 +
                      x1435 + x1436 * x531 + x1437 + x1438 - 16.0 * x1439 -
                      x1440 * x522 + x1441 - x271 * x563 - x411 + 16.0 * x434 +
                      x534 * x851 - 10.0 * x539) +
             x1117 * (-AXG * x1115 - AXG * x151 * x966 - MXs * x973 +
                      10.0 * x1067 + 20.0 * x1068 - 26.0 * x1070 +
                      32.0 * x1073 - 18.0 * x1076 + x1079 * x522 -
                      12.0 * x1079 + 20.0 * x11 * x95 + x1110 * x389 + x1111 +
                      x1112 + x1113 * x611 + x1114 * x218 + x1115 -
                      x1116 * x218 - x132 * x238 - x145 * x968 - x412 * x520 +
                      x521 + 10.0 * x542 + x548 * x966 - 46.0 * x977) +
             x1136 * (-11.0 * x1044 + x1052 + x1062 - x1120 + x1121 * x211 +
                      3.0 * x1123 + x1127 + x1131 + x1135 - 5.0 * x434) +
             x1136 * (MXs * x1137 + x1057 + 9.0 * x1068 - 11.0 * x1070 + x1075 +
                      x1091 + x1105 - x1138 + x1139 - x393) +
             x1148 * (MQis * x602 + x1053 - x1129 + x1133 - x1134 +
                      x1140 * x24 - x1141 - x1142 * x88 + x1147) +
             x1148 * (9.0 * x1067 + x1077 + x1101 + x1106 + x1150 + x1151 +
                      3.0 * x1153 + x1154 - x1155 - x1157 + x1158 * x1168 -
                      x1160 + x1162 + x1164 + x1166 - x1170 - 7.0 * x1171 -
                      x1172 - x1173 + x1179) -
             x1148 * (-AXG * x1094 - x1066 + x1069 + 17.0 * x1070 + x1082 +
                      x1085 - x1111 - x1168 * x1251 - x1233 - x1235 - x1244 +
                      x1245 + x1247 + x1248 + x1249 - x1250 + x1252 + x1253 +
                      x1441 + x1452 + x1453 + x4 * x527 - x528) +
             x1148 * (x1058 + x1127 + x1181 + x1182 - x1183 + x1185 + x1187 +
                      x1190 - x1192 - x1194 - x1196 - x1197 + x1198 + x1200 +
                      x1202 + x1204 - x1206 - x1207 - x1211 - x1213 - x1214 +
                      x1217 + x1219 + x1220 - x1221 + x781) +
             x1238 * (-AXG * x1107 + AXG * x1224 - AXG * x392 + x1104 + x1108 +
                      x1167 * x1228 + x1179 + x1222 - x1224 + x1226 - x1229 -
                      x1231 + x1236 + x1237 - x134 * x153 - x145 * x991 -
                      x174 * x415 - x174 * x416 + x392 + x442) -
             2.0 * x1238 *
                 (AXG * x1184 + AXG * x1205 - AXG * x1440 + AXG * x1529 +
                  AXG * x1530 - MQis * x1124 - t * x1158 + t * x611 + x1044 +
                  x1055 - x1123 - x1180 - x1184 + x1191 - x1205 + x1212 -
                  x1218 - x1439 + x1440 - x1529 - x1530 + x194 * x47 * x71 +
                  x434 - x539 + x559) +
             x1242 * (-AXG * x1240 - AXG * x1241 + x1042 - x1043 + x1046 +
                      x1048 + x1049 + x1050 + x1131 + x1147 + x1240 -
                      x142 * x355 - x271 * x631 + x519) +
             x1242 * (AXG * x1252 - MXs * x1201 + x1054 - 13.0 * x1070 + x1087 +
                      x1098 + x1139 + x1236 + x1241 + x1243 + x1244 - x1245 +
                      10.0 * x1246 - x1247 - x1248 - x1249 + x1250 - x1253 +
                      x154 * x389 - x457) -
             x1242 * (-AXG * x1160 + MXs * x1186 - MXs * x1193 + x1074 + x1080 +
                      x1112 + x1138 - x1150 - x1151 - x1153 - x1154 + x1155 +
                      x1157 - x1162 - x1164 - x1166 + x1170 + 5.0 * x1171 +
                      x1172 + x1173 + x1174 - x1175 - x1178 + x1222 + x1452) -
             x1242 * (t * x475 + x1120 + x1126 - x1181 - x1182 + x1183 - x1185 -
                      x1187 - x1190 + x1192 + x1194 + x1196 + x1197 - x1198 -
                      x1200 - x1202 - x1204 + x1206 + x1207 + x1211 + x1213 +
                      x1214 - x1217 - x1219 - x1220 + x1221 - x781) +
             x1281 * (-x1102 + x1125 + x1268 + x1273 - x1274 + x1275 + x1276 -
                      x1277 + x1278 + x1280 - x271 * x694 - x431 + x602) -
             x1281 *
                 (-AXG * x1473 + AXG * x1474 + MUs * x1442 - MUs * x1443 +
                  x103 * x1364 + x1149 - x1223 + x1269 - x1271 - x1445 * x204 +
                  x1445 * x91 + x1456 + x1467 + x1473 - x1474 - x206 * x384 -
                  x47 * x510 + x474 + 11.0 * x482 - x486 - x597) +
             x1301 * (-x1282 + x1283 + x1284 - x1285 - x1286 + x1287 + x1288 -
                      x1289 * x859 + x1300 - 11.0 * x743 + x857 + 7.0 * x859) +
             x1301 * (-x1289 * x1312 - x1302 + 10.0 * x1303 - 3.0 * x1304 +
                      7.0 * x1305 + x1306 * x217 + 11.0 * x1307 + x1309 -
                      26.0 * x1310 - 11.0 * x1311 + x1318 + x1323) +
             x1329 * x846 + x1329 * x892 +
             x1339 * (AXG * x1330 + AXG * x1331 + AXG * x1335 + AXG * x1336 -
                      AXG * x761 + x1011 * x534 + x1114 * x235 - x1116 * x235 +
                      x1292 * x540 - 24.0 * x1292 - x1330 - x1331 +
                      x1332 * x743 + x1333 * x522 - 8.0 * x1333 - x1335 -
                      x1336 - x1337 * x522 + 8.0 * x1337 - x1338 * x859 -
                      x451 * x534 - 16.0 * x743 + x761 + 10.0 * x859) +
             x1339 *
                 (-20.0 * AXG * x1315 - x1011 * x242 * x347 + x1110 * x753 -
                  x116 * x145 * x851 + x116 * x856 + x1195 * x217 +
                  14.0 * x1303 + 10.0 * x1305 + 20.0 * x1307 + 20.0 * x1308 -
                  46.0 * x1310 - 18.0 * x1311 - x1312 * x1338 + 32.0 * x1315 +
                  x1340 + x1341 * x87 - 8.0 * x1341 + x242 * x666 -
                  x520 * x556 + x522 * x869 + x847 - 4.0 * x863 + 20.0 * x864 -
                  12.0 * x869 - x968 * x978) +
             x135 * x351 *
                 (x112 - x121 + x128 + x166 + x179 - x181 + x342 - x343 + x40 -
                  x43 - x58) +
             x1352 *
                 (x1121 * x235 + x1300 + x1342 - x1343 + x1344 + x1346 - x1347 -
                  x1348 + x1349 + x1350 - 11.0 * x1351 - 5.0 * x743) +
             x1352 * (3.0 * x1305 + 9.0 * x1307 - 20.0 * x1310 - 7.0 * x1311 +
                      x1314 + x1317 + x1323 + x1353 + x1354 - x1355 * x888 +
                      x1356 + x367 * x860 + x857) -
             2.0 * x1372 *
                 (AXG * x1532 - AXG * x1533 - MQis * x258 + MXs * x1414 -
                  MXs * x1415 - MXs * x1416 + MXs * x785 + Mqis * x258 -
                  x1167 * x1414 + x1167 * x1415 + x1167 * x1416 - x1308 +
                  x1494 + x1495 + x150 * x37 - x1520 + x1521 + x1528 +
                  x1531 * x37 - x1531 * x766 - x1532 + x1533 - x299 * x304) +
             x1372 * (-AXG * x1362 + AXG * x1366 + AXG * x1367 - MQis * x38 -
                      MXs * x1369 + Mqis * x38 + x1031 * x1358 + x1113 * x1368 -
                      x1167 * x1371 + x1176 * x37 + x1262 * x242 -
                      x1264 * x242 + 2.0 * x1303 + x1326 + x1357 - x1358 * x36 -
                      x1359 + x1360 + x1362 + x1363 - x1365 - x1366 - x1367 -
                      x1370 - x151 * x51 + x857) +
             x1396 * (x1347 - x1349 + x1412 + 4.0 * x729 - x730) -
             x1396 * (-MQis * x77 + MQis * x89 + x1027 - x1140 * x360 + x1282 -
                      x1283 - x1284 + x1286 - x1288 + x1345 - x1480 + x1486 +
                      x236 * x4 + x854) +
             x1396 * (MQis * x1392 + MUs * x799 + 13.0 * x1308 + x1318 + x1324 +
                      x1328 + x1374 + 3.0 * x1375 + x1376 - x1377 - x1380 -
                      x1381 - x1384 + x1385 + x1388 + x1389 + x1390 - x1393 -
                      7.0 * x1394 - x1395 - 3.0 * x778 - 5.0 * x863) +
             x1413 * (-AXG * x1342 + x1285 - x1287 + x1342 + x1412) -
             x1413 * (AXG * x1493 - x10 * x1492 + x10 * x206 * x931 + x1343 -
                      x1344 - x1346 + x1348 - x1350 - x1476 + x1486 - x1493 +
                      x37 * x631 + 7.0 * x743 + x751) -
             x1413 *
                 (AXG * x1324 - AXG * x1384 + x1302 - x1309 + x1313 - x1314 +
                  x1316 + x1321 + x1340 - x1354 - x1374 - x1376 + x1377 +
                  x1380 + x1381 - x1385 - x1388 - x1389 - x1390 + x1393 +
                  5.0 * x1394 + x1395 + x1494 + 3.0 * x1495 + x778) -
             x1448 * (MQis * x1444 - MQis * x1447 + MXs * x274 - Mqis * x1444 +
                      Mqis * x1447 + Mqis * x274 - x103 * x1152 - x1103 +
                      2.0 * x1209 * x767 + x1224 + x1231 + x1237 - x136 * x991 -
                      x1386 * x242 + x1387 * x242 - x1442 * x1445 -
                      x1442 * x304 + x1443 * x1445 + x1443 * x304 + x1446 +
                      x151 * x238 + x271 * x306 + x304 * x471 - x304 * x991) -
             x1448 * (-AXG * x1449 + MQis * x1216 - MQis * x1451 + MQis * x272 +
                      MXs * x1279 - MXs * x373 - Mqis * x1216 + Mqis * x1451 -
                      Mqis * x272 - t * x983 + x1059 - x1088 + x1092 - x1097 -
                      x1113 * x1215 + x1135 + x1146 - x1176 * x271 + x1229 +
                      x1449 + x1450 * x83 - x206 * x271 * x874 - x222 * x386 +
                      x306 * x369 + x502 * x874) -
             x1479 * x846 - x1479 * x892 +
             2.0 * x1490 *
                 (AXG * x1522 + AXG * x1523 + AXG * x1524 - AXG * x1526 -
                  AXG * x1527 - MQis * x798 + x1015 + x1307 + x1308 + x1341 +
                  x1356 + x1375 - x1394 + x1520 - x1521 - x1522 - x1523 +
                  x1528 - x735 - x778 + x864 + x869) -
             x1490 * (AXG * x1327 - AXG * x1360 + AXG * x1488 - AXG * x1489 +
                      x1 * x153 + x1032 * x145 - x1167 * x1427 - x1326 - x1327 +
                      x1360 - x1363 + x1365 + x1370 + x1478 + x1480 + x1482 +
                      x1487 * x766 - x1487 * x767 - x1488 + x1489 +
                      x174 * x766 + x174 * x767 + x737 - x754) -
             x1490 *
                 (-AXG * x1357 + AXG * x1491 - MQis * x183 * x79 - t * x1424 +
                  t * x1425 + t * x1428 - t * x827 - x1032 * x183 +
                  x1095 * x145 - x1261 * x242 - x1325 * x170 + x1359 + x1398 -
                  x1407 + x1408 + x1450 * x79 + x1485 - x1491 + x1492 * x170 -
                  x1492 * x183 + x170 * x833 + x242 * x373 + x37 * x386 -
                  x565) -
             x1548 * x1550 * x309 -
             x1548 * (-x107 - x110 - x113 + x114 + x115 + x119 + x122 + x124 -
                      x126 - x127 - x129 + x130 + x1545 + x59 - x61) -
             x1550 * x361 * x982 -
             x341 * x9 *
                 (AXG * x1013 + AXG * x1525 - AXG * x1534 - AXG * x248 +
                  AXG * x258 - AXG * x920 - x1013 + x1019 + x1401 - x1404 -
                  x1525 + x1534 + x248 - x258 + x68 + x72 - x74 - x775 - x785 -
                  x792 + x798 + x900 + x920) +
             x341 * (AXG * x17 + AXG * x297 + AXG * x32 - AXG * x324 -
                     AXG * x327 - AXG * x337 + AXG * x338 + kL * x334 -
                     kL * x340 * x47 - x116 * x328 + x116 * x330 -
                     x145 * x173 * x8 - x17 + x173 * x323 + x315 - x32 + x324 +
                     x327 + x328 * x47 - x330 * x47 + x337 - x338) -
             x351 * (-AXG * x345 - AXG * x348 - x116 * x125 * x2 - x116 * x350 -
                     x117 * x349 + x123 * x135 + x125 * x349 - x128 * x2 +
                     x135 * x216 + x135 * x346 - x162 + x163 - x167 - x168 +
                     x169 - x180 + x182 - x192 * x47 + x2 * x216 - x2 * x342 +
                     x2 * x343 + x2 * x346 - x244 + x246 + x345 + x348 +
                     x350 * x47 - x54 + x64) -
             x383 * x8 *
                 (-AXG * x1454 + AXG * x311 + x1137 - x1199 + 13.0 * x1255 -
                  x1257 + x1454 + x1458 + x325 * x710 + x384 * x710 +
                  10.0 * x488 - x508 - x587 - x589) +
             x383 * (-x385 + x421) - x383 * (x949 + x953) -
             x383 * (x957 + x96) + x383 * (x441 + x445 + x449 + x458) +
             x383 *
                 (x122 * x404 - x127 * x404 + x404 * x59 + x460 - x461 + x485) +
             x383 *
                 (-AXG * x357 + AXG * x362 + AXG * x364 + x354 + x357 - x359 -
                  x362 - x364 - 7.0 * x365 + 7.0 * x366 + x372 + x382) -
             x518 * x8 *
                 (-AXG * x1459 + 6.0 * x1215 - 6.0 * x1260 + x1277 - x1278 +
                  x1459 + x1460 + x1462 - x480 + 11.0 * x593) +
             x518 *
                 (7.0 * x326 - x482 * x487 + x486 * x8 - x487 * x488 + x517) -
             x543 * x8 *
                 (-AXG * x1193 + AXG * x1463 - x1007 + x1193 - x1255 * x1465 +
                  20.0 * x1255 + 11.0 * x1258 + x1461 - x1463 - x1464 + x1467 -
                  9.0 * x463 - x47 * x550 + 17.0 * x482 + 17.0 * x488 +
                  17.0 * x593 - x598 + x872 + x969) -
             x543 *
                 (-AXG * x190 * x969 + x12 * x217 * x522 - x195 * x970 -
                  x211 * x963 + x218 * x563 - x218 * x963 - x238 * x958 -
                  x242 * x732 + x242 * x971 + x271 * x47 * x88 - x389 * x959 -
                  x425 * x550 + x434 * x965 + x47 * x966 * x967 - x522 * x948 +
                  x530 * x962 - x554 * x960 + x611 * x88 + 24.0 * x634 -
                  x689 * x964 + 12.0 * x948 + x960 * x961) +
             x543 *
                 (-AXG * x569 - AXG * x570 + AXG * x572 + AXG * x576 +
                  AXG * x580 + x166 * x571 + x238 * x567 + x313 * x581 +
                  6.0 * x353 - 11.0 * x365 + 11.0 * x366 - 6.0 * x377 +
                  x40 * x571 - x43 * x571 + x569 + x570 - x573 * x574 -
                  6.0 * x575 - x576 + x577 * x87 - 4.0 * x577 - x580 + x583) +
             x543 *
                 (MUs * x216 * x530 + MUs * x465 * x529 + t * x519 + t * x521 -
                  t * x523 + t * x528 + x12 * x242 * x529 + x138 * x369 -
                  x202 * x238 + x209 * x451 * x47 - x231 * x49 + x347 * x532 -
                  x425 * x524 - x428 * x451 + 24.0 * x435 - 24.0 * x444 -
                  12.0 * x455 + x456 * x536 - 16.0 * x525 - x531 * x532 -
                  x534 * x535 + x537 * x538 - x539 * x541 - x541 * x542) -
             x543 *
                 (MUs * x451 * x548 - x136 * x973 + x177 * x611 + x217 * x868 +
                  x300 * x49 - 16.0 * x347 * x632 + x360 * x974 + x378 * x974 -
                  x412 * x909 + x412 * x972 + x451 * x909 - x451 * x972 -
                  24.0 * x478 * x648 - x49 * x537 + x522 * x645 + x554 * x975 -
                  x556 * x656 - x566 * x977 + x656 * x851 + 28.0 * x659 +
                  12.0 * x702 - x802 * x976 - x923 * x978 - x961 * x975) +
             x543 * (AXG * x546 + MXs * x523 - t * x561 - x134 * x545 -
                     x140 * x369 - x183 * x478 * x566 - x211 * x563 +
                     x238 * x544 - x238 * x553 - x242 * x551 + x242 * x557 +
                     x242 * x565 - x384 * x536 + x394 * x522 +
                     x395 * x49 * x522 + x400 * x530 - x409 * x550 +
                     x409 * x555 - x415 * x547 + x416 * x549 + x524 * x548 +
                     x535 * x558 - x535 * x562 - x546 + x559 * x560) +
             x543 *
                 (-AXG * x585 + x111 * x587 + x111 * x588 + x111 * x589 -
                  x112 * x510 - x120 * x250 + x120 * x261 + x121 * x590 -
                  x179 * x590 + x181 * x590 + x271 * x595 + x313 * x595 -
                  x399 * x495 + x451 * x459 - 14.0 * x465 * x592 - x473 * x530 +
                  28.0 * x479 - x482 * x594 + x497 - 14.0 * x499 - x505 * x511 +
                  x509 + x531 * x586 + x585 + x591 - x593 * x594) -
             x599 * x8 *
                 (AXG * x1464 + x1462 + x1468 - x1469 + x1470 - 5.0 * x384) +
             x599 * (-AXG * x598 * x8 - 9.0 * x488 * x8 + x517 + 5.0 * x596 +
                     x597 * x8) -
             x610 * x8 *
                 (7.0 * x1255 - x1276 + x1458 + x1470 - x464 - x468 + x474 +
                  x481 + x483) -
             x610 * (x385 + x953) + x610 * (x421 + x445) -
             x610 * (x949 + x957) + x610 * (3.0 * x326 + x485 + x516 - x615) +
             x610 * (-x118 * x611 - x134 * x202 + x440 - x446 + 9.0 * x448 +
                     x458 + x614) +
             x610 * (x371 + x382 - x431 * x8 + x601 + x602 * x8 - x603 - x604 -
                     x605 + x606 + x609) +
             x626 * x8 *
                 (-AXG * x1201 + AXG * x1259 - x1137 + x1199 + x1201 + x1215 -
                  x1254 - x1256 + x1257 - x1259 - x1260 + x1268 + x356 - x475 -
                  x721) +
             x626 * (x628 + x637) +
             x626 * (-x202 * x568 - x222 * x422 + x222 * x447 + x441 -
                     5.0 * x446 + 13.0 * x448 + 11.0 * x612 + x640) +
             x626 * (x372 + x583 + x618 + x619 - x621 - x622 + x623 + x624 +
                     x625) -
             x626 * (x692 + x693 - x697 - x699 - x701 - x704 - x713 + x715 +
                     x716 + x717 + x952 + x955) -
             x626 * (AXG * x987 + MUs * x620 - MUs * x720 - x222 * x309 +
                     x222 * x312 + x222 * x59 - x222 * x61 + x224 * x59 -
                     x224 * x61 + x326 * x988 - 8.0 * x326 + x340 * x93 + x461 +
                     x476 + x488 * x896 - x491 - x568 * x897 - 10.0 * x596 -
                     x984 + x987 + x989) +
             x626 * (x322 + x642 + x643 - x646 + x647 - x649 - x650 - x653 -
                     x655 - x657 - x660 + x662 + x663 + x665 + x667 + x668 +
                     x669 - x670 - x671 - x672 - x674 + x675 + x677 - x678) +
             x684 * x8 *
                 (-AXG * x468 + x1167 * x991 - x1209 * x504 + x1223 -
                  6.0 * x1255 + x1258 * x206 + x1267 - x1269 - x1270 * x71 +
                  x1271 - x1272 + x1273 - x170 * x204 + x183 * x83 - x272 +
                  x274 + x438 * x87 - 4.0 * x438 + x468) +
             x684 * (x609 + x618 - x679 + x681 - x682 + x683) -
             x684 * (AXG * x990 + MQis * x500 - MQis * x515 + x181 * x213 -
                     x181 * x502 - x267 * x271 + x30 * x937 - x477 * x991 +
                     x503 + x506 - x587 * x8 + x986 + x989 - x990) +
             x718 * (x387 + x637) -
             x718 * (-x647 - x663 - x669 + x670 + x671 + x672 - x675 + x678 +
                     x954 + x956) -
             x718 * (x370 - x601 + x604 - x607 + x617 - x619 + x621 + x622 -
                     x623 - x624 + x679 + x682 - x683) +
             x718 * (-x117 * x391 + x125 * x391 + x161 * x222 - x222 * x613 +
                     x423 + x424 + x429 + x430 + x433 + x436 + x437 + x439 +
                     x449 + x614 + x640) +
             x718 * (x628 + x685 + x686 + x688 + x691 - x692 - x693 - x696 +
                     x697 + x698 + x699 + x701 + x703 + x704 - x707 - x708 -
                     x709 - x712 + x713 - x714 - x715 - x716 - x717) +
             x759 * x783 -
             x759 * x8 *
                 (-AXG * x1498 + 10.0 * x1019 + 13.0 * x1382 + 3.0 * x1497 +
                  x1498 + x1501 + x1505 + 10.0 * x798 - x836 - x837 +
                  10.0 * x900 - x933) +
             x759 * (x422 * x725 - x726 - 3.0 * x727 + x758) +
             x759 *
                 (x158 * x231 - x198 * x763 + x59 * x725 + x784 + x787 + x788 +
                  x789 + x791 + x794 - x795 - x796 - x797 - x800 + x801 + x804 -
                  x805 - x806 - x809 - x812 + x813 + x814 + x819 + x823) +
             x783 * x893 -
             x8 * x846 *
                 (-AXG * x1506 - AXG * x1507 + 11.0 * x1019 + x1429 - x1432 +
                  x1506 + x1507 + x1512 + 11.0 * x798 - x807 + 11.0 * x900 -
                  x946) -
             x8 * x861 *
                 (-AXG * x1403 + AXG * x1513 + 17.0 * x1019 + x1209 * x888 -
                  x1382 * x1465 + 20.0 * x1382 + x1403 - x1410 + x1466 * x412 +
                  x1509 + x1511 - x1513 - x1514 + x1515 + x251 - x262 + x564 -
                  x676 - 9.0 * x785 + 17.0 * x798 + 8.0 * x920) -
             x8 * x892 *
                 (AXG * x1514 - AXG * x1516 - x1397 + x1411 + x1433 + x1512 +
                  x1516 + x1517 - 5.0 * x775 + 9.0 * x798) -
             x8 * x893 *
                 (x1383 - x1392 + x1499 + x1505 + x1508 + x1517 + x38 - x786 -
                  x790 - x793 + x799 + x808 + x811) -
             x8 * x908 *
                 (AXG * x1368 + AXG * x1399 - AXG * x1518 + AXG * x1519 +
                  x1019 + x1225 - x1368 + x1386 - x1399 - x1430 + x1500 +
                  x1501 + x1503 + x1510 + x1518 - x1519 - x45 - x73 + x775) +
             x8 * x947 *
                 (-AXG * x1427 + AXG * x829 - x1225 + x1230 - x1373 - x1386 +
                  x1387 - x1391 + x1392 + x1424 - x1425 + x1426 + x1427 -
                  x1428 - x1429 + x1430 + x1431 + x1432 + x1433 - x251 -
                  x37 * x694 + x45 + x73) +
             x846 * (x8 * x810 + x826 + x845) -
             x861 * (x1010 * x88 + x1012 * x531 + x1013 * x88 -
                     x23 * x48 * x970 - x235 * x963 - x242 * x426 + x258 * x88 -
                     x299 * x964 + 20.0 * x299 * x967 - x425 * x98 -
                     x51 * x958 - x522 * x921 - x522 * x994 - x531 * x856 +
                     x536 * x792 + x743 * x965 - x753 * x959 - x87 * x992 +
                     x912 * x962 + 8.0 * x992 + 12.0 * x994 + 24.0 * x999) +
             x861 * (AXG * x858 - t * x10 * x537 + t * x551 + t * x847 +
                     t * x848 - t * x850 + t * x852 - t * x854 - t * x855 +
                     t * x857 - x108 * x849 - x202 * x51 + x209 * x537 +
                     12.0 * x216 * x725 - x255 * x856 - x369 * x540 * x860 -
                     x428 * x851 + x529 * x853 + x536 * x68 + x536 * x72 -
                     x541 * x859 - 16.0 * x736 + 24.0 * x744 - 24.0 * x755) -
             x861 * (MUs * x865 - 16.0 * u * x1015 - 24.0 * x1 * x687 +
                     x1003 * x87 - 8.0 * x1003 + x1004 * x522 + 28.0 * x1006 +
                     x1010 * x177 - x1012 * x347 - x1014 * x552 + x131 * x300 -
                     x131 * x537 - x145 * x876 + x24 * x974 + x31 * x974 +
                     x347 * x856 + x412 * x50 - x412 * x656 - x451 * x50 -
                     x478 * x976 + x49 * x968 - x556 * x909 + x851 * x909 +
                     12.0 * x921) +
             x861 * (AXG * x281 + AXG * x290 - AXG * x878 - AXG * x879 +
                     AXG * x880 + AXG * x885 - t * x886 - x116 * x23 * x883 +
                     x166 * x882 + x23 * x331 * x887 - x281 - x290 +
                     x37 * x581 + x375 * x511 + x402 * x886 - 11.0 * x41 -
                     x43 * x882 + x51 * x567 - x573 * x881 - x58 * x882 + x878 +
                     x879 - x880 - x885 + x888 * x889) +
             x861 * (AXG * x862 + MXs * x850 + t * x560 * x798 - t * x738 +
                     t * x865 + t * x867 - t * x870 - t * x871 + t * x873 -
                     x1 * x545 + x100 * x409 + x170 * x868 - x183 * x876 -
                     x235 * x563 - x242 * x519 + x242 * x866 - x409 * x98 +
                     x51 * x544 - x51 * x553 + x522 * x773 + x529 * x778 -
                     x536 * x775 - x547 * x766 + x549 * x767 - x862) -
             x861 *
                 (AXG * x1016 - x1016 + x1017 + x1018 * x112 - x1018 * x121 +
                  x1018 * x179 - x1018 * x181 + x1019 * x594 - x111 * x836 -
                  x111 * x837 + x111 * x901 - x111 * x933 + x112 * x838 +
                  x116 * x120 * x51 - x120 * x492 - x245 * x531 - x37 * x595 +
                  x399 * x505 + x473 * x912 - x537 * x57 - x595 * x746 -
                  28.0 * x803 + x832 + 14.0 * x835 + x903 * x988 - 6.0 * x903) +
             x892 * (x845 + x891) +
             x893 * (-x1 * x202 + x727 * x895 - x727 + x758 + x894) -
             x893 * (x158 * x198 * x206 - x784 - x787 - x788 - x789 - x791 -
                     x794 + x795 + x796 + x797 + x800 - x801 - x804 + x805 +
                     x806 + x809 + x812 - x813 - x814 + x815 + x817 + x818 -
                     x821 + x822 + x825 - x842 + x936) +
             x908 * x929 + x908 * x932 -
             x908 * (x1028 + x202 * x877 + x422 * x898 - x447 * x898) +
             x908 * (x257 * x883 - x59 * x898 - x59 * x899 - x61 * x725 +
                     x61 * x898 + x61 * x899 + x785 * x896 - x798 * x896 -
                     x8 * x901 + x877 * x897 + 8.0 * x890 + x907) +
             x929 * x947 + x932 * x947 +
             x945 * (-AXG * x935 - x11 * x937 - x181 * x938 + x267 * x37 +
                     x8 * x933 + x816 + x844 + x905 + x906 + x934 + x935 +
                     x936 + x939 + x940 + x944) -
             x945 * (-AXG * x1029 - AXG * x1030 + AXG * x934 + x1029 + x1030 -
                     x1032 * x329 - x116 * x839 + 4.0 * x16 - x185 * x889 +
                     x185 * x981 + x204 * x579 + x250 * x335 -
                     x31 * x331 * x779 + x316 * x616 - x316 * x943 - x317 * x8 -
                     x329 * x661 - x375 * x504 - x40 * x938 - x40 * x942 +
                     x43 * x938 + x43 * x942 + x511 * x979 + x535 * x58) -
             x947 * (x1028 - x161 * x898 + x613 * x898 + x726 - x894) +
             x947 * (x257 * x884 - x487 * x900 + x51 * x897 - x799 * x8 +
                     x8 * x946 + x891 + x907 + x944) -
             x982 * (MQis * x352 - MQis * x376 - MQis * x843 + x121 * x213 +
                     x222 * x267 + x224 * x267 - x224 * x268 +
                     x242 * x556 * x58 - x313 * x841 - x323 * x83 - x467 +
                     x472 + x504 * x985 + x615 - x8 * x983 + x984 + x986) -
             x982 * (-AXG * x980 - MUs * x352 + MUs * x376 - Mqis * x352 +
                     Mqis * x376 + t * x513 + t * x719 + x238 * x979 -
                     x271 * x309 + x271 * x312 + x271 * x59 - x271 * x61 +
                     x333 * x36 - x333 * x93 - x402 * x719 - x471 * x57 -
                     x504 * x889 + x504 * x981 + x603 + x605 - x606 + x617 -
                     x681 - x780 * x8 + x980)) *
            Denom(768.0 * Nc * (x111 - x117 + x8) * pow(Pi, 2)));
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