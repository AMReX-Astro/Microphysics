module actual_network

  use physical_constants, only: ERG_PER_MeV
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

  real(rt), parameter :: avo = 6.0221417930d23
  real(rt), parameter :: c_light = 2.99792458d10
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740d-12
  real(rt), parameter :: mev2erg = ev2erg*1.0d6
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184d-24
  real(rt), parameter :: mass_proton   = 1.67262163783d-24
  real(rt), parameter :: mass_electron = 9.10938215450d-28

  integer, parameter :: nrates = 506
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 189
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 189

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 506
  integer, parameter :: number_reaclib_sets = 762

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jp   = 1
  integer, parameter :: jd   = 2
  integer, parameter :: jhe3   = 3
  integer, parameter :: jhe4   = 4
  integer, parameter :: jli7   = 5
  integer, parameter :: jbe7   = 6
  integer, parameter :: jbe8   = 7
  integer, parameter :: jb8   = 8
  integer, parameter :: jc12   = 9
  integer, parameter :: jc13   = 10
  integer, parameter :: jn13   = 11
  integer, parameter :: jn14   = 12
  integer, parameter :: jn15   = 13
  integer, parameter :: jo14   = 14
  integer, parameter :: jo15   = 15
  integer, parameter :: jo16   = 16
  integer, parameter :: jo17   = 17
  integer, parameter :: jo18   = 18
  integer, parameter :: jf16   = 19
  integer, parameter :: jf17   = 20
  integer, parameter :: jf18   = 21
  integer, parameter :: jf19   = 22
  integer, parameter :: jf20   = 23
  integer, parameter :: jne18   = 24
  integer, parameter :: jne19   = 25
  integer, parameter :: jne20   = 26
  integer, parameter :: jne21   = 27
  integer, parameter :: jne22   = 28
  integer, parameter :: jna19   = 29
  integer, parameter :: jna20   = 30
  integer, parameter :: jna21   = 31
  integer, parameter :: jna22   = 32
  integer, parameter :: jna23   = 33
  integer, parameter :: jna24   = 34
  integer, parameter :: jmg21   = 35
  integer, parameter :: jmg22   = 36
  integer, parameter :: jmg23   = 37
  integer, parameter :: jmg24   = 38
  integer, parameter :: jmg25   = 39
  integer, parameter :: jmg26   = 40
  integer, parameter :: jal22   = 41
  integer, parameter :: jal23   = 42
  integer, parameter :: jal24   = 43
  integer, parameter :: jal25   = 44
  integer, parameter :: jal26   = 45
  integer, parameter :: jal27   = 46
  integer, parameter :: jal28   = 47
  integer, parameter :: jsi24   = 48
  integer, parameter :: jsi25   = 49
  integer, parameter :: jsi26   = 50
  integer, parameter :: jsi27   = 51
  integer, parameter :: jsi28   = 52
  integer, parameter :: jsi29   = 53
  integer, parameter :: jsi30   = 54
  integer, parameter :: jp26   = 55
  integer, parameter :: jp27   = 56
  integer, parameter :: jp28   = 57
  integer, parameter :: jp29   = 58
  integer, parameter :: jp30   = 59
  integer, parameter :: jp31   = 60
  integer, parameter :: jp32   = 61
  integer, parameter :: js28   = 62
  integer, parameter :: js29   = 63
  integer, parameter :: js30   = 64
  integer, parameter :: js31   = 65
  integer, parameter :: js32   = 66
  integer, parameter :: js33   = 67
  integer, parameter :: js34   = 68
  integer, parameter :: jcl29   = 69
  integer, parameter :: jcl30   = 70
  integer, parameter :: jcl31   = 71
  integer, parameter :: jcl32   = 72
  integer, parameter :: jcl33   = 73
  integer, parameter :: jcl34   = 74
  integer, parameter :: jcl35   = 75
  integer, parameter :: jcl36   = 76
  integer, parameter :: jcl37   = 77
  integer, parameter :: jar31   = 78
  integer, parameter :: jar32   = 79
  integer, parameter :: jar33   = 80
  integer, parameter :: jar34   = 81
  integer, parameter :: jar35   = 82
  integer, parameter :: jar36   = 83
  integer, parameter :: jar37   = 84
  integer, parameter :: jar38   = 85
  integer, parameter :: jk33   = 86
  integer, parameter :: jk34   = 87
  integer, parameter :: jk35   = 88
  integer, parameter :: jk36   = 89
  integer, parameter :: jk37   = 90
  integer, parameter :: jk38   = 91
  integer, parameter :: jk39   = 92
  integer, parameter :: jk40   = 93
  integer, parameter :: jk41   = 94
  integer, parameter :: jca34   = 95
  integer, parameter :: jca35   = 96
  integer, parameter :: jca36   = 97
  integer, parameter :: jca37   = 98
  integer, parameter :: jca38   = 99
  integer, parameter :: jca39   = 100
  integer, parameter :: jca40   = 101
  integer, parameter :: jca41   = 102
  integer, parameter :: jca42   = 103
  integer, parameter :: jca43   = 104
  integer, parameter :: jca44   = 105
  integer, parameter :: jsc36   = 106
  integer, parameter :: jsc37   = 107
  integer, parameter :: jsc38   = 108
  integer, parameter :: jsc39   = 109
  integer, parameter :: jsc40   = 110
  integer, parameter :: jsc41   = 111
  integer, parameter :: jsc42   = 112
  integer, parameter :: jsc43   = 113
  integer, parameter :: jsc44   = 114
  integer, parameter :: jsc45   = 115
  integer, parameter :: jti38   = 116
  integer, parameter :: jti39   = 117
  integer, parameter :: jti40   = 118
  integer, parameter :: jti41   = 119
  integer, parameter :: jti42   = 120
  integer, parameter :: jti43   = 121
  integer, parameter :: jti44   = 122
  integer, parameter :: jti45   = 123
  integer, parameter :: jti46   = 124
  integer, parameter :: jti47   = 125
  integer, parameter :: jti48   = 126
  integer, parameter :: jv40   = 127
  integer, parameter :: jv41   = 128
  integer, parameter :: jv42   = 129
  integer, parameter :: jv43   = 130
  integer, parameter :: jv44   = 131
  integer, parameter :: jv45   = 132
  integer, parameter :: jv46   = 133
  integer, parameter :: jv47   = 134
  integer, parameter :: jv48   = 135
  integer, parameter :: jv49   = 136
  integer, parameter :: jcr42   = 137
  integer, parameter :: jcr43   = 138
  integer, parameter :: jcr44   = 139
  integer, parameter :: jcr45   = 140
  integer, parameter :: jcr46   = 141
  integer, parameter :: jcr47   = 142
  integer, parameter :: jcr48   = 143
  integer, parameter :: jcr49   = 144
  integer, parameter :: jcr50   = 145
  integer, parameter :: jcr51   = 146
  integer, parameter :: jcr52   = 147
  integer, parameter :: jmn44   = 148
  integer, parameter :: jmn45   = 149
  integer, parameter :: jmn46   = 150
  integer, parameter :: jmn47   = 151
  integer, parameter :: jmn48   = 152
  integer, parameter :: jmn49   = 153
  integer, parameter :: jmn50   = 154
  integer, parameter :: jmn51   = 155
  integer, parameter :: jmn52   = 156
  integer, parameter :: jmn53   = 157
  integer, parameter :: jmn55   = 158
  integer, parameter :: jfe45   = 159
  integer, parameter :: jfe46   = 160
  integer, parameter :: jfe47   = 161
  integer, parameter :: jfe48   = 162
  integer, parameter :: jfe49   = 163
  integer, parameter :: jfe50   = 164
  integer, parameter :: jfe51   = 165
  integer, parameter :: jfe52   = 166
  integer, parameter :: jfe53   = 167
  integer, parameter :: jfe54   = 168
  integer, parameter :: jfe55   = 169
  integer, parameter :: jfe56   = 170
  integer, parameter :: jco47   = 171
  integer, parameter :: jco48   = 172
  integer, parameter :: jco49   = 173
  integer, parameter :: jco50   = 174
  integer, parameter :: jco51   = 175
  integer, parameter :: jco52   = 176
  integer, parameter :: jco53   = 177
  integer, parameter :: jco54   = 178
  integer, parameter :: jco55   = 179
  integer, parameter :: jco56   = 180
  integer, parameter :: jni48   = 181
  integer, parameter :: jni49   = 182
  integer, parameter :: jni50   = 183
  integer, parameter :: jni51   = 184
  integer, parameter :: jni52   = 185
  integer, parameter :: jni53   = 186
  integer, parameter :: jni54   = 187
  integer, parameter :: jni55   = 188
  integer, parameter :: jni56   = 189

  ! Reactions
  integer, parameter :: k_be7__li7__weak__electron_capture   = 1
  integer, parameter :: k_b8__be8__weak__wc17   = 2
  integer, parameter :: k_o14__n14__weak__wc12   = 3
  integer, parameter :: k_o15__n15__weak__wc12   = 4
  integer, parameter :: k_f17__o17__weak__wc12   = 5
  integer, parameter :: k_f18__o18__weak__wc12   = 6
  integer, parameter :: k_f20__ne20__weak__wc12   = 7
  integer, parameter :: k_ne18__f18__weak__wc12   = 8
  integer, parameter :: k_ne19__f19__weak__wc12   = 9
  integer, parameter :: k_he3__p_d   = 10
  integer, parameter :: k_he4__d_d   = 11
  integer, parameter :: k_be7__he4_he3   = 12
  integer, parameter :: k_b8__p_be7   = 13
  integer, parameter :: k_b8__he4_he4__weak__wc12   = 14
  integer, parameter :: k_n13__p_c12   = 15
  integer, parameter :: k_o14__p_n13   = 16
  integer, parameter :: k_o15__p_n14   = 17
  integer, parameter :: k_o16__p_n15   = 18
  integer, parameter :: k_o16__he4_c12   = 19
  integer, parameter :: k_f17__p_o16   = 20
  integer, parameter :: k_f18__p_o17   = 21
  integer, parameter :: k_f18__he4_n14   = 22
  integer, parameter :: k_f19__p_o18   = 23
  integer, parameter :: k_f19__he4_n15   = 24
  integer, parameter :: k_ne18__p_f17   = 25
  integer, parameter :: k_ne18__he4_o14   = 26
  integer, parameter :: k_ne19__p_f18   = 27
  integer, parameter :: k_ne19__he4_o15   = 28
  integer, parameter :: k_ne20__p_f19   = 29
  integer, parameter :: k_ne20__he4_o16   = 30
  integer, parameter :: k_ne21__p_f20   = 31
  integer, parameter :: k_ne21__he4_o17   = 32
  integer, parameter :: k_c12__he4_he4_he4   = 33
  integer, parameter :: k_p_p__d__weak__bet_pos_   = 34
  integer, parameter :: k_p_p__d__weak__electron_capture   = 35
  integer, parameter :: k_p_d__he3   = 36
  integer, parameter :: k_d_d__he4   = 37
  integer, parameter :: k_p_he3__he4__weak__bet_pos_   = 38
  integer, parameter :: k_he4_he3__be7   = 39
  integer, parameter :: k_p_be7__b8   = 40
  integer, parameter :: k_p_c12__n13   = 41
  integer, parameter :: k_he4_c12__o16   = 42
  integer, parameter :: k_p_n13__o14   = 43
  integer, parameter :: k_p_n14__o15   = 44
  integer, parameter :: k_he4_n14__f18   = 45
  integer, parameter :: k_p_n15__o16   = 46
  integer, parameter :: k_he4_n15__f19   = 47
  integer, parameter :: k_he4_o14__ne18   = 48
  integer, parameter :: k_he4_o15__ne19   = 49
  integer, parameter :: k_p_o16__f17   = 50
  integer, parameter :: k_he4_o16__ne20   = 51
  integer, parameter :: k_p_o17__f18   = 52
  integer, parameter :: k_he4_o17__ne21   = 53
  integer, parameter :: k_p_o18__f19   = 54
  integer, parameter :: k_p_f17__ne18   = 55
  integer, parameter :: k_p_f18__ne19   = 56
  integer, parameter :: k_p_f19__ne20   = 57
  integer, parameter :: k_p_f20__ne21   = 58
  integer, parameter :: k_d_he3__p_he4   = 59
  integer, parameter :: k_p_he4__d_he3   = 60
  integer, parameter :: k_he4_he4__p_li7   = 61
  integer, parameter :: k_p_li7__he4_he4   = 62
  integer, parameter :: k_he4_c12__p_n15   = 63
  integer, parameter :: k_c12_c12__he4_ne20   = 64
  integer, parameter :: k_he4_n13__p_o16   = 65
  integer, parameter :: k_he4_n14__p_o17   = 66
  integer, parameter :: k_p_n15__he4_c12   = 67
  integer, parameter :: k_he4_n15__p_o18   = 68
  integer, parameter :: k_he4_o14__p_f17   = 69
  integer, parameter :: k_he4_o15__p_f18   = 70
  integer, parameter :: k_p_o16__he4_n13   = 71
  integer, parameter :: k_he4_o16__p_f19   = 72
  integer, parameter :: k_p_o17__he4_n14   = 73
  integer, parameter :: k_he4_o17__p_f20   = 74
  integer, parameter :: k_p_o18__he4_n15   = 75
  integer, parameter :: k_p_f17__he4_o14   = 76
  integer, parameter :: k_he4_f17__p_ne20   = 77
  integer, parameter :: k_p_f18__he4_o15   = 78
  integer, parameter :: k_he4_f18__p_ne21   = 79
  integer, parameter :: k_p_f19__he4_o16   = 80
  integer, parameter :: k_p_f20__he4_o17   = 81
  integer, parameter :: k_p_ne20__he4_f17   = 82
  integer, parameter :: k_he4_ne20__c12_c12   = 83
  integer, parameter :: k_p_ne21__he4_f18   = 84
  integer, parameter :: k_he3_he3__p_p_he4   = 85
  integer, parameter :: k_d_be7__p_he4_he4   = 86
  integer, parameter :: k_he3_be7__p_p_he4_he4   = 87
  integer, parameter :: k_he4_he4_he4__c12   = 88
  integer, parameter :: k_p_p_he4__he3_he3   = 89
  integer, parameter :: k_p_he4_he4__d_be7   = 90
  integer, parameter :: k_p_p_he4_he4__he3_be7   = 91
  integer, parameter :: k_n13__c13__weak__wc12   = 92
  integer, parameter :: k_p_o15__f16   = 93
  integer, parameter :: k_he4_o18__ne22   = 94
  integer, parameter :: k_he4_f17__na21   = 95
  integer, parameter :: k_he4_f18__na22   = 96
  integer, parameter :: k_he4_f19__na23   = 97
  integer, parameter :: k_he4_f20__na24   = 98
  integer, parameter :: k_p_ne18__na19   = 99
  integer, parameter :: k_he4_ne18__mg22   = 100
  integer, parameter :: k_p_ne19__na20   = 101
  integer, parameter :: k_he4_ne19__mg23   = 102
  integer, parameter :: k_p_ne20__na21   = 103
  integer, parameter :: k_he4_ne20__mg24   = 104
  integer, parameter :: k_p_ne21__na22   = 105
  integer, parameter :: k_he4_ne21__mg25   = 106
  integer, parameter :: k_p_c13__n14   = 107
  integer, parameter :: k_he4_f16__na20   = 108
  integer, parameter :: k_p_ne22__na23   = 109
  integer, parameter :: k_he4_ne22__mg26   = 110
  integer, parameter :: k_na21__ne21__weak__wc12   = 111
  integer, parameter :: k_p_na21__mg22   = 112
  integer, parameter :: k_he4_na21__al25   = 113
  integer, parameter :: k_na22__ne22__weak__wc12   = 114
  integer, parameter :: k_p_na22__mg23   = 115
  integer, parameter :: k_he4_na22__al26   = 116
  integer, parameter :: k_p_na23__mg24   = 117
  integer, parameter :: k_he4_na23__al27   = 118
  integer, parameter :: k_p_na24__mg25   = 119
  integer, parameter :: k_he4_na24__al28   = 120
  integer, parameter :: k_mg22__na22__weak__wc12   = 121
  integer, parameter :: k_p_mg22__al23   = 122
  integer, parameter :: k_he4_mg22__si26   = 123
  integer, parameter :: k_na19__ne19__weak__bqa_pos_   = 124
  integer, parameter :: k_he4_na19__al23   = 125
  integer, parameter :: k_na20__ne20__weak__wc12   = 126
  integer, parameter :: k_p_na20__mg21   = 127
  integer, parameter :: k_he4_na20__al24   = 128
  integer, parameter :: k_mg23__na23__weak__wc12   = 129
  integer, parameter :: k_p_mg23__al24   = 130
  integer, parameter :: k_he4_mg23__si27   = 131
  integer, parameter :: k_p_mg24__al25   = 132
  integer, parameter :: k_he4_mg24__si28   = 133
  integer, parameter :: k_p_mg25__al26   = 134
  integer, parameter :: k_he4_mg25__si29   = 135
  integer, parameter :: k_p_mg26__al27   = 136
  integer, parameter :: k_he4_mg26__si30   = 137
  integer, parameter :: k_al25__mg25__weak__wc12   = 138
  integer, parameter :: k_p_al25__si26   = 139
  integer, parameter :: k_he4_al25__p29   = 140
  integer, parameter :: k_al26__mg26__weak__wc12   = 141
  integer, parameter :: k_p_al26__si27   = 142
  integer, parameter :: k_he4_al26__p30   = 143
  integer, parameter :: k_p_al27__si28   = 144
  integer, parameter :: k_he4_al27__p31   = 145
  integer, parameter :: k_p_al28__si29   = 146
  integer, parameter :: k_he4_al28__p32   = 147
  integer, parameter :: k_si26__al26__weak__wc12   = 148
  integer, parameter :: k_p_si26__p27   = 149
  integer, parameter :: k_he4_si26__s30   = 150
  integer, parameter :: k_al23__mg23__weak__wc12   = 151
  integer, parameter :: k_p_al23__si24   = 152
  integer, parameter :: k_he4_al23__p27   = 153
  integer, parameter :: k_al24__mg24__weak__wc12   = 154
  integer, parameter :: k_p_al24__si25   = 155
  integer, parameter :: k_he4_al24__p28   = 156
  integer, parameter :: k_mg21__na21__weak__wc12   = 157
  integer, parameter :: k_p_mg21__al22   = 158
  integer, parameter :: k_he4_mg21__si25   = 159
  integer, parameter :: k_si27__al27__weak__wc12   = 160
  integer, parameter :: k_p_si27__p28   = 161
  integer, parameter :: k_he4_si27__s31   = 162
  integer, parameter :: k_p_si28__p29   = 163
  integer, parameter :: k_he4_si28__s32   = 164
  integer, parameter :: k_p_si29__p30   = 165
  integer, parameter :: k_he4_si29__s33   = 166
  integer, parameter :: k_p_si30__p31   = 167
  integer, parameter :: k_he4_si30__s34   = 168
  integer, parameter :: k_p29__si29__weak__wc12   = 169
  integer, parameter :: k_p_p29__s30   = 170
  integer, parameter :: k_he4_p29__cl33   = 171
  integer, parameter :: k_p30__si30__weak__wc12   = 172
  integer, parameter :: k_p_p30__s31   = 173
  integer, parameter :: k_he4_p30__cl34   = 174
  integer, parameter :: k_p_p31__s32   = 175
  integer, parameter :: k_he4_p31__cl35   = 176
  integer, parameter :: k_p_p32__s33   = 177
  integer, parameter :: k_he4_p32__cl36   = 178
  integer, parameter :: k_s30__p30__weak__wc12   = 179
  integer, parameter :: k_p_s30__cl31   = 180
  integer, parameter :: k_he4_s30__ar34   = 181
  integer, parameter :: k_p27__si27__weak__wc12   = 182
  integer, parameter :: k_p_p27__s28   = 183
  integer, parameter :: k_he4_p27__cl31   = 184
  integer, parameter :: k_si24__al24__weak__wc12   = 185
  integer, parameter :: k_he4_si24__s28   = 186
  integer, parameter :: k_p28__si28__weak__wc12   = 187
  integer, parameter :: k_p_p28__s29   = 188
  integer, parameter :: k_he4_p28__cl32   = 189
  integer, parameter :: k_si25__al25__weak__wc12   = 190
  integer, parameter :: k_p_si25__p26   = 191
  integer, parameter :: k_he4_si25__s29   = 192
  integer, parameter :: k_al22__mg22__weak__wc12   = 193
  integer, parameter :: k_he4_al22__p26   = 194
  integer, parameter :: k_s31__p31__weak__wc12   = 195
  integer, parameter :: k_p_s31__cl32   = 196
  integer, parameter :: k_he4_s31__ar35   = 197
  integer, parameter :: k_p_s32__cl33   = 198
  integer, parameter :: k_he4_s32__ar36   = 199
  integer, parameter :: k_p_s33__cl34   = 200
  integer, parameter :: k_he4_s33__ar37   = 201
  integer, parameter :: k_p_s34__cl35   = 202
  integer, parameter :: k_he4_s34__ar38   = 203
  integer, parameter :: k_cl33__s33__weak__wc12   = 204
  integer, parameter :: k_p_cl33__ar34   = 205
  integer, parameter :: k_he4_cl33__k37   = 206
  integer, parameter :: k_cl34__s34__weak__wc12   = 207
  integer, parameter :: k_p_cl34__ar35   = 208
  integer, parameter :: k_he4_cl34__k38   = 209
  integer, parameter :: k_p_cl35__ar36   = 210
  integer, parameter :: k_he4_cl35__k39   = 211
  integer, parameter :: k_p_cl36__ar37   = 212
  integer, parameter :: k_he4_cl36__k40   = 213
  integer, parameter :: k_cl31__s31__weak__wc12   = 214
  integer, parameter :: k_p_cl31__ar32   = 215
  integer, parameter :: k_he4_cl31__k35   = 216
  integer, parameter :: k_ar34__cl34__weak__wc12   = 217
  integer, parameter :: k_p_ar34__k35   = 218
  integer, parameter :: k_he4_ar34__ca38   = 219
  integer, parameter :: k_s28__p28__weak__wc12   = 220
  integer, parameter :: k_p_s28__cl29   = 221
  integer, parameter :: k_he4_s28__ar32   = 222
  integer, parameter :: k_cl32__s32__weak__wc12   = 223
  integer, parameter :: k_p_cl32__ar33   = 224
  integer, parameter :: k_he4_cl32__k36   = 225
  integer, parameter :: k_s29__p29__weak__wc12   = 226
  integer, parameter :: k_p_s29__cl30   = 227
  integer, parameter :: k_he4_s29__ar33   = 228
  integer, parameter :: k_p26__si26__weak__wc12   = 229
  integer, parameter :: k_he4_p26__cl30   = 230
  integer, parameter :: k_ar35__cl35__weak__wc12   = 231
  integer, parameter :: k_p_ar35__k36   = 232
  integer, parameter :: k_he4_ar35__ca39   = 233
  integer, parameter :: k_p_ar36__k37   = 234
  integer, parameter :: k_he4_ar36__ca40   = 235
  integer, parameter :: k_ar37__cl37__weak__wc12   = 236
  integer, parameter :: k_p_ar37__k38   = 237
  integer, parameter :: k_he4_ar37__ca41   = 238
  integer, parameter :: k_p_ar38__k39   = 239
  integer, parameter :: k_he4_ar38__ca42   = 240
  integer, parameter :: k_k37__ar37__weak__wc12   = 241
  integer, parameter :: k_p_k37__ca38   = 242
  integer, parameter :: k_he4_k37__sc41   = 243
  integer, parameter :: k_k38__ar38__weak__wc12   = 244
  integer, parameter :: k_p_k38__ca39   = 245
  integer, parameter :: k_he4_k38__sc42   = 246
  integer, parameter :: k_p_k39__ca40   = 247
  integer, parameter :: k_he4_k39__sc43   = 248
  integer, parameter :: k_p_k40__ca41   = 249
  integer, parameter :: k_he4_k40__sc44   = 250
  integer, parameter :: k_ar32__cl32__weak__wc12   = 251
  integer, parameter :: k_p_ar32__k33   = 252
  integer, parameter :: k_he4_ar32__ca36   = 253
  integer, parameter :: k_k35__ar35__weak__wc12   = 254
  integer, parameter :: k_p_k35__ca36   = 255
  integer, parameter :: k_he4_k35__sc39   = 256
  integer, parameter :: k_ca38__k38__weak__wc12   = 257
  integer, parameter :: k_p_ca38__sc39   = 258
  integer, parameter :: k_he4_ca38__ti42   = 259
  integer, parameter :: k_cl29__s29__weak__bqa_pos_   = 260
  integer, parameter :: k_he4_cl29__k33   = 261
  integer, parameter :: k_ar33__cl33__weak__wc12   = 262
  integer, parameter :: k_p_ar33__k34   = 263
  integer, parameter :: k_he4_ar33__ca37   = 264
  integer, parameter :: k_k36__ar36__weak__wc12   = 265
  integer, parameter :: k_p_k36__ca37   = 266
  integer, parameter :: k_he4_k36__sc40   = 267
  integer, parameter :: k_cl30__s30__weak__bqa_pos_   = 268
  integer, parameter :: k_p_cl30__ar31   = 269
  integer, parameter :: k_he4_cl30__k34   = 270
  integer, parameter :: k_ca39__k39__weak__wc12   = 271
  integer, parameter :: k_p_ca39__sc40   = 272
  integer, parameter :: k_he4_ca39__ti43   = 273
  integer, parameter :: k_p_ca40__sc41   = 274
  integer, parameter :: k_he4_ca40__ti44   = 275
  integer, parameter :: k_p_cl37__ar38   = 276
  integer, parameter :: k_he4_cl37__k41   = 277
  integer, parameter :: k_ca41__k41__weak__wc12   = 278
  integer, parameter :: k_p_ca41__sc42   = 279
  integer, parameter :: k_he4_ca41__ti45   = 280
  integer, parameter :: k_p_ca42__sc43   = 281
  integer, parameter :: k_he4_ca42__ti46   = 282
  integer, parameter :: k_sc41__ca41__weak__wc12   = 283
  integer, parameter :: k_p_sc41__ti42   = 284
  integer, parameter :: k_he4_sc41__v45   = 285
  integer, parameter :: k_sc42__ca42__weak__wc12   = 286
  integer, parameter :: k_p_sc42__ti43   = 287
  integer, parameter :: k_he4_sc42__v46   = 288
  integer, parameter :: k_sc43__ca43__weak__wc12   = 289
  integer, parameter :: k_p_sc43__ti44   = 290
  integer, parameter :: k_he4_sc43__v47   = 291
  integer, parameter :: k_sc44__ca44__weak__wc12   = 292
  integer, parameter :: k_p_sc44__ti45   = 293
  integer, parameter :: k_he4_sc44__v48   = 294
  integer, parameter :: k_ca36__k36__weak__wc12   = 295
  integer, parameter :: k_p_ca36__sc37   = 296
  integer, parameter :: k_he4_ca36__ti40   = 297
  integer, parameter :: k_k33__ar33__weak__bqa_pos_   = 298
  integer, parameter :: k_p_k33__ca34   = 299
  integer, parameter :: k_he4_k33__sc37   = 300
  integer, parameter :: k_sc39__ca39__weak__mo97   = 301
  integer, parameter :: k_p_sc39__ti40   = 302
  integer, parameter :: k_he4_sc39__v43   = 303
  integer, parameter :: k_ti42__sc42__weak__wc12   = 304
  integer, parameter :: k_p_ti42__v43   = 305
  integer, parameter :: k_he4_ti42__cr46   = 306
  integer, parameter :: k_k34__ar34__weak__bqa_pos_   = 307
  integer, parameter :: k_p_k34__ca35   = 308
  integer, parameter :: k_he4_k34__sc38   = 309
  integer, parameter :: k_ca37__k37__weak__wc12   = 310
  integer, parameter :: k_p_ca37__sc38   = 311
  integer, parameter :: k_he4_ca37__ti41   = 312
  integer, parameter :: k_sc40__ca40__weak__wc12   = 313
  integer, parameter :: k_p_sc40__ti41   = 314
  integer, parameter :: k_he4_sc40__v44   = 315
  integer, parameter :: k_ar31__cl31__weak__wc12   = 316
  integer, parameter :: k_he4_ar31__ca35   = 317
  integer, parameter :: k_ti43__sc43__weak__wc12   = 318
  integer, parameter :: k_p_ti43__v44   = 319
  integer, parameter :: k_he4_ti43__cr47   = 320
  integer, parameter :: k_ti44__sc44__weak__wc12   = 321
  integer, parameter :: k_p_ti44__v45   = 322
  integer, parameter :: k_he4_ti44__cr48   = 323
  integer, parameter :: k_p_k41__ca42   = 324
  integer, parameter :: k_he4_k41__sc45   = 325
  integer, parameter :: k_ti45__sc45__weak__wc12   = 326
  integer, parameter :: k_p_ti45__v46   = 327
  integer, parameter :: k_he4_ti45__cr49   = 328
  integer, parameter :: k_p_ti46__v47   = 329
  integer, parameter :: k_he4_ti46__cr50   = 330
  integer, parameter :: k_v45__ti45__weak__wc12   = 331
  integer, parameter :: k_p_v45__cr46   = 332
  integer, parameter :: k_he4_v45__mn49   = 333
  integer, parameter :: k_v46__ti46__weak__wc12   = 334
  integer, parameter :: k_p_v46__cr47   = 335
  integer, parameter :: k_he4_v46__mn50   = 336
  integer, parameter :: k_p_ca43__sc44   = 337
  integer, parameter :: k_he4_ca43__ti47   = 338
  integer, parameter :: k_v47__ti47__weak__wc12   = 339
  integer, parameter :: k_p_v47__cr48   = 340
  integer, parameter :: k_he4_v47__mn51   = 341
  integer, parameter :: k_v48__ti48__weak__wc12   = 342
  integer, parameter :: k_p_v48__cr49   = 343
  integer, parameter :: k_he4_v48__mn52   = 344
  integer, parameter :: k_p_ca44__sc45   = 345
  integer, parameter :: k_he4_ca44__ti48   = 346
  integer, parameter :: k_ti40__sc40__weak__wc17   = 347
  integer, parameter :: k_p_ti40__v41   = 348
  integer, parameter :: k_he4_ti40__cr44   = 349
  integer, parameter :: k_sc37__ca37__weak__bqa_pos_   = 350
  integer, parameter :: k_p_sc37__ti38   = 351
  integer, parameter :: k_he4_sc37__v41   = 352
  integer, parameter :: k_ca34__k34__weak__bqa_pos_   = 353
  integer, parameter :: k_he4_ca34__ti38   = 354
  integer, parameter :: k_v43__ti43__weak__wc12   = 355
  integer, parameter :: k_p_v43__cr44   = 356
  integer, parameter :: k_he4_v43__mn47   = 357
  integer, parameter :: k_cr46__v46__weak__wc12   = 358
  integer, parameter :: k_p_cr46__mn47   = 359
  integer, parameter :: k_he4_cr46__fe50   = 360
  integer, parameter :: k_ca35__k35__weak__wc17   = 361
  integer, parameter :: k_p_ca35__sc36   = 362
  integer, parameter :: k_he4_ca35__ti39   = 363
  integer, parameter :: k_sc38__ca38__weak__mo97   = 364
  integer, parameter :: k_p_sc38__ti39   = 365
  integer, parameter :: k_he4_sc38__v42   = 366
  integer, parameter :: k_ti41__sc41__weak__wc17   = 367
  integer, parameter :: k_p_ti41__v42   = 368
  integer, parameter :: k_he4_ti41__cr45   = 369
  integer, parameter :: k_v44__ti44__weak__wc12   = 370
  integer, parameter :: k_p_v44__cr45   = 371
  integer, parameter :: k_he4_v44__mn48   = 372
  integer, parameter :: k_cr47__v47__weak__wc12   = 373
  integer, parameter :: k_p_cr47__mn48   = 374
  integer, parameter :: k_he4_cr47__fe51   = 375
  integer, parameter :: k_cr48__v48__weak__wc12   = 376
  integer, parameter :: k_p_cr48__mn49   = 377
  integer, parameter :: k_he4_cr48__fe52   = 378
  integer, parameter :: k_p_sc45__ti46   = 379
  integer, parameter :: k_he4_sc45__v49   = 380
  integer, parameter :: k_cr49__v49__weak__wc12   = 381
  integer, parameter :: k_p_cr49__mn50   = 382
  integer, parameter :: k_he4_cr49__fe53   = 383
  integer, parameter :: k_p_cr50__mn51   = 384
  integer, parameter :: k_he4_cr50__fe54   = 385
  integer, parameter :: k_mn49__cr49__weak__wc12   = 386
  integer, parameter :: k_p_mn49__fe50   = 387
  integer, parameter :: k_he4_mn49__co53   = 388
  integer, parameter :: k_mn50__cr50__weak__wc12   = 389
  integer, parameter :: k_p_mn50__fe51   = 390
  integer, parameter :: k_he4_mn50__co54   = 391
  integer, parameter :: k_p_ti47__v48   = 392
  integer, parameter :: k_he4_ti47__cr51   = 393
  integer, parameter :: k_mn51__cr51__weak__wc12   = 394
  integer, parameter :: k_p_mn51__fe52   = 395
  integer, parameter :: k_he4_mn51__co55   = 396
  integer, parameter :: k_p_ti48__v49   = 397
  integer, parameter :: k_he4_ti48__cr52   = 398
  integer, parameter :: k_mn52__cr52__weak__wc12   = 399
  integer, parameter :: k_p_mn52__fe53   = 400
  integer, parameter :: k_he4_mn52__co56   = 401
  integer, parameter :: k_cr44__v44__weak__wc12   = 402
  integer, parameter :: k_p_cr44__mn45   = 403
  integer, parameter :: k_he4_cr44__fe48   = 404
  integer, parameter :: k_v41__ti41__weak__bqa_pos_   = 405
  integer, parameter :: k_p_v41__cr42   = 406
  integer, parameter :: k_he4_v41__mn45   = 407
  integer, parameter :: k_ti38__sc38__weak__mo97   = 408
  integer, parameter :: k_he4_ti38__cr42   = 409
  integer, parameter :: k_mn47__cr47__weak__wc17   = 410
  integer, parameter :: k_p_mn47__fe48   = 411
  integer, parameter :: k_he4_mn47__co51   = 412
  integer, parameter :: k_fe50__mn50__weak__wc12   = 413
  integer, parameter :: k_p_fe50__co51   = 414
  integer, parameter :: k_he4_fe50__ni54   = 415
  integer, parameter :: k_sc36__ca36__weak__bqa_pos_   = 416
  integer, parameter :: k_he4_sc36__v40   = 417
  integer, parameter :: k_ti39__sc39__weak__wc17   = 418
  integer, parameter :: k_p_ti39__v40   = 419
  integer, parameter :: k_he4_ti39__cr43   = 420
  integer, parameter :: k_v42__ti42__weak__mo97   = 421
  integer, parameter :: k_p_v42__cr43   = 422
  integer, parameter :: k_he4_v42__mn46   = 423
  integer, parameter :: k_cr45__v45__weak__wc12   = 424
  integer, parameter :: k_p_cr45__mn46   = 425
  integer, parameter :: k_he4_cr45__fe49   = 426
  integer, parameter :: k_mn48__cr48__weak__wc12   = 427
  integer, parameter :: k_p_mn48__fe49   = 428
  integer, parameter :: k_he4_mn48__co52   = 429
  integer, parameter :: k_fe51__mn51__weak__wc12   = 430
  integer, parameter :: k_p_fe51__co52   = 431
  integer, parameter :: k_he4_fe51__ni55   = 432
  integer, parameter :: k_fe52__mn52__weak__wc12   = 433
  integer, parameter :: k_p_fe52__co53   = 434
  integer, parameter :: k_he4_fe52__ni56   = 435
  integer, parameter :: k_p_v49__cr50   = 436
  integer, parameter :: k_he4_v49__mn53   = 437
  integer, parameter :: k_fe53__mn53__weak__wc12   = 438
  integer, parameter :: k_p_fe53__co54   = 439
  integer, parameter :: k_p_fe54__co55   = 440
  integer, parameter :: k_co53__fe53__weak__wc12   = 441
  integer, parameter :: k_p_co53__ni54   = 442
  integer, parameter :: k_co54__fe54__weak__wc12   = 443
  integer, parameter :: k_p_co54__ni55   = 444
  integer, parameter :: k_p_cr51__mn52   = 445
  integer, parameter :: k_he4_cr51__fe55   = 446
  integer, parameter :: k_co55__fe55__weak__wc12   = 447
  integer, parameter :: k_p_co55__ni56   = 448
  integer, parameter :: k_p_cr52__mn53   = 449
  integer, parameter :: k_he4_cr52__fe56   = 450
  integer, parameter :: k_co56__fe56__weak__wc12   = 451
  integer, parameter :: k_mn45__cr45__weak__bqa_pos_   = 452
  integer, parameter :: k_p_mn45__fe46   = 453
  integer, parameter :: k_he4_mn45__co49   = 454
  integer, parameter :: k_fe48__mn48__weak__wc12   = 455
  integer, parameter :: k_p_fe48__co49   = 456
  integer, parameter :: k_he4_fe48__ni52   = 457
  integer, parameter :: k_cr42__v42__weak__wc12   = 458
  integer, parameter :: k_he4_cr42__fe46   = 459
  integer, parameter :: k_co51__fe51__weak__mo97   = 460
  integer, parameter :: k_p_co51__ni52   = 461
  integer, parameter :: k_ni54__co54__weak__wc12   = 462
  integer, parameter :: k_v40__ti40__weak__bqa_pos_   = 463
  integer, parameter :: k_he4_v40__mn44   = 464
  integer, parameter :: k_cr43__v43__weak__wc12   = 465
  integer, parameter :: k_p_cr43__mn44   = 466
  integer, parameter :: k_he4_cr43__fe47   = 467
  integer, parameter :: k_mn46__cr46__weak__wc12   = 468
  integer, parameter :: k_p_mn46__fe47   = 469
  integer, parameter :: k_he4_mn46__co50   = 470
  integer, parameter :: k_fe49__mn49__weak__wc12   = 471
  integer, parameter :: k_p_fe49__co50   = 472
  integer, parameter :: k_he4_fe49__ni53   = 473
  integer, parameter :: k_co52__fe52__weak__wc12   = 474
  integer, parameter :: k_p_co52__ni53   = 475
  integer, parameter :: k_ni55__co55__weak__wc12   = 476
  integer, parameter :: k_ni56__co56__weak__wc12   = 477
  integer, parameter :: k_p_mn53__fe54   = 478
  integer, parameter :: k_fe55__mn55__weak__wc12   = 479
  integer, parameter :: k_p_fe55__co56   = 480
  integer, parameter :: k_fe46__mn46__weak__wc12   = 481
  integer, parameter :: k_p_fe46__co47   = 482
  integer, parameter :: k_he4_fe46__ni50   = 483
  integer, parameter :: k_co49__fe49__weak__bqa_pos_   = 484
  integer, parameter :: k_p_co49__ni50   = 485
  integer, parameter :: k_ni52__co52__weak__wc12   = 486
  integer, parameter :: k_mn44__cr44__weak__bqa_pos_   = 487
  integer, parameter :: k_p_mn44__fe45   = 488
  integer, parameter :: k_he4_mn44__co48   = 489
  integer, parameter :: k_fe47__mn47__weak__wc12   = 490
  integer, parameter :: k_p_fe47__co48   = 491
  integer, parameter :: k_he4_fe47__ni51   = 492
  integer, parameter :: k_co50__fe50__weak__wc12   = 493
  integer, parameter :: k_p_co50__ni51   = 494
  integer, parameter :: k_ni53__co53__weak__wc12   = 495
  integer, parameter :: k_p_mn55__fe56   = 496
  integer, parameter :: k_ni50__co50__weak__wc12   = 497
  integer, parameter :: k_co47__fe47__weak__bqa_pos_   = 498
  integer, parameter :: k_p_co47__ni48   = 499
  integer, parameter :: k_fe45__mn45__weak__wc17   = 500
  integer, parameter :: k_he4_fe45__ni49   = 501
  integer, parameter :: k_co48__fe48__weak__bqa_pos_   = 502
  integer, parameter :: k_p_co48__ni49   = 503
  integer, parameter :: k_ni51__co51__weak__wc17   = 504
  integer, parameter :: k_ni48__co48__weak__wc17   = 505
  integer, parameter :: k_ni49__co49__weak__wc12   = 506

  ! reactvec indices
  integer, parameter :: i_rate        = 1
  integer, parameter :: i_drate_dt    = 2
  integer, parameter :: i_scor        = 3
  integer, parameter :: i_dscor_dt    = 4
  integer, parameter :: i_dqweak      = 5
  integer, parameter :: i_epart       = 6

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(rt), allocatable, save :: aion(:), zion(:), bion(:)
  real(rt), allocatable, save :: nion(:), mion(:), wion(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: aion, zion, bion, nion, mion, wion
#endif

  !$acc declare create(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
  ! Shape of Jacobian in Compressed Sparse Row format
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 1940
  integer, allocatable :: csr_jac_col_index(:), csr_jac_row_count(:)

#ifdef AMREX_USE_CUDA
  attributes(managed) :: csr_jac_col_index, csr_jac_row_count
#endif
#endif

contains

  subroutine actual_network_init()

    implicit none

    integer :: i

    ! Allocate ion info arrays
    allocate(aion(nspec))
    allocate(zion(nspec))
    allocate(bion(nspec))
    allocate(nion(nspec))
    allocate(mion(nspec))
    allocate(wion(nspec))

    spec_names(jp)   = "hydrogen-1"
    spec_names(jd)   = "hydrogen-2"
    spec_names(jhe3)   = "helium-3"
    spec_names(jhe4)   = "helium-4"
    spec_names(jli7)   = "lithium-7"
    spec_names(jbe7)   = "beryllium-7"
    spec_names(jbe8)   = "beryllium-8"
    spec_names(jb8)   = "boron-8"
    spec_names(jc12)   = "carbon-12"
    spec_names(jc13)   = "carbon-13"
    spec_names(jn13)   = "nitrogen-13"
    spec_names(jn14)   = "nitrogen-14"
    spec_names(jn15)   = "nitrogen-15"
    spec_names(jo14)   = "oxygen-14"
    spec_names(jo15)   = "oxygen-15"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jo17)   = "oxygen-17"
    spec_names(jo18)   = "oxygen-18"
    spec_names(jf16)   = "fluorine-16"
    spec_names(jf17)   = "fluorine-17"
    spec_names(jf18)   = "fluorine-18"
    spec_names(jf19)   = "fluorine-19"
    spec_names(jf20)   = "fluorine-20"
    spec_names(jne18)   = "neon-18"
    spec_names(jne19)   = "neon-19"
    spec_names(jne20)   = "neon-20"
    spec_names(jne21)   = "neon-21"
    spec_names(jne22)   = "neon-22"
    spec_names(jna19)   = "sodium-19"
    spec_names(jna20)   = "sodium-20"
    spec_names(jna21)   = "sodium-21"
    spec_names(jna22)   = "sodium-22"
    spec_names(jna23)   = "sodium-23"
    spec_names(jna24)   = "sodium-24"
    spec_names(jmg21)   = "magnesium-21"
    spec_names(jmg22)   = "magnesium-22"
    spec_names(jmg23)   = "magnesium-23"
    spec_names(jmg24)   = "magnesium-24"
    spec_names(jmg25)   = "magnesium-25"
    spec_names(jmg26)   = "magnesium-26"
    spec_names(jal22)   = "aluminum-22"
    spec_names(jal23)   = "aluminum-23"
    spec_names(jal24)   = "aluminum-24"
    spec_names(jal25)   = "aluminum-25"
    spec_names(jal26)   = "aluminum-26"
    spec_names(jal27)   = "aluminum-27"
    spec_names(jal28)   = "aluminum-28"
    spec_names(jsi24)   = "silicon-24"
    spec_names(jsi25)   = "silicon-25"
    spec_names(jsi26)   = "silicon-26"
    spec_names(jsi27)   = "silicon-27"
    spec_names(jsi28)   = "silicon-28"
    spec_names(jsi29)   = "silicon-29"
    spec_names(jsi30)   = "silicon-30"
    spec_names(jp26)   = "phosphorus-26"
    spec_names(jp27)   = "phosphorus-27"
    spec_names(jp28)   = "phosphorus-28"
    spec_names(jp29)   = "phosphorus-29"
    spec_names(jp30)   = "phosphorus-30"
    spec_names(jp31)   = "phosphorus-31"
    spec_names(jp32)   = "phosphorus-32"
    spec_names(js28)   = "sulfur-28"
    spec_names(js29)   = "sulfur-29"
    spec_names(js30)   = "sulfur-30"
    spec_names(js31)   = "sulfur-31"
    spec_names(js32)   = "sulfur-32"
    spec_names(js33)   = "sulfur-33"
    spec_names(js34)   = "sulfur-34"
    spec_names(jcl29)   = "chlorine-29"
    spec_names(jcl30)   = "chlorine-30"
    spec_names(jcl31)   = "chlorine-31"
    spec_names(jcl32)   = "chlorine-32"
    spec_names(jcl33)   = "chlorine-33"
    spec_names(jcl34)   = "chlorine-34"
    spec_names(jcl35)   = "chlorine-35"
    spec_names(jcl36)   = "chlorine-36"
    spec_names(jcl37)   = "chlorine-37"
    spec_names(jar31)   = "argon-31"
    spec_names(jar32)   = "argon-32"
    spec_names(jar33)   = "argon-33"
    spec_names(jar34)   = "argon-34"
    spec_names(jar35)   = "argon-35"
    spec_names(jar36)   = "argon-36"
    spec_names(jar37)   = "argon-37"
    spec_names(jar38)   = "argon-38"
    spec_names(jk33)   = "potassium-33"
    spec_names(jk34)   = "potassium-34"
    spec_names(jk35)   = "potassium-35"
    spec_names(jk36)   = "potassium-36"
    spec_names(jk37)   = "potassium-37"
    spec_names(jk38)   = "potassium-38"
    spec_names(jk39)   = "potassium-39"
    spec_names(jk40)   = "potassium-40"
    spec_names(jk41)   = "potassium-41"
    spec_names(jca34)   = "calcium-34"
    spec_names(jca35)   = "calcium-35"
    spec_names(jca36)   = "calcium-36"
    spec_names(jca37)   = "calcium-37"
    spec_names(jca38)   = "calcium-38"
    spec_names(jca39)   = "calcium-39"
    spec_names(jca40)   = "calcium-40"
    spec_names(jca41)   = "calcium-41"
    spec_names(jca42)   = "calcium-42"
    spec_names(jca43)   = "calcium-43"
    spec_names(jca44)   = "calcium-44"
    spec_names(jsc36)   = "scandium-36"
    spec_names(jsc37)   = "scandium-37"
    spec_names(jsc38)   = "scandium-38"
    spec_names(jsc39)   = "scandium-39"
    spec_names(jsc40)   = "scandium-40"
    spec_names(jsc41)   = "scandium-41"
    spec_names(jsc42)   = "scandium-42"
    spec_names(jsc43)   = "scandium-43"
    spec_names(jsc44)   = "scandium-44"
    spec_names(jsc45)   = "scandium-45"
    spec_names(jti38)   = "titanium-38"
    spec_names(jti39)   = "titanium-39"
    spec_names(jti40)   = "titanium-40"
    spec_names(jti41)   = "titanium-41"
    spec_names(jti42)   = "titanium-42"
    spec_names(jti43)   = "titanium-43"
    spec_names(jti44)   = "titanium-44"
    spec_names(jti45)   = "titanium-45"
    spec_names(jti46)   = "titanium-46"
    spec_names(jti47)   = "titanium-47"
    spec_names(jti48)   = "titanium-48"
    spec_names(jv40)   = "vanadium-40"
    spec_names(jv41)   = "vanadium-41"
    spec_names(jv42)   = "vanadium-42"
    spec_names(jv43)   = "vanadium-43"
    spec_names(jv44)   = "vanadium-44"
    spec_names(jv45)   = "vanadium-45"
    spec_names(jv46)   = "vanadium-46"
    spec_names(jv47)   = "vanadium-47"
    spec_names(jv48)   = "vanadium-48"
    spec_names(jv49)   = "vanadium-49"
    spec_names(jcr42)   = "chromium-42"
    spec_names(jcr43)   = "chromium-43"
    spec_names(jcr44)   = "chromium-44"
    spec_names(jcr45)   = "chromium-45"
    spec_names(jcr46)   = "chromium-46"
    spec_names(jcr47)   = "chromium-47"
    spec_names(jcr48)   = "chromium-48"
    spec_names(jcr49)   = "chromium-49"
    spec_names(jcr50)   = "chromium-50"
    spec_names(jcr51)   = "chromium-51"
    spec_names(jcr52)   = "chromium-52"
    spec_names(jmn44)   = "manganese-44"
    spec_names(jmn45)   = "manganese-45"
    spec_names(jmn46)   = "manganese-46"
    spec_names(jmn47)   = "manganese-47"
    spec_names(jmn48)   = "manganese-48"
    spec_names(jmn49)   = "manganese-49"
    spec_names(jmn50)   = "manganese-50"
    spec_names(jmn51)   = "manganese-51"
    spec_names(jmn52)   = "manganese-52"
    spec_names(jmn53)   = "manganese-53"
    spec_names(jmn55)   = "manganese-55"
    spec_names(jfe45)   = "iron-45"
    spec_names(jfe46)   = "iron-46"
    spec_names(jfe47)   = "iron-47"
    spec_names(jfe48)   = "iron-48"
    spec_names(jfe49)   = "iron-49"
    spec_names(jfe50)   = "iron-50"
    spec_names(jfe51)   = "iron-51"
    spec_names(jfe52)   = "iron-52"
    spec_names(jfe53)   = "iron-53"
    spec_names(jfe54)   = "iron-54"
    spec_names(jfe55)   = "iron-55"
    spec_names(jfe56)   = "iron-56"
    spec_names(jco47)   = "cobalt-47"
    spec_names(jco48)   = "cobalt-48"
    spec_names(jco49)   = "cobalt-49"
    spec_names(jco50)   = "cobalt-50"
    spec_names(jco51)   = "cobalt-51"
    spec_names(jco52)   = "cobalt-52"
    spec_names(jco53)   = "cobalt-53"
    spec_names(jco54)   = "cobalt-54"
    spec_names(jco55)   = "cobalt-55"
    spec_names(jco56)   = "cobalt-56"
    spec_names(jni48)   = "nickel-48"
    spec_names(jni49)   = "nickel-49"
    spec_names(jni50)   = "nickel-50"
    spec_names(jni51)   = "nickel-51"
    spec_names(jni52)   = "nickel-52"
    spec_names(jni53)   = "nickel-53"
    spec_names(jni54)   = "nickel-54"
    spec_names(jni55)   = "nickel-55"
    spec_names(jni56)   = "nickel-56"

    short_spec_names(jp)   = "h1"
    short_spec_names(jd)   = "h2"
    short_spec_names(jhe3)   = "he3"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jli7)   = "li7"
    short_spec_names(jbe7)   = "be7"
    short_spec_names(jbe8)   = "be8"
    short_spec_names(jb8)   = "b8"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jc13)   = "c13"
    short_spec_names(jn13)   = "n13"
    short_spec_names(jn14)   = "n14"
    short_spec_names(jn15)   = "n15"
    short_spec_names(jo14)   = "o14"
    short_spec_names(jo15)   = "o15"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jo17)   = "o17"
    short_spec_names(jo18)   = "o18"
    short_spec_names(jf16)   = "f16"
    short_spec_names(jf17)   = "f17"
    short_spec_names(jf18)   = "f18"
    short_spec_names(jf19)   = "f19"
    short_spec_names(jf20)   = "f20"
    short_spec_names(jne18)   = "ne18"
    short_spec_names(jne19)   = "ne19"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jne21)   = "ne21"
    short_spec_names(jne22)   = "ne22"
    short_spec_names(jna19)   = "na19"
    short_spec_names(jna20)   = "na20"
    short_spec_names(jna21)   = "na21"
    short_spec_names(jna22)   = "na22"
    short_spec_names(jna23)   = "na23"
    short_spec_names(jna24)   = "na24"
    short_spec_names(jmg21)   = "mg21"
    short_spec_names(jmg22)   = "mg22"
    short_spec_names(jmg23)   = "mg23"
    short_spec_names(jmg24)   = "mg24"
    short_spec_names(jmg25)   = "mg25"
    short_spec_names(jmg26)   = "mg26"
    short_spec_names(jal22)   = "al22"
    short_spec_names(jal23)   = "al23"
    short_spec_names(jal24)   = "al24"
    short_spec_names(jal25)   = "al25"
    short_spec_names(jal26)   = "al26"
    short_spec_names(jal27)   = "al27"
    short_spec_names(jal28)   = "al28"
    short_spec_names(jsi24)   = "si24"
    short_spec_names(jsi25)   = "si25"
    short_spec_names(jsi26)   = "si26"
    short_spec_names(jsi27)   = "si27"
    short_spec_names(jsi28)   = "si28"
    short_spec_names(jsi29)   = "si29"
    short_spec_names(jsi30)   = "si30"
    short_spec_names(jp26)   = "p26"
    short_spec_names(jp27)   = "p27"
    short_spec_names(jp28)   = "p28"
    short_spec_names(jp29)   = "p29"
    short_spec_names(jp30)   = "p30"
    short_spec_names(jp31)   = "p31"
    short_spec_names(jp32)   = "p32"
    short_spec_names(js28)   = "s28"
    short_spec_names(js29)   = "s29"
    short_spec_names(js30)   = "s30"
    short_spec_names(js31)   = "s31"
    short_spec_names(js32)   = "s32"
    short_spec_names(js33)   = "s33"
    short_spec_names(js34)   = "s34"
    short_spec_names(jcl29)   = "cl29"
    short_spec_names(jcl30)   = "cl30"
    short_spec_names(jcl31)   = "cl31"
    short_spec_names(jcl32)   = "cl32"
    short_spec_names(jcl33)   = "cl33"
    short_spec_names(jcl34)   = "cl34"
    short_spec_names(jcl35)   = "cl35"
    short_spec_names(jcl36)   = "cl36"
    short_spec_names(jcl37)   = "cl37"
    short_spec_names(jar31)   = "ar31"
    short_spec_names(jar32)   = "ar32"
    short_spec_names(jar33)   = "ar33"
    short_spec_names(jar34)   = "ar34"
    short_spec_names(jar35)   = "ar35"
    short_spec_names(jar36)   = "ar36"
    short_spec_names(jar37)   = "ar37"
    short_spec_names(jar38)   = "ar38"
    short_spec_names(jk33)   = "k33"
    short_spec_names(jk34)   = "k34"
    short_spec_names(jk35)   = "k35"
    short_spec_names(jk36)   = "k36"
    short_spec_names(jk37)   = "k37"
    short_spec_names(jk38)   = "k38"
    short_spec_names(jk39)   = "k39"
    short_spec_names(jk40)   = "k40"
    short_spec_names(jk41)   = "k41"
    short_spec_names(jca34)   = "ca34"
    short_spec_names(jca35)   = "ca35"
    short_spec_names(jca36)   = "ca36"
    short_spec_names(jca37)   = "ca37"
    short_spec_names(jca38)   = "ca38"
    short_spec_names(jca39)   = "ca39"
    short_spec_names(jca40)   = "ca40"
    short_spec_names(jca41)   = "ca41"
    short_spec_names(jca42)   = "ca42"
    short_spec_names(jca43)   = "ca43"
    short_spec_names(jca44)   = "ca44"
    short_spec_names(jsc36)   = "sc36"
    short_spec_names(jsc37)   = "sc37"
    short_spec_names(jsc38)   = "sc38"
    short_spec_names(jsc39)   = "sc39"
    short_spec_names(jsc40)   = "sc40"
    short_spec_names(jsc41)   = "sc41"
    short_spec_names(jsc42)   = "sc42"
    short_spec_names(jsc43)   = "sc43"
    short_spec_names(jsc44)   = "sc44"
    short_spec_names(jsc45)   = "sc45"
    short_spec_names(jti38)   = "ti38"
    short_spec_names(jti39)   = "ti39"
    short_spec_names(jti40)   = "ti40"
    short_spec_names(jti41)   = "ti41"
    short_spec_names(jti42)   = "ti42"
    short_spec_names(jti43)   = "ti43"
    short_spec_names(jti44)   = "ti44"
    short_spec_names(jti45)   = "ti45"
    short_spec_names(jti46)   = "ti46"
    short_spec_names(jti47)   = "ti47"
    short_spec_names(jti48)   = "ti48"
    short_spec_names(jv40)   = "v40"
    short_spec_names(jv41)   = "v41"
    short_spec_names(jv42)   = "v42"
    short_spec_names(jv43)   = "v43"
    short_spec_names(jv44)   = "v44"
    short_spec_names(jv45)   = "v45"
    short_spec_names(jv46)   = "v46"
    short_spec_names(jv47)   = "v47"
    short_spec_names(jv48)   = "v48"
    short_spec_names(jv49)   = "v49"
    short_spec_names(jcr42)   = "cr42"
    short_spec_names(jcr43)   = "cr43"
    short_spec_names(jcr44)   = "cr44"
    short_spec_names(jcr45)   = "cr45"
    short_spec_names(jcr46)   = "cr46"
    short_spec_names(jcr47)   = "cr47"
    short_spec_names(jcr48)   = "cr48"
    short_spec_names(jcr49)   = "cr49"
    short_spec_names(jcr50)   = "cr50"
    short_spec_names(jcr51)   = "cr51"
    short_spec_names(jcr52)   = "cr52"
    short_spec_names(jmn44)   = "mn44"
    short_spec_names(jmn45)   = "mn45"
    short_spec_names(jmn46)   = "mn46"
    short_spec_names(jmn47)   = "mn47"
    short_spec_names(jmn48)   = "mn48"
    short_spec_names(jmn49)   = "mn49"
    short_spec_names(jmn50)   = "mn50"
    short_spec_names(jmn51)   = "mn51"
    short_spec_names(jmn52)   = "mn52"
    short_spec_names(jmn53)   = "mn53"
    short_spec_names(jmn55)   = "mn55"
    short_spec_names(jfe45)   = "fe45"
    short_spec_names(jfe46)   = "fe46"
    short_spec_names(jfe47)   = "fe47"
    short_spec_names(jfe48)   = "fe48"
    short_spec_names(jfe49)   = "fe49"
    short_spec_names(jfe50)   = "fe50"
    short_spec_names(jfe51)   = "fe51"
    short_spec_names(jfe52)   = "fe52"
    short_spec_names(jfe53)   = "fe53"
    short_spec_names(jfe54)   = "fe54"
    short_spec_names(jfe55)   = "fe55"
    short_spec_names(jfe56)   = "fe56"
    short_spec_names(jco47)   = "co47"
    short_spec_names(jco48)   = "co48"
    short_spec_names(jco49)   = "co49"
    short_spec_names(jco50)   = "co50"
    short_spec_names(jco51)   = "co51"
    short_spec_names(jco52)   = "co52"
    short_spec_names(jco53)   = "co53"
    short_spec_names(jco54)   = "co54"
    short_spec_names(jco55)   = "co55"
    short_spec_names(jco56)   = "co56"
    short_spec_names(jni48)   = "ni48"
    short_spec_names(jni49)   = "ni49"
    short_spec_names(jni50)   = "ni50"
    short_spec_names(jni51)   = "ni51"
    short_spec_names(jni52)   = "ni52"
    short_spec_names(jni53)   = "ni53"
    short_spec_names(jni54)   = "ni54"
    short_spec_names(jni55)   = "ni55"
    short_spec_names(jni56)   = "ni56"

    ebind_per_nucleon(jp)   = 0.00000000000000d+00
    ebind_per_nucleon(jd)   = 1.11228300000000d+00
    ebind_per_nucleon(jhe3)   = 2.57268000000000d+00
    ebind_per_nucleon(jhe4)   = 7.07391500000000d+00
    ebind_per_nucleon(jli7)   = 5.60643900000000d+00
    ebind_per_nucleon(jbe7)   = 5.37154800000000d+00
    ebind_per_nucleon(jbe8)   = 7.06243500000000d+00
    ebind_per_nucleon(jb8)   = 4.71715500000000d+00
    ebind_per_nucleon(jc12)   = 7.68014400000000d+00
    ebind_per_nucleon(jc13)   = 7.46984900000000d+00
    ebind_per_nucleon(jn13)   = 7.23886300000000d+00
    ebind_per_nucleon(jn14)   = 7.47561400000000d+00
    ebind_per_nucleon(jn15)   = 7.69946000000000d+00
    ebind_per_nucleon(jo14)   = 7.05227800000000d+00
    ebind_per_nucleon(jo15)   = 7.46369200000000d+00
    ebind_per_nucleon(jo16)   = 7.97620600000000d+00
    ebind_per_nucleon(jo17)   = 7.75072800000000d+00
    ebind_per_nucleon(jo18)   = 7.76709700000000d+00
    ebind_per_nucleon(jf16)   = 6.96373100000000d+00
    ebind_per_nucleon(jf17)   = 7.54232800000000d+00
    ebind_per_nucleon(jf18)   = 7.63163800000000d+00
    ebind_per_nucleon(jf19)   = 7.77901800000000d+00
    ebind_per_nucleon(jf20)   = 7.72013400000000d+00
    ebind_per_nucleon(jne18)   = 7.34125700000000d+00
    ebind_per_nucleon(jne19)   = 7.56734300000000d+00
    ebind_per_nucleon(jne20)   = 8.03224000000000d+00
    ebind_per_nucleon(jne21)   = 7.97171300000000d+00
    ebind_per_nucleon(jne22)   = 8.08046500000000d+00
    ebind_per_nucleon(jna19)   = 6.93788500000000d+00
    ebind_per_nucleon(jna20)   = 7.29849600000000d+00
    ebind_per_nucleon(jna21)   = 7.76554700000000d+00
    ebind_per_nucleon(jna22)   = 7.91566700000000d+00
    ebind_per_nucleon(jna23)   = 8.11149300000000d+00
    ebind_per_nucleon(jna24)   = 8.06348800000000d+00
    ebind_per_nucleon(jmg21)   = 7.10503100000000d+00
    ebind_per_nucleon(jmg22)   = 7.66276100000000d+00
    ebind_per_nucleon(jmg23)   = 7.90111500000000d+00
    ebind_per_nucleon(jmg24)   = 8.26070900000000d+00
    ebind_per_nucleon(jmg25)   = 8.22350200000000d+00
    ebind_per_nucleon(jmg26)   = 8.33387000000000d+00
    ebind_per_nucleon(jal22)   = 6.78200000000000d+00
    ebind_per_nucleon(jal23)   = 7.33572700000000d+00
    ebind_per_nucleon(jal24)   = 7.64958200000000d+00
    ebind_per_nucleon(jal25)   = 8.02113600000000d+00
    ebind_per_nucleon(jal26)   = 8.14976500000000d+00
    ebind_per_nucleon(jal27)   = 8.33155300000000d+00
    ebind_per_nucleon(jal28)   = 8.30989400000000d+00
    ebind_per_nucleon(jsi24)   = 7.16723200000000d+00
    ebind_per_nucleon(jsi25)   = 7.48011000000000d+00
    ebind_per_nucleon(jsi26)   = 7.92470800000000d+00
    ebind_per_nucleon(jsi27)   = 8.12434100000000d+00
    ebind_per_nucleon(jsi28)   = 8.44774400000000d+00
    ebind_per_nucleon(jsi29)   = 8.44863500000000d+00
    ebind_per_nucleon(jsi30)   = 8.52065400000000d+00
    ebind_per_nucleon(jp26)   = 7.19800000000000d+00
    ebind_per_nucleon(jp27)   = 7.66343800000000d+00
    ebind_per_nucleon(jp28)   = 7.90747900000000d+00
    ebind_per_nucleon(jp29)   = 8.25123600000000d+00
    ebind_per_nucleon(jp30)   = 8.35350600000000d+00
    ebind_per_nucleon(jp31)   = 8.48116700000000d+00
    ebind_per_nucleon(jp32)   = 8.46412000000000d+00
    ebind_per_nucleon(js28)   = 7.47879000000000d+00
    ebind_per_nucleon(js29)   = 7.74852000000000d+00
    ebind_per_nucleon(js30)   = 8.12270700000000d+00
    ebind_per_nucleon(js31)   = 8.28180000000000d+00
    ebind_per_nucleon(js32)   = 8.49312900000000d+00
    ebind_per_nucleon(js33)   = 8.49763000000000d+00
    ebind_per_nucleon(js34)   = 8.58349800000000d+00
    ebind_per_nucleon(jcl29)   = 7.15883200000000d+00
    ebind_per_nucleon(jcl30)   = 7.48000000000000d+00
    ebind_per_nucleon(jcl31)   = 7.86920900000000d+00
    ebind_per_nucleon(jcl32)   = 8.07240400000000d+00
    ebind_per_nucleon(jcl33)   = 8.30475500000000d+00
    ebind_per_nucleon(jcl34)   = 8.39897000000000d+00
    ebind_per_nucleon(jcl35)   = 8.52027800000000d+00
    ebind_per_nucleon(jcl36)   = 8.52193100000000d+00
    ebind_per_nucleon(jcl37)   = 8.57028100000000d+00
    ebind_per_nucleon(jar31)   = 7.25200000000000d+00
    ebind_per_nucleon(jar32)   = 7.70000800000000d+00
    ebind_per_nucleon(jar33)   = 7.92895500000000d+00
    ebind_per_nucleon(jar34)   = 8.19767200000000d+00
    ebind_per_nucleon(jar35)   = 8.32746100000000d+00
    ebind_per_nucleon(jar36)   = 8.51990900000000d+00
    ebind_per_nucleon(jar37)   = 8.52713900000000d+00
    ebind_per_nucleon(jar38)   = 8.61428000000000d+00
    ebind_per_nucleon(jk33)   = 7.40700000000000d+00
    ebind_per_nucleon(jk34)   = 7.67000000000000d+00
    ebind_per_nucleon(jk35)   = 7.96584000000000d+00
    ebind_per_nucleon(jk36)   = 8.14221900000000d+00
    ebind_per_nucleon(jk37)   = 8.33984700000000d+00
    ebind_per_nucleon(jk38)   = 8.43805800000000d+00
    ebind_per_nucleon(jk39)   = 8.55702500000000d+00
    ebind_per_nucleon(jk40)   = 8.53809000000000d+00
    ebind_per_nucleon(jk41)   = 8.57607200000000d+00
    ebind_per_nucleon(jca34)   = 7.20400000000000d+00
    ebind_per_nucleon(jca35)   = 7.48700000000000d+00
    ebind_per_nucleon(jca36)   = 7.81587900000000d+00
    ebind_per_nucleon(jca37)   = 8.00345600000000d+00
    ebind_per_nucleon(jca38)   = 8.24004300000000d+00
    ebind_per_nucleon(jca39)   = 8.36967000000000d+00
    ebind_per_nucleon(jca40)   = 8.55130300000000d+00
    ebind_per_nucleon(jca41)   = 8.54670600000000d+00
    ebind_per_nucleon(jca42)   = 8.61656300000000d+00
    ebind_per_nucleon(jca43)   = 8.60066300000000d+00
    ebind_per_nucleon(jca44)   = 8.65817500000000d+00
    ebind_per_nucleon(jsc36)   = 7.18900000000000d+00
    ebind_per_nucleon(jsc37)   = 7.53200000000000d+00
    ebind_per_nucleon(jsc38)   = 7.75100000000000d+00
    ebind_per_nucleon(jsc39)   = 8.01345600000000d+00
    ebind_per_nucleon(jsc40)   = 8.17366900000000d+00
    ebind_per_nucleon(jsc41)   = 8.36919800000000d+00
    ebind_per_nucleon(jsc42)   = 8.44493300000000d+00
    ebind_per_nucleon(jsc43)   = 8.53082500000000d+00
    ebind_per_nucleon(jsc44)   = 8.55737900000000d+00
    ebind_per_nucleon(jsc45)   = 8.61893100000000d+00
    ebind_per_nucleon(jti38)   = 7.33200000000000d+00
    ebind_per_nucleon(jti39)   = 7.57400000000000d+00
    ebind_per_nucleon(jti40)   = 7.86228600000000d+00
    ebind_per_nucleon(jti41)   = 8.03438800000000d+00
    ebind_per_nucleon(jti42)   = 8.25924700000000d+00
    ebind_per_nucleon(jti43)   = 8.35293200000000d+00
    ebind_per_nucleon(jti44)   = 8.53352000000000d+00
    ebind_per_nucleon(jti45)   = 8.55572200000000d+00
    ebind_per_nucleon(jti46)   = 8.65645100000000d+00
    ebind_per_nucleon(jti47)   = 8.66122700000000d+00
    ebind_per_nucleon(jti48)   = 8.72300600000000d+00
    ebind_per_nucleon(jv40)   = 7.31700000000000d+00
    ebind_per_nucleon(jv41)   = 7.62500000000000d+00
    ebind_per_nucleon(jv42)   = 7.82400000000000d+00
    ebind_per_nucleon(jv43)   = 8.06951200000000d+00
    ebind_per_nucleon(jv44)   = 8.21046300000000d+00
    ebind_per_nucleon(jv45)   = 8.38002900000000d+00
    ebind_per_nucleon(jv46)   = 8.48613000000000d+00
    ebind_per_nucleon(jv47)   = 8.58222500000000d+00
    ebind_per_nucleon(jv48)   = 8.62306100000000d+00
    ebind_per_nucleon(jv49)   = 8.68290800000000d+00
    ebind_per_nucleon(jcr42)   = 7.46400000000000d+00
    ebind_per_nucleon(jcr43)   = 7.68000000000000d+00
    ebind_per_nucleon(jcr44)   = 7.94800000000000d+00
    ebind_per_nucleon(jcr45)   = 8.08772800000000d+00
    ebind_per_nucleon(jcr46)   = 8.30382300000000d+00
    ebind_per_nucleon(jcr47)   = 8.40719500000000d+00
    ebind_per_nucleon(jcr48)   = 8.57226900000000d+00
    ebind_per_nucleon(jcr49)   = 8.61329100000000d+00
    ebind_per_nucleon(jcr50)   = 8.70103200000000d+00
    ebind_per_nucleon(jcr51)   = 8.71200500000000d+00
    ebind_per_nucleon(jcr52)   = 8.77598900000000d+00
    ebind_per_nucleon(jmn44)   = 7.46700000000000d+00
    ebind_per_nucleon(jmn45)   = 7.75300000000000d+00
    ebind_per_nucleon(jmn46)   = 7.91900000000000d+00
    ebind_per_nucleon(jmn47)   = 8.13531100000000d+00
    ebind_per_nucleon(jmn48)   = 8.27418500000000d+00
    ebind_per_nucleon(jmn49)   = 8.43992900000000d+00
    ebind_per_nucleon(jmn50)   = 8.53269600000000d+00
    ebind_per_nucleon(jmn51)   = 8.63377200000000d+00
    ebind_per_nucleon(jmn52)   = 8.67032900000000d+00
    ebind_per_nucleon(jmn53)   = 8.73417500000000d+00
    ebind_per_nucleon(jmn55)   = 8.76502200000000d+00
    ebind_per_nucleon(jfe45)   = 7.31300000000000d+00
    ebind_per_nucleon(jfe46)   = 7.60900000000000d+00
    ebind_per_nucleon(jfe47)   = 7.78500000000000d+00
    ebind_per_nucleon(jfe48)   = 8.02300000000000d+00
    ebind_per_nucleon(jfe49)   = 8.16131100000000d+00
    ebind_per_nucleon(jfe50)   = 8.35402600000000d+00
    ebind_per_nucleon(jfe51)   = 8.46075900000000d+00
    ebind_per_nucleon(jfe52)   = 8.60957400000000d+00
    ebind_per_nucleon(jfe53)   = 8.64879900000000d+00
    ebind_per_nucleon(jfe54)   = 8.73638200000000d+00
    ebind_per_nucleon(jfe55)   = 8.74659500000000d+00
    ebind_per_nucleon(jfe56)   = 8.79035400000000d+00
    ebind_per_nucleon(jco47)   = 7.40100000000000d+00
    ebind_per_nucleon(jco48)   = 7.60000000000000d+00
    ebind_per_nucleon(jco49)   = 7.84200000000000d+00
    ebind_per_nucleon(jco50)   = 8.00100000000000d+00
    ebind_per_nucleon(jco51)   = 8.19325400000000d+00
    ebind_per_nucleon(jco52)   = 8.32588600000000d+00
    ebind_per_nucleon(jco53)   = 8.47765800000000d+00
    ebind_per_nucleon(jco54)   = 8.56921700000000d+00
    ebind_per_nucleon(jco55)   = 8.66961800000000d+00
    ebind_per_nucleon(jco56)   = 8.69483600000000d+00
    ebind_per_nucleon(jni48)   = 7.26500000000000d+00
    ebind_per_nucleon(jni49)   = 7.45700000000000d+00
    ebind_per_nucleon(jni50)   = 7.71600000000000d+00
    ebind_per_nucleon(jni51)   = 7.87500000000000d+00
    ebind_per_nucleon(jni52)   = 8.07900000000000d+00
    ebind_per_nucleon(jni53)   = 8.21707400000000d+00
    ebind_per_nucleon(jni54)   = 8.39303200000000d+00
    ebind_per_nucleon(jni55)   = 8.49732000000000d+00
    ebind_per_nucleon(jni56)   = 8.64277900000000d+00

    aion(jp)   = 1.00000000000000d+00
    aion(jd)   = 2.00000000000000d+00
    aion(jhe3)   = 3.00000000000000d+00
    aion(jhe4)   = 4.00000000000000d+00
    aion(jli7)   = 7.00000000000000d+00
    aion(jbe7)   = 7.00000000000000d+00
    aion(jbe8)   = 8.00000000000000d+00
    aion(jb8)   = 8.00000000000000d+00
    aion(jc12)   = 1.20000000000000d+01
    aion(jc13)   = 1.30000000000000d+01
    aion(jn13)   = 1.30000000000000d+01
    aion(jn14)   = 1.40000000000000d+01
    aion(jn15)   = 1.50000000000000d+01
    aion(jo14)   = 1.40000000000000d+01
    aion(jo15)   = 1.50000000000000d+01
    aion(jo16)   = 1.60000000000000d+01
    aion(jo17)   = 1.70000000000000d+01
    aion(jo18)   = 1.80000000000000d+01
    aion(jf16)   = 1.60000000000000d+01
    aion(jf17)   = 1.70000000000000d+01
    aion(jf18)   = 1.80000000000000d+01
    aion(jf19)   = 1.90000000000000d+01
    aion(jf20)   = 2.00000000000000d+01
    aion(jne18)   = 1.80000000000000d+01
    aion(jne19)   = 1.90000000000000d+01
    aion(jne20)   = 2.00000000000000d+01
    aion(jne21)   = 2.10000000000000d+01
    aion(jne22)   = 2.20000000000000d+01
    aion(jna19)   = 1.90000000000000d+01
    aion(jna20)   = 2.00000000000000d+01
    aion(jna21)   = 2.10000000000000d+01
    aion(jna22)   = 2.20000000000000d+01
    aion(jna23)   = 2.30000000000000d+01
    aion(jna24)   = 2.40000000000000d+01
    aion(jmg21)   = 2.10000000000000d+01
    aion(jmg22)   = 2.20000000000000d+01
    aion(jmg23)   = 2.30000000000000d+01
    aion(jmg24)   = 2.40000000000000d+01
    aion(jmg25)   = 2.50000000000000d+01
    aion(jmg26)   = 2.60000000000000d+01
    aion(jal22)   = 2.20000000000000d+01
    aion(jal23)   = 2.30000000000000d+01
    aion(jal24)   = 2.40000000000000d+01
    aion(jal25)   = 2.50000000000000d+01
    aion(jal26)   = 2.60000000000000d+01
    aion(jal27)   = 2.70000000000000d+01
    aion(jal28)   = 2.80000000000000d+01
    aion(jsi24)   = 2.40000000000000d+01
    aion(jsi25)   = 2.50000000000000d+01
    aion(jsi26)   = 2.60000000000000d+01
    aion(jsi27)   = 2.70000000000000d+01
    aion(jsi28)   = 2.80000000000000d+01
    aion(jsi29)   = 2.90000000000000d+01
    aion(jsi30)   = 3.00000000000000d+01
    aion(jp26)   = 2.60000000000000d+01
    aion(jp27)   = 2.70000000000000d+01
    aion(jp28)   = 2.80000000000000d+01
    aion(jp29)   = 2.90000000000000d+01
    aion(jp30)   = 3.00000000000000d+01
    aion(jp31)   = 3.10000000000000d+01
    aion(jp32)   = 3.20000000000000d+01
    aion(js28)   = 2.80000000000000d+01
    aion(js29)   = 2.90000000000000d+01
    aion(js30)   = 3.00000000000000d+01
    aion(js31)   = 3.10000000000000d+01
    aion(js32)   = 3.20000000000000d+01
    aion(js33)   = 3.30000000000000d+01
    aion(js34)   = 3.40000000000000d+01
    aion(jcl29)   = 2.90000000000000d+01
    aion(jcl30)   = 3.00000000000000d+01
    aion(jcl31)   = 3.10000000000000d+01
    aion(jcl32)   = 3.20000000000000d+01
    aion(jcl33)   = 3.30000000000000d+01
    aion(jcl34)   = 3.40000000000000d+01
    aion(jcl35)   = 3.50000000000000d+01
    aion(jcl36)   = 3.60000000000000d+01
    aion(jcl37)   = 3.70000000000000d+01
    aion(jar31)   = 3.10000000000000d+01
    aion(jar32)   = 3.20000000000000d+01
    aion(jar33)   = 3.30000000000000d+01
    aion(jar34)   = 3.40000000000000d+01
    aion(jar35)   = 3.50000000000000d+01
    aion(jar36)   = 3.60000000000000d+01
    aion(jar37)   = 3.70000000000000d+01
    aion(jar38)   = 3.80000000000000d+01
    aion(jk33)   = 3.30000000000000d+01
    aion(jk34)   = 3.40000000000000d+01
    aion(jk35)   = 3.50000000000000d+01
    aion(jk36)   = 3.60000000000000d+01
    aion(jk37)   = 3.70000000000000d+01
    aion(jk38)   = 3.80000000000000d+01
    aion(jk39)   = 3.90000000000000d+01
    aion(jk40)   = 4.00000000000000d+01
    aion(jk41)   = 4.10000000000000d+01
    aion(jca34)   = 3.40000000000000d+01
    aion(jca35)   = 3.50000000000000d+01
    aion(jca36)   = 3.60000000000000d+01
    aion(jca37)   = 3.70000000000000d+01
    aion(jca38)   = 3.80000000000000d+01
    aion(jca39)   = 3.90000000000000d+01
    aion(jca40)   = 4.00000000000000d+01
    aion(jca41)   = 4.10000000000000d+01
    aion(jca42)   = 4.20000000000000d+01
    aion(jca43)   = 4.30000000000000d+01
    aion(jca44)   = 4.40000000000000d+01
    aion(jsc36)   = 3.60000000000000d+01
    aion(jsc37)   = 3.70000000000000d+01
    aion(jsc38)   = 3.80000000000000d+01
    aion(jsc39)   = 3.90000000000000d+01
    aion(jsc40)   = 4.00000000000000d+01
    aion(jsc41)   = 4.10000000000000d+01
    aion(jsc42)   = 4.20000000000000d+01
    aion(jsc43)   = 4.30000000000000d+01
    aion(jsc44)   = 4.40000000000000d+01
    aion(jsc45)   = 4.50000000000000d+01
    aion(jti38)   = 3.80000000000000d+01
    aion(jti39)   = 3.90000000000000d+01
    aion(jti40)   = 4.00000000000000d+01
    aion(jti41)   = 4.10000000000000d+01
    aion(jti42)   = 4.20000000000000d+01
    aion(jti43)   = 4.30000000000000d+01
    aion(jti44)   = 4.40000000000000d+01
    aion(jti45)   = 4.50000000000000d+01
    aion(jti46)   = 4.60000000000000d+01
    aion(jti47)   = 4.70000000000000d+01
    aion(jti48)   = 4.80000000000000d+01
    aion(jv40)   = 4.00000000000000d+01
    aion(jv41)   = 4.10000000000000d+01
    aion(jv42)   = 4.20000000000000d+01
    aion(jv43)   = 4.30000000000000d+01
    aion(jv44)   = 4.40000000000000d+01
    aion(jv45)   = 4.50000000000000d+01
    aion(jv46)   = 4.60000000000000d+01
    aion(jv47)   = 4.70000000000000d+01
    aion(jv48)   = 4.80000000000000d+01
    aion(jv49)   = 4.90000000000000d+01
    aion(jcr42)   = 4.20000000000000d+01
    aion(jcr43)   = 4.30000000000000d+01
    aion(jcr44)   = 4.40000000000000d+01
    aion(jcr45)   = 4.50000000000000d+01
    aion(jcr46)   = 4.60000000000000d+01
    aion(jcr47)   = 4.70000000000000d+01
    aion(jcr48)   = 4.80000000000000d+01
    aion(jcr49)   = 4.90000000000000d+01
    aion(jcr50)   = 5.00000000000000d+01
    aion(jcr51)   = 5.10000000000000d+01
    aion(jcr52)   = 5.20000000000000d+01
    aion(jmn44)   = 4.40000000000000d+01
    aion(jmn45)   = 4.50000000000000d+01
    aion(jmn46)   = 4.60000000000000d+01
    aion(jmn47)   = 4.70000000000000d+01
    aion(jmn48)   = 4.80000000000000d+01
    aion(jmn49)   = 4.90000000000000d+01
    aion(jmn50)   = 5.00000000000000d+01
    aion(jmn51)   = 5.10000000000000d+01
    aion(jmn52)   = 5.20000000000000d+01
    aion(jmn53)   = 5.30000000000000d+01
    aion(jmn55)   = 5.50000000000000d+01
    aion(jfe45)   = 4.50000000000000d+01
    aion(jfe46)   = 4.60000000000000d+01
    aion(jfe47)   = 4.70000000000000d+01
    aion(jfe48)   = 4.80000000000000d+01
    aion(jfe49)   = 4.90000000000000d+01
    aion(jfe50)   = 5.00000000000000d+01
    aion(jfe51)   = 5.10000000000000d+01
    aion(jfe52)   = 5.20000000000000d+01
    aion(jfe53)   = 5.30000000000000d+01
    aion(jfe54)   = 5.40000000000000d+01
    aion(jfe55)   = 5.50000000000000d+01
    aion(jfe56)   = 5.60000000000000d+01
    aion(jco47)   = 4.70000000000000d+01
    aion(jco48)   = 4.80000000000000d+01
    aion(jco49)   = 4.90000000000000d+01
    aion(jco50)   = 5.00000000000000d+01
    aion(jco51)   = 5.10000000000000d+01
    aion(jco52)   = 5.20000000000000d+01
    aion(jco53)   = 5.30000000000000d+01
    aion(jco54)   = 5.40000000000000d+01
    aion(jco55)   = 5.50000000000000d+01
    aion(jco56)   = 5.60000000000000d+01
    aion(jni48)   = 4.80000000000000d+01
    aion(jni49)   = 4.90000000000000d+01
    aion(jni50)   = 5.00000000000000d+01
    aion(jni51)   = 5.10000000000000d+01
    aion(jni52)   = 5.20000000000000d+01
    aion(jni53)   = 5.30000000000000d+01
    aion(jni54)   = 5.40000000000000d+01
    aion(jni55)   = 5.50000000000000d+01
    aion(jni56)   = 5.60000000000000d+01

    zion(jp)   = 1.00000000000000d+00
    zion(jd)   = 1.00000000000000d+00
    zion(jhe3)   = 2.00000000000000d+00
    zion(jhe4)   = 2.00000000000000d+00
    zion(jli7)   = 3.00000000000000d+00
    zion(jbe7)   = 4.00000000000000d+00
    zion(jbe8)   = 4.00000000000000d+00
    zion(jb8)   = 5.00000000000000d+00
    zion(jc12)   = 6.00000000000000d+00
    zion(jc13)   = 6.00000000000000d+00
    zion(jn13)   = 7.00000000000000d+00
    zion(jn14)   = 7.00000000000000d+00
    zion(jn15)   = 7.00000000000000d+00
    zion(jo14)   = 8.00000000000000d+00
    zion(jo15)   = 8.00000000000000d+00
    zion(jo16)   = 8.00000000000000d+00
    zion(jo17)   = 8.00000000000000d+00
    zion(jo18)   = 8.00000000000000d+00
    zion(jf16)   = 9.00000000000000d+00
    zion(jf17)   = 9.00000000000000d+00
    zion(jf18)   = 9.00000000000000d+00
    zion(jf19)   = 9.00000000000000d+00
    zion(jf20)   = 9.00000000000000d+00
    zion(jne18)   = 1.00000000000000d+01
    zion(jne19)   = 1.00000000000000d+01
    zion(jne20)   = 1.00000000000000d+01
    zion(jne21)   = 1.00000000000000d+01
    zion(jne22)   = 1.00000000000000d+01
    zion(jna19)   = 1.10000000000000d+01
    zion(jna20)   = 1.10000000000000d+01
    zion(jna21)   = 1.10000000000000d+01
    zion(jna22)   = 1.10000000000000d+01
    zion(jna23)   = 1.10000000000000d+01
    zion(jna24)   = 1.10000000000000d+01
    zion(jmg21)   = 1.20000000000000d+01
    zion(jmg22)   = 1.20000000000000d+01
    zion(jmg23)   = 1.20000000000000d+01
    zion(jmg24)   = 1.20000000000000d+01
    zion(jmg25)   = 1.20000000000000d+01
    zion(jmg26)   = 1.20000000000000d+01
    zion(jal22)   = 1.30000000000000d+01
    zion(jal23)   = 1.30000000000000d+01
    zion(jal24)   = 1.30000000000000d+01
    zion(jal25)   = 1.30000000000000d+01
    zion(jal26)   = 1.30000000000000d+01
    zion(jal27)   = 1.30000000000000d+01
    zion(jal28)   = 1.30000000000000d+01
    zion(jsi24)   = 1.40000000000000d+01
    zion(jsi25)   = 1.40000000000000d+01
    zion(jsi26)   = 1.40000000000000d+01
    zion(jsi27)   = 1.40000000000000d+01
    zion(jsi28)   = 1.40000000000000d+01
    zion(jsi29)   = 1.40000000000000d+01
    zion(jsi30)   = 1.40000000000000d+01
    zion(jp26)   = 1.50000000000000d+01
    zion(jp27)   = 1.50000000000000d+01
    zion(jp28)   = 1.50000000000000d+01
    zion(jp29)   = 1.50000000000000d+01
    zion(jp30)   = 1.50000000000000d+01
    zion(jp31)   = 1.50000000000000d+01
    zion(jp32)   = 1.50000000000000d+01
    zion(js28)   = 1.60000000000000d+01
    zion(js29)   = 1.60000000000000d+01
    zion(js30)   = 1.60000000000000d+01
    zion(js31)   = 1.60000000000000d+01
    zion(js32)   = 1.60000000000000d+01
    zion(js33)   = 1.60000000000000d+01
    zion(js34)   = 1.60000000000000d+01
    zion(jcl29)   = 1.70000000000000d+01
    zion(jcl30)   = 1.70000000000000d+01
    zion(jcl31)   = 1.70000000000000d+01
    zion(jcl32)   = 1.70000000000000d+01
    zion(jcl33)   = 1.70000000000000d+01
    zion(jcl34)   = 1.70000000000000d+01
    zion(jcl35)   = 1.70000000000000d+01
    zion(jcl36)   = 1.70000000000000d+01
    zion(jcl37)   = 1.70000000000000d+01
    zion(jar31)   = 1.80000000000000d+01
    zion(jar32)   = 1.80000000000000d+01
    zion(jar33)   = 1.80000000000000d+01
    zion(jar34)   = 1.80000000000000d+01
    zion(jar35)   = 1.80000000000000d+01
    zion(jar36)   = 1.80000000000000d+01
    zion(jar37)   = 1.80000000000000d+01
    zion(jar38)   = 1.80000000000000d+01
    zion(jk33)   = 1.90000000000000d+01
    zion(jk34)   = 1.90000000000000d+01
    zion(jk35)   = 1.90000000000000d+01
    zion(jk36)   = 1.90000000000000d+01
    zion(jk37)   = 1.90000000000000d+01
    zion(jk38)   = 1.90000000000000d+01
    zion(jk39)   = 1.90000000000000d+01
    zion(jk40)   = 1.90000000000000d+01
    zion(jk41)   = 1.90000000000000d+01
    zion(jca34)   = 2.00000000000000d+01
    zion(jca35)   = 2.00000000000000d+01
    zion(jca36)   = 2.00000000000000d+01
    zion(jca37)   = 2.00000000000000d+01
    zion(jca38)   = 2.00000000000000d+01
    zion(jca39)   = 2.00000000000000d+01
    zion(jca40)   = 2.00000000000000d+01
    zion(jca41)   = 2.00000000000000d+01
    zion(jca42)   = 2.00000000000000d+01
    zion(jca43)   = 2.00000000000000d+01
    zion(jca44)   = 2.00000000000000d+01
    zion(jsc36)   = 2.10000000000000d+01
    zion(jsc37)   = 2.10000000000000d+01
    zion(jsc38)   = 2.10000000000000d+01
    zion(jsc39)   = 2.10000000000000d+01
    zion(jsc40)   = 2.10000000000000d+01
    zion(jsc41)   = 2.10000000000000d+01
    zion(jsc42)   = 2.10000000000000d+01
    zion(jsc43)   = 2.10000000000000d+01
    zion(jsc44)   = 2.10000000000000d+01
    zion(jsc45)   = 2.10000000000000d+01
    zion(jti38)   = 2.20000000000000d+01
    zion(jti39)   = 2.20000000000000d+01
    zion(jti40)   = 2.20000000000000d+01
    zion(jti41)   = 2.20000000000000d+01
    zion(jti42)   = 2.20000000000000d+01
    zion(jti43)   = 2.20000000000000d+01
    zion(jti44)   = 2.20000000000000d+01
    zion(jti45)   = 2.20000000000000d+01
    zion(jti46)   = 2.20000000000000d+01
    zion(jti47)   = 2.20000000000000d+01
    zion(jti48)   = 2.20000000000000d+01
    zion(jv40)   = 2.30000000000000d+01
    zion(jv41)   = 2.30000000000000d+01
    zion(jv42)   = 2.30000000000000d+01
    zion(jv43)   = 2.30000000000000d+01
    zion(jv44)   = 2.30000000000000d+01
    zion(jv45)   = 2.30000000000000d+01
    zion(jv46)   = 2.30000000000000d+01
    zion(jv47)   = 2.30000000000000d+01
    zion(jv48)   = 2.30000000000000d+01
    zion(jv49)   = 2.30000000000000d+01
    zion(jcr42)   = 2.40000000000000d+01
    zion(jcr43)   = 2.40000000000000d+01
    zion(jcr44)   = 2.40000000000000d+01
    zion(jcr45)   = 2.40000000000000d+01
    zion(jcr46)   = 2.40000000000000d+01
    zion(jcr47)   = 2.40000000000000d+01
    zion(jcr48)   = 2.40000000000000d+01
    zion(jcr49)   = 2.40000000000000d+01
    zion(jcr50)   = 2.40000000000000d+01
    zion(jcr51)   = 2.40000000000000d+01
    zion(jcr52)   = 2.40000000000000d+01
    zion(jmn44)   = 2.50000000000000d+01
    zion(jmn45)   = 2.50000000000000d+01
    zion(jmn46)   = 2.50000000000000d+01
    zion(jmn47)   = 2.50000000000000d+01
    zion(jmn48)   = 2.50000000000000d+01
    zion(jmn49)   = 2.50000000000000d+01
    zion(jmn50)   = 2.50000000000000d+01
    zion(jmn51)   = 2.50000000000000d+01
    zion(jmn52)   = 2.50000000000000d+01
    zion(jmn53)   = 2.50000000000000d+01
    zion(jmn55)   = 2.50000000000000d+01
    zion(jfe45)   = 2.60000000000000d+01
    zion(jfe46)   = 2.60000000000000d+01
    zion(jfe47)   = 2.60000000000000d+01
    zion(jfe48)   = 2.60000000000000d+01
    zion(jfe49)   = 2.60000000000000d+01
    zion(jfe50)   = 2.60000000000000d+01
    zion(jfe51)   = 2.60000000000000d+01
    zion(jfe52)   = 2.60000000000000d+01
    zion(jfe53)   = 2.60000000000000d+01
    zion(jfe54)   = 2.60000000000000d+01
    zion(jfe55)   = 2.60000000000000d+01
    zion(jfe56)   = 2.60000000000000d+01
    zion(jco47)   = 2.70000000000000d+01
    zion(jco48)   = 2.70000000000000d+01
    zion(jco49)   = 2.70000000000000d+01
    zion(jco50)   = 2.70000000000000d+01
    zion(jco51)   = 2.70000000000000d+01
    zion(jco52)   = 2.70000000000000d+01
    zion(jco53)   = 2.70000000000000d+01
    zion(jco54)   = 2.70000000000000d+01
    zion(jco55)   = 2.70000000000000d+01
    zion(jco56)   = 2.70000000000000d+01
    zion(jni48)   = 2.80000000000000d+01
    zion(jni49)   = 2.80000000000000d+01
    zion(jni50)   = 2.80000000000000d+01
    zion(jni51)   = 2.80000000000000d+01
    zion(jni52)   = 2.80000000000000d+01
    zion(jni53)   = 2.80000000000000d+01
    zion(jni54)   = 2.80000000000000d+01
    zion(jni55)   = 2.80000000000000d+01
    zion(jni56)   = 2.80000000000000d+01

    nion(jp)   = 0.00000000000000d+00
    nion(jd)   = 1.00000000000000d+00
    nion(jhe3)   = 1.00000000000000d+00
    nion(jhe4)   = 2.00000000000000d+00
    nion(jli7)   = 4.00000000000000d+00
    nion(jbe7)   = 3.00000000000000d+00
    nion(jbe8)   = 4.00000000000000d+00
    nion(jb8)   = 3.00000000000000d+00
    nion(jc12)   = 6.00000000000000d+00
    nion(jc13)   = 7.00000000000000d+00
    nion(jn13)   = 6.00000000000000d+00
    nion(jn14)   = 7.00000000000000d+00
    nion(jn15)   = 8.00000000000000d+00
    nion(jo14)   = 6.00000000000000d+00
    nion(jo15)   = 7.00000000000000d+00
    nion(jo16)   = 8.00000000000000d+00
    nion(jo17)   = 9.00000000000000d+00
    nion(jo18)   = 1.00000000000000d+01
    nion(jf16)   = 7.00000000000000d+00
    nion(jf17)   = 8.00000000000000d+00
    nion(jf18)   = 9.00000000000000d+00
    nion(jf19)   = 1.00000000000000d+01
    nion(jf20)   = 1.10000000000000d+01
    nion(jne18)   = 8.00000000000000d+00
    nion(jne19)   = 9.00000000000000d+00
    nion(jne20)   = 1.00000000000000d+01
    nion(jne21)   = 1.10000000000000d+01
    nion(jne22)   = 1.20000000000000d+01
    nion(jna19)   = 8.00000000000000d+00
    nion(jna20)   = 9.00000000000000d+00
    nion(jna21)   = 1.00000000000000d+01
    nion(jna22)   = 1.10000000000000d+01
    nion(jna23)   = 1.20000000000000d+01
    nion(jna24)   = 1.30000000000000d+01
    nion(jmg21)   = 9.00000000000000d+00
    nion(jmg22)   = 1.00000000000000d+01
    nion(jmg23)   = 1.10000000000000d+01
    nion(jmg24)   = 1.20000000000000d+01
    nion(jmg25)   = 1.30000000000000d+01
    nion(jmg26)   = 1.40000000000000d+01
    nion(jal22)   = 9.00000000000000d+00
    nion(jal23)   = 1.00000000000000d+01
    nion(jal24)   = 1.10000000000000d+01
    nion(jal25)   = 1.20000000000000d+01
    nion(jal26)   = 1.30000000000000d+01
    nion(jal27)   = 1.40000000000000d+01
    nion(jal28)   = 1.50000000000000d+01
    nion(jsi24)   = 1.00000000000000d+01
    nion(jsi25)   = 1.10000000000000d+01
    nion(jsi26)   = 1.20000000000000d+01
    nion(jsi27)   = 1.30000000000000d+01
    nion(jsi28)   = 1.40000000000000d+01
    nion(jsi29)   = 1.50000000000000d+01
    nion(jsi30)   = 1.60000000000000d+01
    nion(jp26)   = 1.10000000000000d+01
    nion(jp27)   = 1.20000000000000d+01
    nion(jp28)   = 1.30000000000000d+01
    nion(jp29)   = 1.40000000000000d+01
    nion(jp30)   = 1.50000000000000d+01
    nion(jp31)   = 1.60000000000000d+01
    nion(jp32)   = 1.70000000000000d+01
    nion(js28)   = 1.20000000000000d+01
    nion(js29)   = 1.30000000000000d+01
    nion(js30)   = 1.40000000000000d+01
    nion(js31)   = 1.50000000000000d+01
    nion(js32)   = 1.60000000000000d+01
    nion(js33)   = 1.70000000000000d+01
    nion(js34)   = 1.80000000000000d+01
    nion(jcl29)   = 1.20000000000000d+01
    nion(jcl30)   = 1.30000000000000d+01
    nion(jcl31)   = 1.40000000000000d+01
    nion(jcl32)   = 1.50000000000000d+01
    nion(jcl33)   = 1.60000000000000d+01
    nion(jcl34)   = 1.70000000000000d+01
    nion(jcl35)   = 1.80000000000000d+01
    nion(jcl36)   = 1.90000000000000d+01
    nion(jcl37)   = 2.00000000000000d+01
    nion(jar31)   = 1.30000000000000d+01
    nion(jar32)   = 1.40000000000000d+01
    nion(jar33)   = 1.50000000000000d+01
    nion(jar34)   = 1.60000000000000d+01
    nion(jar35)   = 1.70000000000000d+01
    nion(jar36)   = 1.80000000000000d+01
    nion(jar37)   = 1.90000000000000d+01
    nion(jar38)   = 2.00000000000000d+01
    nion(jk33)   = 1.40000000000000d+01
    nion(jk34)   = 1.50000000000000d+01
    nion(jk35)   = 1.60000000000000d+01
    nion(jk36)   = 1.70000000000000d+01
    nion(jk37)   = 1.80000000000000d+01
    nion(jk38)   = 1.90000000000000d+01
    nion(jk39)   = 2.00000000000000d+01
    nion(jk40)   = 2.10000000000000d+01
    nion(jk41)   = 2.20000000000000d+01
    nion(jca34)   = 1.40000000000000d+01
    nion(jca35)   = 1.50000000000000d+01
    nion(jca36)   = 1.60000000000000d+01
    nion(jca37)   = 1.70000000000000d+01
    nion(jca38)   = 1.80000000000000d+01
    nion(jca39)   = 1.90000000000000d+01
    nion(jca40)   = 2.00000000000000d+01
    nion(jca41)   = 2.10000000000000d+01
    nion(jca42)   = 2.20000000000000d+01
    nion(jca43)   = 2.30000000000000d+01
    nion(jca44)   = 2.40000000000000d+01
    nion(jsc36)   = 1.50000000000000d+01
    nion(jsc37)   = 1.60000000000000d+01
    nion(jsc38)   = 1.70000000000000d+01
    nion(jsc39)   = 1.80000000000000d+01
    nion(jsc40)   = 1.90000000000000d+01
    nion(jsc41)   = 2.00000000000000d+01
    nion(jsc42)   = 2.10000000000000d+01
    nion(jsc43)   = 2.20000000000000d+01
    nion(jsc44)   = 2.30000000000000d+01
    nion(jsc45)   = 2.40000000000000d+01
    nion(jti38)   = 1.60000000000000d+01
    nion(jti39)   = 1.70000000000000d+01
    nion(jti40)   = 1.80000000000000d+01
    nion(jti41)   = 1.90000000000000d+01
    nion(jti42)   = 2.00000000000000d+01
    nion(jti43)   = 2.10000000000000d+01
    nion(jti44)   = 2.20000000000000d+01
    nion(jti45)   = 2.30000000000000d+01
    nion(jti46)   = 2.40000000000000d+01
    nion(jti47)   = 2.50000000000000d+01
    nion(jti48)   = 2.60000000000000d+01
    nion(jv40)   = 1.70000000000000d+01
    nion(jv41)   = 1.80000000000000d+01
    nion(jv42)   = 1.90000000000000d+01
    nion(jv43)   = 2.00000000000000d+01
    nion(jv44)   = 2.10000000000000d+01
    nion(jv45)   = 2.20000000000000d+01
    nion(jv46)   = 2.30000000000000d+01
    nion(jv47)   = 2.40000000000000d+01
    nion(jv48)   = 2.50000000000000d+01
    nion(jv49)   = 2.60000000000000d+01
    nion(jcr42)   = 1.80000000000000d+01
    nion(jcr43)   = 1.90000000000000d+01
    nion(jcr44)   = 2.00000000000000d+01
    nion(jcr45)   = 2.10000000000000d+01
    nion(jcr46)   = 2.20000000000000d+01
    nion(jcr47)   = 2.30000000000000d+01
    nion(jcr48)   = 2.40000000000000d+01
    nion(jcr49)   = 2.50000000000000d+01
    nion(jcr50)   = 2.60000000000000d+01
    nion(jcr51)   = 2.70000000000000d+01
    nion(jcr52)   = 2.80000000000000d+01
    nion(jmn44)   = 1.90000000000000d+01
    nion(jmn45)   = 2.00000000000000d+01
    nion(jmn46)   = 2.10000000000000d+01
    nion(jmn47)   = 2.20000000000000d+01
    nion(jmn48)   = 2.30000000000000d+01
    nion(jmn49)   = 2.40000000000000d+01
    nion(jmn50)   = 2.50000000000000d+01
    nion(jmn51)   = 2.60000000000000d+01
    nion(jmn52)   = 2.70000000000000d+01
    nion(jmn53)   = 2.80000000000000d+01
    nion(jmn55)   = 3.00000000000000d+01
    nion(jfe45)   = 1.90000000000000d+01
    nion(jfe46)   = 2.00000000000000d+01
    nion(jfe47)   = 2.10000000000000d+01
    nion(jfe48)   = 2.20000000000000d+01
    nion(jfe49)   = 2.30000000000000d+01
    nion(jfe50)   = 2.40000000000000d+01
    nion(jfe51)   = 2.50000000000000d+01
    nion(jfe52)   = 2.60000000000000d+01
    nion(jfe53)   = 2.70000000000000d+01
    nion(jfe54)   = 2.80000000000000d+01
    nion(jfe55)   = 2.90000000000000d+01
    nion(jfe56)   = 3.00000000000000d+01
    nion(jco47)   = 2.00000000000000d+01
    nion(jco48)   = 2.10000000000000d+01
    nion(jco49)   = 2.20000000000000d+01
    nion(jco50)   = 2.30000000000000d+01
    nion(jco51)   = 2.40000000000000d+01
    nion(jco52)   = 2.50000000000000d+01
    nion(jco53)   = 2.60000000000000d+01
    nion(jco54)   = 2.70000000000000d+01
    nion(jco55)   = 2.80000000000000d+01
    nion(jco56)   = 2.90000000000000d+01
    nion(jni48)   = 2.00000000000000d+01
    nion(jni49)   = 2.10000000000000d+01
    nion(jni50)   = 2.20000000000000d+01
    nion(jni51)   = 2.30000000000000d+01
    nion(jni52)   = 2.40000000000000d+01
    nion(jni53)   = 2.50000000000000d+01
    nion(jni54)   = 2.60000000000000d+01
    nion(jni55)   = 2.70000000000000d+01
    nion(jni56)   = 2.80000000000000d+01

    do i = 1, nspec
       bion(i) = ebind_per_nucleon(i) * aion(i) * ERG_PER_MeV
    end do

    ! Set the mass
    mion(:) = nion(:) * mass_neutron + zion(:) * (mass_proton + mass_electron) &
         - bion(:)/(c_light**2)

    ! Molar mass
    wion(:) = avo * mion(:)

    ! Common approximation
    !wion(:) = aion(:)

    !$acc update device(aion, zion, bion, nion, mion, wion)

#ifdef REACT_SPARSE_JACOBIAN
    ! Set CSR format metadata for Jacobian
    allocate(csr_jac_col_index(NETWORK_SPARSE_JAC_NNZ))
    allocate(csr_jac_row_count(nspec_evolve + 3)) ! neq + 1

    csr_jac_col_index = [ &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      190, &
      1, &
      2, &
      3, &
      4, &
      6, &
      190, &
      1, &
      2, &
      3, &
      4, &
      6, &
      190, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      8, &
      9, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      190, &
      1, &
      4, &
      5, &
      6, &
      190, &
      1, &
      2, &
      3, &
      4, &
      6, &
      8, &
      190, &
      7, &
      8, &
      190, &
      1, &
      6, &
      8, &
      190, &
      1, &
      4, &
      9, &
      11, &
      13, &
      16, &
      26, &
      190, &
      1, &
      10, &
      11, &
      190, &
      1, &
      4, &
      9, &
      11, &
      14, &
      16, &
      190, &
      1, &
      4, &
      10, &
      12, &
      14, &
      15, &
      17, &
      21, &
      190, &
      1, &
      4, &
      9, &
      13, &
      15, &
      16, &
      18, &
      22, &
      190, &
      1, &
      4, &
      11, &
      14, &
      20, &
      24, &
      190, &
      1, &
      4, &
      12, &
      15, &
      21, &
      25, &
      190, &
      1, &
      4, &
      9, &
      11, &
      13, &
      16, &
      20, &
      22, &
      26, &
      190, &
      1, &
      4, &
      12, &
      17, &
      20, &
      21, &
      23, &
      27, &
      190, &
      1, &
      4, &
      13, &
      18, &
      21, &
      22, &
      190, &
      1, &
      4, &
      15, &
      19, &
      190, &
      1, &
      4, &
      14, &
      16, &
      20, &
      24, &
      26, &
      190, &
      1, &
      4, &
      12, &
      15, &
      17, &
      21, &
      24, &
      25, &
      27, &
      190, &
      1, &
      4, &
      13, &
      16, &
      18, &
      22, &
      25, &
      26, &
      190, &
      1, &
      4, &
      17, &
      23, &
      27, &
      190, &
      1, &
      4, &
      14, &
      20, &
      24, &
      190, &
      1, &
      4, &
      15, &
      21, &
      25, &
      29, &
      190, &
      1, &
      4, &
      9, &
      16, &
      20, &
      22, &
      23, &
      26, &
      30, &
      190, &
      1, &
      4, &
      17, &
      21, &
      23, &
      27, &
      31, &
      190, &
      1, &
      4, &
      18, &
      28, &
      32, &
      190, &
      1, &
      4, &
      24, &
      29, &
      190, &
      1, &
      4, &
      19, &
      25, &
      30, &
      190, &
      1, &
      4, &
      20, &
      26, &
      31, &
      35, &
      190, &
      1, &
      4, &
      21, &
      27, &
      32, &
      36, &
      190, &
      1, &
      4, &
      22, &
      28, &
      33, &
      37, &
      190, &
      1, &
      4, &
      23, &
      34, &
      190, &
      1, &
      4, &
      30, &
      35, &
      190, &
      1, &
      4, &
      24, &
      31, &
      36, &
      41, &
      190, &
      1, &
      4, &
      25, &
      32, &
      37, &
      42, &
      190, &
      1, &
      4, &
      26, &
      33, &
      38, &
      43, &
      190, &
      1, &
      4, &
      27, &
      34, &
      39, &
      44, &
      190, &
      1, &
      4, &
      28, &
      40, &
      45, &
      190, &
      1, &
      4, &
      35, &
      41, &
      190, &
      1, &
      4, &
      29, &
      36, &
      42, &
      190, &
      1, &
      4, &
      30, &
      37, &
      43, &
      48, &
      190, &
      1, &
      4, &
      31, &
      38, &
      44, &
      49, &
      190, &
      1, &
      4, &
      32, &
      39, &
      45, &
      50, &
      190, &
      1, &
      4, &
      33, &
      40, &
      46, &
      51, &
      190, &
      1, &
      4, &
      34, &
      47, &
      190, &
      1, &
      4, &
      42, &
      48, &
      190, &
      1, &
      4, &
      35, &
      43, &
      49, &
      190, &
      1, &
      4, &
      36, &
      44, &
      50, &
      55, &
      190, &
      1, &
      4, &
      37, &
      45, &
      51, &
      56, &
      190, &
      1, &
      4, &
      38, &
      46, &
      52, &
      57, &
      190, &
      1, &
      4, &
      39, &
      47, &
      53, &
      58, &
      190, &
      1, &
      4, &
      40, &
      54, &
      59, &
      190, &
      1, &
      4, &
      41, &
      49, &
      55, &
      190, &
      1, &
      4, &
      42, &
      50, &
      56, &
      190, &
      1, &
      4, &
      43, &
      51, &
      57, &
      62, &
      190, &
      1, &
      4, &
      44, &
      52, &
      58, &
      63, &
      190, &
      1, &
      4, &
      45, &
      53, &
      59, &
      64, &
      190, &
      1, &
      4, &
      46, &
      54, &
      60, &
      65, &
      190, &
      1, &
      4, &
      47, &
      61, &
      190, &
      1, &
      4, &
      48, &
      56, &
      62, &
      190, &
      1, &
      4, &
      49, &
      57, &
      63, &
      69, &
      190, &
      1, &
      4, &
      50, &
      58, &
      64, &
      70, &
      190, &
      1, &
      4, &
      51, &
      59, &
      65, &
      71, &
      190, &
      1, &
      4, &
      52, &
      60, &
      66, &
      72, &
      190, &
      1, &
      4, &
      53, &
      61, &
      67, &
      73, &
      190, &
      1, &
      4, &
      54, &
      68, &
      74, &
      190, &
      1, &
      4, &
      62, &
      69, &
      190, &
      1, &
      4, &
      55, &
      63, &
      70, &
      190, &
      1, &
      4, &
      56, &
      64, &
      71, &
      78, &
      190, &
      1, &
      4, &
      57, &
      65, &
      72, &
      79, &
      190, &
      1, &
      4, &
      58, &
      66, &
      73, &
      80, &
      190, &
      1, &
      4, &
      59, &
      67, &
      74, &
      81, &
      190, &
      1, &
      4, &
      60, &
      68, &
      75, &
      82, &
      190, &
      1, &
      4, &
      61, &
      76, &
      190, &
      1, &
      4, &
      77, &
      84, &
      190, &
      1, &
      4, &
      70, &
      78, &
      190, &
      1, &
      4, &
      62, &
      71, &
      79, &
      190, &
      1, &
      4, &
      63, &
      72, &
      80, &
      86, &
      190, &
      1, &
      4, &
      64, &
      73, &
      81, &
      87, &
      190, &
      1, &
      4, &
      65, &
      74, &
      82, &
      88, &
      190, &
      1, &
      4, &
      66, &
      75, &
      83, &
      89, &
      190, &
      1, &
      4, &
      67, &
      76, &
      84, &
      90, &
      190, &
      1, &
      4, &
      68, &
      77, &
      85, &
      91, &
      190, &
      1, &
      4, &
      69, &
      79, &
      86, &
      190, &
      1, &
      4, &
      70, &
      80, &
      87, &
      95, &
      190, &
      1, &
      4, &
      71, &
      81, &
      88, &
      96, &
      190, &
      1, &
      4, &
      72, &
      82, &
      89, &
      97, &
      190, &
      1, &
      4, &
      73, &
      83, &
      90, &
      98, &
      190, &
      1, &
      4, &
      74, &
      84, &
      91, &
      99, &
      190, &
      1, &
      4, &
      75, &
      85, &
      92, &
      100, &
      190, &
      1, &
      4, &
      76, &
      93, &
      190, &
      1, &
      4, &
      77, &
      94, &
      102, &
      190, &
      1, &
      4, &
      86, &
      95, &
      190, &
      1, &
      4, &
      78, &
      87, &
      96, &
      190, &
      1, &
      4, &
      79, &
      88, &
      97, &
      106, &
      190, &
      1, &
      4, &
      80, &
      89, &
      98, &
      107, &
      190, &
      1, &
      4, &
      81, &
      90, &
      99, &
      108, &
      190, &
      1, &
      4, &
      82, &
      91, &
      100, &
      109, &
      190, &
      1, &
      4, &
      83, &
      92, &
      101, &
      110, &
      190, &
      1, &
      4, &
      84, &
      93, &
      102, &
      111, &
      190, &
      1, &
      4, &
      85, &
      94, &
      103, &
      112, &
      190, &
      1, &
      4, &
      104, &
      113, &
      190, &
      1, &
      4, &
      105, &
      114, &
      190, &
      1, &
      4, &
      96, &
      106, &
      190, &
      1, &
      4, &
      86, &
      97, &
      107, &
      190, &
      1, &
      4, &
      87, &
      98, &
      108, &
      116, &
      190, &
      1, &
      4, &
      88, &
      99, &
      109, &
      117, &
      190, &
      1, &
      4, &
      89, &
      100, &
      110, &
      118, &
      190, &
      1, &
      4, &
      90, &
      101, &
      111, &
      119, &
      190, &
      1, &
      4, &
      91, &
      102, &
      112, &
      120, &
      190, &
      1, &
      4, &
      92, &
      103, &
      113, &
      121, &
      190, &
      1, &
      4, &
      93, &
      104, &
      114, &
      122, &
      190, &
      1, &
      4, &
      94, &
      105, &
      115, &
      123, &
      190, &
      1, &
      4, &
      95, &
      107, &
      116, &
      190, &
      1, &
      4, &
      96, &
      108, &
      117, &
      190, &
      1, &
      4, &
      97, &
      109, &
      118, &
      127, &
      190, &
      1, &
      4, &
      98, &
      110, &
      119, &
      128, &
      190, &
      1, &
      4, &
      99, &
      111, &
      120, &
      129, &
      190, &
      1, &
      4, &
      100, &
      112, &
      121, &
      130, &
      190, &
      1, &
      4, &
      101, &
      113, &
      122, &
      131, &
      190, &
      1, &
      4, &
      102, &
      114, &
      123, &
      132, &
      190, &
      1, &
      4, &
      103, &
      115, &
      124, &
      133, &
      190, &
      1, &
      4, &
      104, &
      125, &
      134, &
      190, &
      1, &
      4, &
      105, &
      126, &
      135, &
      190, &
      1, &
      4, &
      106, &
      117, &
      127, &
      190, &
      1, &
      4, &
      107, &
      118, &
      128, &
      190, &
      1, &
      4, &
      108, &
      119, &
      129, &
      137, &
      190, &
      1, &
      4, &
      109, &
      120, &
      130, &
      138, &
      190, &
      1, &
      4, &
      110, &
      121, &
      131, &
      139, &
      190, &
      1, &
      4, &
      111, &
      122, &
      132, &
      140, &
      190, &
      1, &
      4, &
      112, &
      123, &
      133, &
      141, &
      190, &
      1, &
      4, &
      113, &
      124, &
      134, &
      142, &
      190, &
      1, &
      4, &
      114, &
      125, &
      135, &
      143, &
      190, &
      1, &
      4, &
      115, &
      126, &
      136, &
      144, &
      190, &
      1, &
      4, &
      116, &
      128, &
      137, &
      190, &
      1, &
      4, &
      117, &
      129, &
      138, &
      190, &
      1, &
      4, &
      118, &
      130, &
      139, &
      148, &
      190, &
      1, &
      4, &
      119, &
      131, &
      140, &
      149, &
      190, &
      1, &
      4, &
      120, &
      132, &
      141, &
      150, &
      190, &
      1, &
      4, &
      121, &
      133, &
      142, &
      151, &
      190, &
      1, &
      4, &
      122, &
      134, &
      143, &
      152, &
      190, &
      1, &
      4, &
      123, &
      135, &
      144, &
      153, &
      190, &
      1, &
      4, &
      124, &
      136, &
      145, &
      154, &
      190, &
      1, &
      4, &
      125, &
      146, &
      155, &
      190, &
      1, &
      4, &
      126, &
      147, &
      156, &
      190, &
      1, &
      4, &
      127, &
      138, &
      148, &
      190, &
      1, &
      4, &
      128, &
      139, &
      149, &
      159, &
      190, &
      1, &
      4, &
      129, &
      140, &
      150, &
      160, &
      190, &
      1, &
      4, &
      130, &
      141, &
      151, &
      161, &
      190, &
      1, &
      4, &
      131, &
      142, &
      152, &
      162, &
      190, &
      1, &
      4, &
      132, &
      143, &
      153, &
      163, &
      190, &
      1, &
      4, &
      133, &
      144, &
      154, &
      164, &
      190, &
      1, &
      4, &
      134, &
      145, &
      155, &
      165, &
      190, &
      1, &
      4, &
      135, &
      146, &
      156, &
      166, &
      190, &
      1, &
      4, &
      136, &
      147, &
      157, &
      167, &
      190, &
      1, &
      158, &
      169, &
      190, &
      1, &
      4, &
      148, &
      159, &
      190, &
      1, &
      4, &
      137, &
      149, &
      160, &
      190, &
      1, &
      4, &
      138, &
      150, &
      161, &
      171, &
      190, &
      1, &
      4, &
      139, &
      151, &
      162, &
      172, &
      190, &
      1, &
      4, &
      140, &
      152, &
      163, &
      173, &
      190, &
      1, &
      4, &
      141, &
      153, &
      164, &
      174, &
      190, &
      1, &
      4, &
      142, &
      154, &
      165, &
      175, &
      190, &
      1, &
      4, &
      143, &
      155, &
      166, &
      176, &
      190, &
      1, &
      4, &
      144, &
      156, &
      167, &
      177, &
      190, &
      1, &
      4, &
      145, &
      157, &
      168, &
      178, &
      190, &
      1, &
      4, &
      146, &
      169, &
      179, &
      190, &
      1, &
      4, &
      147, &
      158, &
      170, &
      180, &
      190, &
      1, &
      160, &
      171, &
      190, &
      1, &
      4, &
      148, &
      161, &
      172, &
      181, &
      190, &
      1, &
      4, &
      149, &
      162, &
      173, &
      182, &
      190, &
      1, &
      4, &
      150, &
      163, &
      174, &
      183, &
      190, &
      1, &
      4, &
      151, &
      164, &
      175, &
      184, &
      190, &
      1, &
      4, &
      152, &
      165, &
      176, &
      185, &
      190, &
      1, &
      4, &
      153, &
      166, &
      177, &
      186, &
      190, &
      1, &
      4, &
      154, &
      167, &
      178, &
      187, &
      190, &
      1, &
      4, &
      155, &
      168, &
      179, &
      188, &
      190, &
      1, &
      4, &
      156, &
      169, &
      180, &
      189, &
      190, &
      1, &
      171, &
      181, &
      190, &
      1, &
      4, &
      159, &
      172, &
      182, &
      190, &
      1, &
      4, &
      160, &
      173, &
      183, &
      190, &
      1, &
      4, &
      161, &
      174, &
      184, &
      190, &
      1, &
      4, &
      162, &
      175, &
      185, &
      190, &
      1, &
      4, &
      163, &
      176, &
      186, &
      190, &
      1, &
      4, &
      164, &
      177, &
      187, &
      190, &
      1, &
      4, &
      165, &
      178, &
      188, &
      190, &
      1, &
      4, &
      166, &
      179, &
      189, &
      190, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      189, &
      190, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      10, &
      11, &
      12, &
      13, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      21, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      29, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      36, &
      37, &
      38, &
      39, &
      40, &
      41, &
      42, &
      43, &
      44, &
      45, &
      46, &
      47, &
      48, &
      49, &
      50, &
      51, &
      52, &
      53, &
      54, &
      55, &
      56, &
      57, &
      58, &
      59, &
      60, &
      61, &
      62, &
      63, &
      64, &
      65, &
      66, &
      67, &
      68, &
      69, &
      70, &
      71, &
      72, &
      73, &
      74, &
      75, &
      76, &
      77, &
      78, &
      79, &
      80, &
      81, &
      82, &
      83, &
      84, &
      85, &
      86, &
      87, &
      88, &
      89, &
      90, &
      91, &
      92, &
      93, &
      94, &
      95, &
      96, &
      97, &
      98, &
      99, &
      100, &
      101, &
      102, &
      103, &
      104, &
      105, &
      106, &
      107, &
      108, &
      109, &
      110, &
      111, &
      112, &
      113, &
      114, &
      115, &
      116, &
      117, &
      118, &
      119, &
      120, &
      121, &
      122, &
      123, &
      124, &
      125, &
      126, &
      127, &
      128, &
      129, &
      130, &
      131, &
      132, &
      133, &
      134, &
      135, &
      136, &
      137, &
      138, &
      139, &
      140, &
      141, &
      142, &
      143, &
      144, &
      145, &
      146, &
      147, &
      148, &
      149, &
      150, &
      151, &
      152, &
      153, &
      154, &
      155, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      162, &
      163, &
      164, &
      165, &
      166, &
      167, &
      168, &
      169, &
      170, &
      171, &
      172, &
      173, &
      174, &
      175, &
      176, &
      177, &
      178, &
      179, &
      180, &
      181, &
      182, &
      183, &
      184, &
      185, &
      186, &
      187, &
      188, &
      189, &
      190, &
      191  ]

    csr_jac_row_count = [ &
      1, &
      166, &
      172, &
      178, &
      341, &
      346, &
      353, &
      356, &
      360, &
      368, &
      372, &
      379, &
      388, &
      397, &
      404, &
      411, &
      421, &
      430, &
      437, &
      442, &
      450, &
      460, &
      469, &
      475, &
      481, &
      488, &
      498, &
      506, &
      512, &
      517, &
      523, &
      530, &
      537, &
      544, &
      549, &
      554, &
      561, &
      568, &
      575, &
      582, &
      588, &
      593, &
      599, &
      606, &
      613, &
      620, &
      627, &
      632, &
      637, &
      643, &
      650, &
      657, &
      664, &
      671, &
      677, &
      683, &
      689, &
      696, &
      703, &
      710, &
      717, &
      722, &
      728, &
      735, &
      742, &
      749, &
      756, &
      763, &
      769, &
      774, &
      780, &
      787, &
      794, &
      801, &
      808, &
      815, &
      820, &
      825, &
      830, &
      836, &
      843, &
      850, &
      857, &
      864, &
      871, &
      878, &
      884, &
      891, &
      898, &
      905, &
      912, &
      919, &
      926, &
      931, &
      937, &
      942, &
      948, &
      955, &
      962, &
      969, &
      976, &
      983, &
      990, &
      997, &
      1002, &
      1007, &
      1012, &
      1018, &
      1025, &
      1032, &
      1039, &
      1046, &
      1053, &
      1060, &
      1067, &
      1074, &
      1080, &
      1086, &
      1093, &
      1100, &
      1107, &
      1114, &
      1121, &
      1128, &
      1135, &
      1141, &
      1147, &
      1153, &
      1159, &
      1166, &
      1173, &
      1180, &
      1187, &
      1194, &
      1201, &
      1208, &
      1215, &
      1221, &
      1227, &
      1234, &
      1241, &
      1248, &
      1255, &
      1262, &
      1269, &
      1276, &
      1282, &
      1288, &
      1294, &
      1301, &
      1308, &
      1315, &
      1322, &
      1329, &
      1336, &
      1343, &
      1350, &
      1357, &
      1361, &
      1366, &
      1372, &
      1379, &
      1386, &
      1393, &
      1400, &
      1407, &
      1414, &
      1421, &
      1428, &
      1434, &
      1441, &
      1445, &
      1452, &
      1459, &
      1466, &
      1473, &
      1480, &
      1487, &
      1494, &
      1501, &
      1508, &
      1512, &
      1518, &
      1524, &
      1530, &
      1536, &
      1542, &
      1548, &
      1554, &
      1560, &
      1750, &
      1941  ]
#endif

  end subroutine actual_network_init


  subroutine actual_network_finalize()
    ! Deallocate storage arrays

    if (allocated(aion)) then
       deallocate(aion)
    endif

    if (allocated(zion)) then
       deallocate(zion)
    endif

    if (allocated(bion)) then
       deallocate(bion)
    endif

    if (allocated(nion)) then
       deallocate(nion)
    endif

    if (allocated(mion)) then
       deallocate(mion)
    endif

    if (allocated(wion)) then
       deallocate(wion)
    endif

#ifdef REACT_SPARSE_JACOBIAN
    if (allocated(csr_jac_col_index)) then
       deallocate(csr_jac_col_index)
    endif

    if (allocated(csr_jac_row_count)) then
       deallocate(csr_jac_row_count)
    endif
#endif

  end subroutine actual_network_finalize


  subroutine ener_gener_rate(dydt, enuc)
    ! Computes the instantaneous energy generation rate

    !$acc routine seq

    implicit none

    real(rt) :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * mion(:)) * enuc_conv2

  end subroutine ener_gener_rate

end module actual_network
