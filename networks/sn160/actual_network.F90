module actual_network

  use physical_constants, only: ERG_PER_MeV
  use microphysics_type_module

  implicit none

  public

  real(rt), parameter :: avo = 6.0221417930e23_rt
  real(rt), parameter :: c_light = 2.99792458e10_rt
  real(rt), parameter :: enuc_conv2 = -avo*c_light*c_light

  real(rt), parameter :: ev2erg  = 1.60217648740e-12_rt
  real(rt), parameter :: mev2erg = ev2erg*1.0e6_rt
  real(rt), parameter :: mev2gr  = mev2erg/c_light**2

  real(rt), parameter :: mass_neutron  = 1.67492721184e-24_rt
  real(rt), parameter :: mass_proton   = 1.67262163783e-24_rt
  real(rt), parameter :: mass_electron = 9.10938215450e-28_rt

  integer, parameter :: nrates = 1538
  integer, parameter :: num_rate_groups = 4

  ! Evolution and auxiliary
  integer, parameter :: nspec_evolve = 160
  integer, parameter :: naux  = 0

  ! Number of nuclear species in the network
  integer, parameter :: nspec = 160

  ! Number of reaclib rates
  integer, parameter :: nrat_reaclib = 1538
  integer, parameter :: number_reaclib_sets = 1956

  ! Number of tabular rates
  integer, parameter :: nrat_tabular = 0

  ! Binding Energies Per Nucleon (MeV)
  real(rt) :: ebind_per_nucleon(nspec)

  ! aion: Nucleon mass number A
  ! zion: Nucleon atomic number Z
  ! nion: Nucleon neutron number N
  ! bion: Binding Energies (ergs)

  ! Nuclides
  integer, parameter :: jn   = 1
  integer, parameter :: jp   = 2
  integer, parameter :: jd   = 3
  integer, parameter :: jhe3   = 4
  integer, parameter :: jhe4   = 5
  integer, parameter :: jli6   = 6
  integer, parameter :: jli7   = 7
  integer, parameter :: jbe7   = 8
  integer, parameter :: jbe9   = 9
  integer, parameter :: jb8   = 10
  integer, parameter :: jb10   = 11
  integer, parameter :: jb11   = 12
  integer, parameter :: jc12   = 13
  integer, parameter :: jc13   = 14
  integer, parameter :: jc14   = 15
  integer, parameter :: jn13   = 16
  integer, parameter :: jn14   = 17
  integer, parameter :: jn15   = 18
  integer, parameter :: jo14   = 19
  integer, parameter :: jo15   = 20
  integer, parameter :: jo16   = 21
  integer, parameter :: jo17   = 22
  integer, parameter :: jo18   = 23
  integer, parameter :: jf17   = 24
  integer, parameter :: jf18   = 25
  integer, parameter :: jf19   = 26
  integer, parameter :: jne18   = 27
  integer, parameter :: jne19   = 28
  integer, parameter :: jne20   = 29
  integer, parameter :: jne21   = 30
  integer, parameter :: jne22   = 31
  integer, parameter :: jna21   = 32
  integer, parameter :: jna22   = 33
  integer, parameter :: jna23   = 34
  integer, parameter :: jmg23   = 35
  integer, parameter :: jmg24   = 36
  integer, parameter :: jmg25   = 37
  integer, parameter :: jmg26   = 38
  integer, parameter :: jal25   = 39
  integer, parameter :: jal26   = 40
  integer, parameter :: jal27   = 41
  integer, parameter :: jsi28   = 42
  integer, parameter :: jsi29   = 43
  integer, parameter :: jsi30   = 44
  integer, parameter :: jsi31   = 45
  integer, parameter :: jsi32   = 46
  integer, parameter :: jp29   = 47
  integer, parameter :: jp30   = 48
  integer, parameter :: jp31   = 49
  integer, parameter :: jp32   = 50
  integer, parameter :: jp33   = 51
  integer, parameter :: js32   = 52
  integer, parameter :: js33   = 53
  integer, parameter :: js34   = 54
  integer, parameter :: js35   = 55
  integer, parameter :: js36   = 56
  integer, parameter :: jcl33   = 57
  integer, parameter :: jcl34   = 58
  integer, parameter :: jcl35   = 59
  integer, parameter :: jcl36   = 60
  integer, parameter :: jcl37   = 61
  integer, parameter :: jar36   = 62
  integer, parameter :: jar37   = 63
  integer, parameter :: jar38   = 64
  integer, parameter :: jar39   = 65
  integer, parameter :: jar40   = 66
  integer, parameter :: jk37   = 67
  integer, parameter :: jk38   = 68
  integer, parameter :: jk39   = 69
  integer, parameter :: jk40   = 70
  integer, parameter :: jk41   = 71
  integer, parameter :: jca40   = 72
  integer, parameter :: jca41   = 73
  integer, parameter :: jca42   = 74
  integer, parameter :: jca43   = 75
  integer, parameter :: jca44   = 76
  integer, parameter :: jca45   = 77
  integer, parameter :: jca46   = 78
  integer, parameter :: jca47   = 79
  integer, parameter :: jca48   = 80
  integer, parameter :: jsc43   = 81
  integer, parameter :: jsc44   = 82
  integer, parameter :: jsc45   = 83
  integer, parameter :: jsc46   = 84
  integer, parameter :: jsc47   = 85
  integer, parameter :: jsc48   = 86
  integer, parameter :: jsc49   = 87
  integer, parameter :: jti44   = 88
  integer, parameter :: jti45   = 89
  integer, parameter :: jti46   = 90
  integer, parameter :: jti47   = 91
  integer, parameter :: jti48   = 92
  integer, parameter :: jti49   = 93
  integer, parameter :: jti50   = 94
  integer, parameter :: jti51   = 95
  integer, parameter :: jv46   = 96
  integer, parameter :: jv47   = 97
  integer, parameter :: jv48   = 98
  integer, parameter :: jv49   = 99
  integer, parameter :: jv50   = 100
  integer, parameter :: jv51   = 101
  integer, parameter :: jv52   = 102
  integer, parameter :: jcr48   = 103
  integer, parameter :: jcr49   = 104
  integer, parameter :: jcr50   = 105
  integer, parameter :: jcr51   = 106
  integer, parameter :: jcr52   = 107
  integer, parameter :: jcr53   = 108
  integer, parameter :: jcr54   = 109
  integer, parameter :: jmn50   = 110
  integer, parameter :: jmn51   = 111
  integer, parameter :: jmn52   = 112
  integer, parameter :: jmn53   = 113
  integer, parameter :: jmn54   = 114
  integer, parameter :: jmn55   = 115
  integer, parameter :: jfe52   = 116
  integer, parameter :: jfe53   = 117
  integer, parameter :: jfe54   = 118
  integer, parameter :: jfe55   = 119
  integer, parameter :: jfe56   = 120
  integer, parameter :: jfe57   = 121
  integer, parameter :: jfe58   = 122
  integer, parameter :: jco53   = 123
  integer, parameter :: jco54   = 124
  integer, parameter :: jco55   = 125
  integer, parameter :: jco56   = 126
  integer, parameter :: jco57   = 127
  integer, parameter :: jco58   = 128
  integer, parameter :: jco59   = 129
  integer, parameter :: jni56   = 130
  integer, parameter :: jni57   = 131
  integer, parameter :: jni58   = 132
  integer, parameter :: jni59   = 133
  integer, parameter :: jni60   = 134
  integer, parameter :: jni61   = 135
  integer, parameter :: jni62   = 136
  integer, parameter :: jni63   = 137
  integer, parameter :: jni64   = 138
  integer, parameter :: jcu57   = 139
  integer, parameter :: jcu58   = 140
  integer, parameter :: jcu59   = 141
  integer, parameter :: jcu60   = 142
  integer, parameter :: jcu61   = 143
  integer, parameter :: jcu62   = 144
  integer, parameter :: jcu63   = 145
  integer, parameter :: jcu64   = 146
  integer, parameter :: jcu65   = 147
  integer, parameter :: jzn59   = 148
  integer, parameter :: jzn60   = 149
  integer, parameter :: jzn61   = 150
  integer, parameter :: jzn62   = 151
  integer, parameter :: jzn63   = 152
  integer, parameter :: jzn64   = 153
  integer, parameter :: jzn65   = 154
  integer, parameter :: jzn66   = 155
  integer, parameter :: jga62   = 156
  integer, parameter :: jga63   = 157
  integer, parameter :: jga64   = 158
  integer, parameter :: jge63   = 159
  integer, parameter :: jge64   = 160

  ! Reactions
  integer, parameter :: k_n__p__weak__wc12   = 1
  integer, parameter :: k_be7__li7__weak__electron_capture   = 2
  integer, parameter :: k_c14__n14__weak__wc12   = 3
  integer, parameter :: k_n13__c13__weak__wc12   = 4
  integer, parameter :: k_o14__n14__weak__wc12   = 5
  integer, parameter :: k_o15__n15__weak__wc12   = 6
  integer, parameter :: k_f17__o17__weak__wc12   = 7
  integer, parameter :: k_f18__o18__weak__wc12   = 8
  integer, parameter :: k_ne18__f18__weak__wc12   = 9
  integer, parameter :: k_ne19__f19__weak__wc12   = 10
  integer, parameter :: k_na21__ne21__weak__wc12   = 11
  integer, parameter :: k_na22__ne22__weak__wc12   = 12
  integer, parameter :: k_mg23__na23__weak__wc12   = 13
  integer, parameter :: k_al25__mg25__weak__wc12   = 14
  integer, parameter :: k_al26__mg26__weak__wc12   = 15
  integer, parameter :: k_si31__p31__weak__wc12   = 16
  integer, parameter :: k_si32__p32__weak__wc12   = 17
  integer, parameter :: k_p29__si29__weak__wc12   = 18
  integer, parameter :: k_p30__si30__weak__wc12   = 19
  integer, parameter :: k_p32__s32__weak__wc12   = 20
  integer, parameter :: k_p33__s33__weak__wc12   = 21
  integer, parameter :: k_s35__cl35__weak__wc12   = 22
  integer, parameter :: k_cl33__s33__weak__wc12   = 23
  integer, parameter :: k_cl34__s34__weak__wc12   = 24
  integer, parameter :: k_cl36__ar36__weak__wc12   = 25
  integer, parameter :: k_cl36__s36__weak__wc12   = 26
  integer, parameter :: k_ar37__cl37__weak__wc12   = 27
  integer, parameter :: k_ar39__k39__weak__wc12   = 28
  integer, parameter :: k_k37__ar37__weak__wc12   = 29
  integer, parameter :: k_k38__ar38__weak__wc12   = 30
  integer, parameter :: k_k40__ca40__weak__wc12   = 31
  integer, parameter :: k_k40__ar40__weak__wc12   = 32
  integer, parameter :: k_ca41__k41__weak__wc12   = 33
  integer, parameter :: k_ca45__sc45__weak__wc12   = 34
  integer, parameter :: k_ca47__sc47__weak__wc12   = 35
  integer, parameter :: k_ca48__sc48__weak__mo03   = 36
  integer, parameter :: k_sc43__ca43__weak__wc12   = 37
  integer, parameter :: k_sc44__ca44__weak__wc12   = 38
  integer, parameter :: k_sc46__ti46__weak__wc12   = 39
  integer, parameter :: k_sc47__ti47__weak__wc12   = 40
  integer, parameter :: k_sc48__ti48__weak__wc12   = 41
  integer, parameter :: k_sc49__ti49__weak__wc12   = 42
  integer, parameter :: k_ti44__sc44__weak__wc12   = 43
  integer, parameter :: k_ti45__sc45__weak__wc12   = 44
  integer, parameter :: k_ti51__v51__weak__wc12   = 45
  integer, parameter :: k_v46__ti46__weak__wc12   = 46
  integer, parameter :: k_v47__ti47__weak__wc12   = 47
  integer, parameter :: k_v48__ti48__weak__wc12   = 48
  integer, parameter :: k_v49__ti49__weak__wc12   = 49
  integer, parameter :: k_v50__ti50__weak__mo03   = 50
  integer, parameter :: k_v52__cr52__weak__wc12   = 51
  integer, parameter :: k_cr48__v48__weak__wc12   = 52
  integer, parameter :: k_cr49__v49__weak__wc12   = 53
  integer, parameter :: k_cr51__v51__weak__wc12   = 54
  integer, parameter :: k_mn50__cr50__weak__wc12   = 55
  integer, parameter :: k_mn51__cr51__weak__wc12   = 56
  integer, parameter :: k_mn52__cr52__weak__wc12   = 57
  integer, parameter :: k_mn53__cr53__weak__wc12   = 58
  integer, parameter :: k_mn54__cr54__weak__wc12   = 59
  integer, parameter :: k_fe52__mn52__weak__wc12   = 60
  integer, parameter :: k_fe53__mn53__weak__wc12   = 61
  integer, parameter :: k_fe55__mn55__weak__wc12   = 62
  integer, parameter :: k_co53__fe53__weak__wc12   = 63
  integer, parameter :: k_co54__fe54__weak__wc12   = 64
  integer, parameter :: k_co55__fe55__weak__wc12   = 65
  integer, parameter :: k_co56__fe56__weak__wc12   = 66
  integer, parameter :: k_co57__fe57__weak__wc12   = 67
  integer, parameter :: k_co58__fe58__weak__wc12   = 68
  integer, parameter :: k_ni56__co56__weak__wc12   = 69
  integer, parameter :: k_ni57__co57__weak__wc12   = 70
  integer, parameter :: k_ni59__co59__weak__wc12   = 71
  integer, parameter :: k_ni63__cu63__weak__wc12   = 72
  integer, parameter :: k_cu57__ni57__weak__wc12   = 73
  integer, parameter :: k_cu58__ni58__weak__wc12   = 74
  integer, parameter :: k_cu59__ni59__weak__wc12   = 75
  integer, parameter :: k_cu60__ni60__weak__wc12   = 76
  integer, parameter :: k_cu61__ni61__weak__wc12   = 77
  integer, parameter :: k_cu62__ni62__weak__wc12   = 78
  integer, parameter :: k_cu64__ni64__weak__wc12   = 79
  integer, parameter :: k_cu64__zn64__weak__wc12   = 80
  integer, parameter :: k_zn59__cu59__weak__wc12   = 81
  integer, parameter :: k_zn60__cu60__weak__wc12   = 82
  integer, parameter :: k_zn61__cu61__weak__wc12   = 83
  integer, parameter :: k_zn62__cu62__weak__wc12   = 84
  integer, parameter :: k_zn63__cu63__weak__wc12   = 85
  integer, parameter :: k_zn65__cu65__weak__wc12   = 86
  integer, parameter :: k_ga62__zn62__weak__wc12   = 87
  integer, parameter :: k_ga63__zn63__weak__wc12   = 88
  integer, parameter :: k_ga64__zn64__weak__wc12   = 89
  integer, parameter :: k_ge63__ga63__weak__wc12   = 90
  integer, parameter :: k_ge64__ga64__weak__wc12   = 91
  integer, parameter :: k_d__n_p   = 92
  integer, parameter :: k_he3__p_d   = 93
  integer, parameter :: k_he4__n_he3   = 94
  integer, parameter :: k_he4__d_d   = 95
  integer, parameter :: k_li6__he4_d   = 96
  integer, parameter :: k_li7__n_li6   = 97
  integer, parameter :: k_be7__p_li6   = 98
  integer, parameter :: k_be7__he4_he3   = 99
  integer, parameter :: k_b8__p_be7   = 100
  integer, parameter :: k_b8__he4_he4__weak__wc12   = 101
  integer, parameter :: k_b10__p_be9   = 102
  integer, parameter :: k_b10__he4_li6   = 103
  integer, parameter :: k_b11__n_b10   = 104
  integer, parameter :: k_b11__he4_li7   = 105
  integer, parameter :: k_c12__p_b11   = 106
  integer, parameter :: k_c13__n_c12   = 107
  integer, parameter :: k_c14__n_c13   = 108
  integer, parameter :: k_n13__p_c12   = 109
  integer, parameter :: k_n14__n_n13   = 110
  integer, parameter :: k_n14__p_c13   = 111
  integer, parameter :: k_n15__n_n14   = 112
  integer, parameter :: k_n15__p_c14   = 113
  integer, parameter :: k_o14__p_n13   = 114
  integer, parameter :: k_o15__n_o14   = 115
  integer, parameter :: k_o15__p_n14   = 116
  integer, parameter :: k_o16__n_o15   = 117
  integer, parameter :: k_o16__p_n15   = 118
  integer, parameter :: k_o16__he4_c12   = 119
  integer, parameter :: k_o17__n_o16   = 120
  integer, parameter :: k_o18__n_o17   = 121
  integer, parameter :: k_o18__he4_c14   = 122
  integer, parameter :: k_f17__p_o16   = 123
  integer, parameter :: k_f18__n_f17   = 124
  integer, parameter :: k_f18__p_o17   = 125
  integer, parameter :: k_f18__he4_n14   = 126
  integer, parameter :: k_f19__n_f18   = 127
  integer, parameter :: k_f19__p_o18   = 128
  integer, parameter :: k_f19__he4_n15   = 129
  integer, parameter :: k_ne18__p_f17   = 130
  integer, parameter :: k_ne18__he4_o14   = 131
  integer, parameter :: k_ne19__n_ne18   = 132
  integer, parameter :: k_ne19__p_f18   = 133
  integer, parameter :: k_ne19__he4_o15   = 134
  integer, parameter :: k_ne20__n_ne19   = 135
  integer, parameter :: k_ne20__p_f19   = 136
  integer, parameter :: k_ne20__he4_o16   = 137
  integer, parameter :: k_ne21__n_ne20   = 138
  integer, parameter :: k_ne21__he4_o17   = 139
  integer, parameter :: k_ne22__n_ne21   = 140
  integer, parameter :: k_ne22__he4_o18   = 141
  integer, parameter :: k_na21__p_ne20   = 142
  integer, parameter :: k_na21__he4_f17   = 143
  integer, parameter :: k_na22__n_na21   = 144
  integer, parameter :: k_na22__p_ne21   = 145
  integer, parameter :: k_na22__he4_f18   = 146
  integer, parameter :: k_na23__n_na22   = 147
  integer, parameter :: k_na23__p_ne22   = 148
  integer, parameter :: k_na23__he4_f19   = 149
  integer, parameter :: k_mg23__p_na22   = 150
  integer, parameter :: k_mg23__he4_ne19   = 151
  integer, parameter :: k_mg24__n_mg23   = 152
  integer, parameter :: k_mg24__p_na23   = 153
  integer, parameter :: k_mg24__he4_ne20   = 154
  integer, parameter :: k_mg25__n_mg24   = 155
  integer, parameter :: k_mg25__he4_ne21   = 156
  integer, parameter :: k_mg26__n_mg25   = 157
  integer, parameter :: k_mg26__he4_ne22   = 158
  integer, parameter :: k_al25__p_mg24   = 159
  integer, parameter :: k_al25__he4_na21   = 160
  integer, parameter :: k_al26__n_al25   = 161
  integer, parameter :: k_al26__p_mg25   = 162
  integer, parameter :: k_al26__he4_na22   = 163
  integer, parameter :: k_al27__n_al26   = 164
  integer, parameter :: k_al27__p_mg26   = 165
  integer, parameter :: k_al27__he4_na23   = 166
  integer, parameter :: k_si28__p_al27   = 167
  integer, parameter :: k_si28__he4_mg24   = 168
  integer, parameter :: k_si29__n_si28   = 169
  integer, parameter :: k_si29__he4_mg25   = 170
  integer, parameter :: k_si30__n_si29   = 171
  integer, parameter :: k_si30__he4_mg26   = 172
  integer, parameter :: k_si31__n_si30   = 173
  integer, parameter :: k_si32__n_si31   = 174
  integer, parameter :: k_p29__p_si28   = 175
  integer, parameter :: k_p29__he4_al25   = 176
  integer, parameter :: k_p30__n_p29   = 177
  integer, parameter :: k_p30__p_si29   = 178
  integer, parameter :: k_p30__he4_al26   = 179
  integer, parameter :: k_p31__n_p30   = 180
  integer, parameter :: k_p31__p_si30   = 181
  integer, parameter :: k_p31__he4_al27   = 182
  integer, parameter :: k_p32__n_p31   = 183
  integer, parameter :: k_p32__p_si31   = 184
  integer, parameter :: k_p33__n_p32   = 185
  integer, parameter :: k_p33__p_si32   = 186
  integer, parameter :: k_s32__p_p31   = 187
  integer, parameter :: k_s32__he4_si28   = 188
  integer, parameter :: k_s33__n_s32   = 189
  integer, parameter :: k_s33__p_p32   = 190
  integer, parameter :: k_s33__he4_si29   = 191
  integer, parameter :: k_s34__n_s33   = 192
  integer, parameter :: k_s34__p_p33   = 193
  integer, parameter :: k_s34__he4_si30   = 194
  integer, parameter :: k_s35__n_s34   = 195
  integer, parameter :: k_s35__he4_si31   = 196
  integer, parameter :: k_s36__n_s35   = 197
  integer, parameter :: k_s36__he4_si32   = 198
  integer, parameter :: k_cl33__p_s32   = 199
  integer, parameter :: k_cl33__he4_p29   = 200
  integer, parameter :: k_cl34__n_cl33   = 201
  integer, parameter :: k_cl34__p_s33   = 202
  integer, parameter :: k_cl34__he4_p30   = 203
  integer, parameter :: k_cl35__n_cl34   = 204
  integer, parameter :: k_cl35__p_s34   = 205
  integer, parameter :: k_cl35__he4_p31   = 206
  integer, parameter :: k_cl36__n_cl35   = 207
  integer, parameter :: k_cl36__p_s35   = 208
  integer, parameter :: k_cl36__he4_p32   = 209
  integer, parameter :: k_cl37__n_cl36   = 210
  integer, parameter :: k_cl37__p_s36   = 211
  integer, parameter :: k_cl37__he4_p33   = 212
  integer, parameter :: k_ar36__p_cl35   = 213
  integer, parameter :: k_ar36__he4_s32   = 214
  integer, parameter :: k_ar37__n_ar36   = 215
  integer, parameter :: k_ar37__p_cl36   = 216
  integer, parameter :: k_ar37__he4_s33   = 217
  integer, parameter :: k_ar38__n_ar37   = 218
  integer, parameter :: k_ar38__p_cl37   = 219
  integer, parameter :: k_ar38__he4_s34   = 220
  integer, parameter :: k_ar39__n_ar38   = 221
  integer, parameter :: k_ar39__he4_s35   = 222
  integer, parameter :: k_ar40__n_ar39   = 223
  integer, parameter :: k_ar40__he4_s36   = 224
  integer, parameter :: k_k37__p_ar36   = 225
  integer, parameter :: k_k37__he4_cl33   = 226
  integer, parameter :: k_k38__n_k37   = 227
  integer, parameter :: k_k38__p_ar37   = 228
  integer, parameter :: k_k38__he4_cl34   = 229
  integer, parameter :: k_k39__n_k38   = 230
  integer, parameter :: k_k39__p_ar38   = 231
  integer, parameter :: k_k39__he4_cl35   = 232
  integer, parameter :: k_k40__n_k39   = 233
  integer, parameter :: k_k40__p_ar39   = 234
  integer, parameter :: k_k40__he4_cl36   = 235
  integer, parameter :: k_k41__n_k40   = 236
  integer, parameter :: k_k41__p_ar40   = 237
  integer, parameter :: k_k41__he4_cl37   = 238
  integer, parameter :: k_ca40__p_k39   = 239
  integer, parameter :: k_ca40__he4_ar36   = 240
  integer, parameter :: k_ca41__n_ca40   = 241
  integer, parameter :: k_ca41__p_k40   = 242
  integer, parameter :: k_ca41__he4_ar37   = 243
  integer, parameter :: k_ca42__n_ca41   = 244
  integer, parameter :: k_ca42__p_k41   = 245
  integer, parameter :: k_ca42__he4_ar38   = 246
  integer, parameter :: k_ca43__n_ca42   = 247
  integer, parameter :: k_ca43__he4_ar39   = 248
  integer, parameter :: k_ca44__n_ca43   = 249
  integer, parameter :: k_ca44__he4_ar40   = 250
  integer, parameter :: k_ca45__n_ca44   = 251
  integer, parameter :: k_ca46__n_ca45   = 252
  integer, parameter :: k_ca47__n_ca46   = 253
  integer, parameter :: k_ca48__n_ca47   = 254
  integer, parameter :: k_sc43__p_ca42   = 255
  integer, parameter :: k_sc43__he4_k39   = 256
  integer, parameter :: k_sc44__n_sc43   = 257
  integer, parameter :: k_sc44__p_ca43   = 258
  integer, parameter :: k_sc44__he4_k40   = 259
  integer, parameter :: k_sc45__n_sc44   = 260
  integer, parameter :: k_sc45__p_ca44   = 261
  integer, parameter :: k_sc45__he4_k41   = 262
  integer, parameter :: k_sc46__n_sc45   = 263
  integer, parameter :: k_sc46__p_ca45   = 264
  integer, parameter :: k_sc47__n_sc46   = 265
  integer, parameter :: k_sc47__p_ca46   = 266
  integer, parameter :: k_sc48__n_sc47   = 267
  integer, parameter :: k_sc48__p_ca47   = 268
  integer, parameter :: k_sc49__n_sc48   = 269
  integer, parameter :: k_sc49__p_ca48   = 270
  integer, parameter :: k_ti44__p_sc43   = 271
  integer, parameter :: k_ti44__he4_ca40   = 272
  integer, parameter :: k_ti45__n_ti44   = 273
  integer, parameter :: k_ti45__p_sc44   = 274
  integer, parameter :: k_ti45__he4_ca41   = 275
  integer, parameter :: k_ti46__n_ti45   = 276
  integer, parameter :: k_ti46__p_sc45   = 277
  integer, parameter :: k_ti46__he4_ca42   = 278
  integer, parameter :: k_ti47__n_ti46   = 279
  integer, parameter :: k_ti47__p_sc46   = 280
  integer, parameter :: k_ti47__he4_ca43   = 281
  integer, parameter :: k_ti48__n_ti47   = 282
  integer, parameter :: k_ti48__p_sc47   = 283
  integer, parameter :: k_ti48__he4_ca44   = 284
  integer, parameter :: k_ti49__n_ti48   = 285
  integer, parameter :: k_ti49__p_sc48   = 286
  integer, parameter :: k_ti49__he4_ca45   = 287
  integer, parameter :: k_ti50__n_ti49   = 288
  integer, parameter :: k_ti50__p_sc49   = 289
  integer, parameter :: k_ti50__he4_ca46   = 290
  integer, parameter :: k_ti51__n_ti50   = 291
  integer, parameter :: k_ti51__he4_ca47   = 292
  integer, parameter :: k_v46__p_ti45   = 293
  integer, parameter :: k_v47__n_v46   = 294
  integer, parameter :: k_v47__p_ti46   = 295
  integer, parameter :: k_v47__he4_sc43   = 296
  integer, parameter :: k_v48__n_v47   = 297
  integer, parameter :: k_v48__p_ti47   = 298
  integer, parameter :: k_v48__he4_sc44   = 299
  integer, parameter :: k_v49__n_v48   = 300
  integer, parameter :: k_v49__p_ti48   = 301
  integer, parameter :: k_v49__he4_sc45   = 302
  integer, parameter :: k_v50__n_v49   = 303
  integer, parameter :: k_v50__p_ti49   = 304
  integer, parameter :: k_v50__he4_sc46   = 305
  integer, parameter :: k_v51__n_v50   = 306
  integer, parameter :: k_v51__p_ti50   = 307
  integer, parameter :: k_v51__he4_sc47   = 308
  integer, parameter :: k_v52__n_v51   = 309
  integer, parameter :: k_v52__p_ti51   = 310
  integer, parameter :: k_v52__he4_sc48   = 311
  integer, parameter :: k_cr48__p_v47   = 312
  integer, parameter :: k_cr48__he4_ti44   = 313
  integer, parameter :: k_cr49__n_cr48   = 314
  integer, parameter :: k_cr49__p_v48   = 315
  integer, parameter :: k_cr49__he4_ti45   = 316
  integer, parameter :: k_cr50__n_cr49   = 317
  integer, parameter :: k_cr50__p_v49   = 318
  integer, parameter :: k_cr50__he4_ti46   = 319
  integer, parameter :: k_cr51__n_cr50   = 320
  integer, parameter :: k_cr51__p_v50   = 321
  integer, parameter :: k_cr51__he4_ti47   = 322
  integer, parameter :: k_cr52__n_cr51   = 323
  integer, parameter :: k_cr52__p_v51   = 324
  integer, parameter :: k_cr52__he4_ti48   = 325
  integer, parameter :: k_cr53__n_cr52   = 326
  integer, parameter :: k_cr53__p_v52   = 327
  integer, parameter :: k_cr53__he4_ti49   = 328
  integer, parameter :: k_cr54__n_cr53   = 329
  integer, parameter :: k_cr54__he4_ti50   = 330
  integer, parameter :: k_mn50__p_cr49   = 331
  integer, parameter :: k_mn50__he4_v46   = 332
  integer, parameter :: k_mn51__n_mn50   = 333
  integer, parameter :: k_mn51__p_cr50   = 334
  integer, parameter :: k_mn51__he4_v47   = 335
  integer, parameter :: k_mn52__n_mn51   = 336
  integer, parameter :: k_mn52__p_cr51   = 337
  integer, parameter :: k_mn52__he4_v48   = 338
  integer, parameter :: k_mn53__n_mn52   = 339
  integer, parameter :: k_mn53__p_cr52   = 340
  integer, parameter :: k_mn53__he4_v49   = 341
  integer, parameter :: k_mn54__n_mn53   = 342
  integer, parameter :: k_mn54__p_cr53   = 343
  integer, parameter :: k_mn54__he4_v50   = 344
  integer, parameter :: k_mn55__n_mn54   = 345
  integer, parameter :: k_mn55__p_cr54   = 346
  integer, parameter :: k_mn55__he4_v51   = 347
  integer, parameter :: k_fe52__p_mn51   = 348
  integer, parameter :: k_fe52__he4_cr48   = 349
  integer, parameter :: k_fe53__n_fe52   = 350
  integer, parameter :: k_fe53__p_mn52   = 351
  integer, parameter :: k_fe53__he4_cr49   = 352
  integer, parameter :: k_fe54__n_fe53   = 353
  integer, parameter :: k_fe54__p_mn53   = 354
  integer, parameter :: k_fe54__he4_cr50   = 355
  integer, parameter :: k_fe55__n_fe54   = 356
  integer, parameter :: k_fe55__p_mn54   = 357
  integer, parameter :: k_fe55__he4_cr51   = 358
  integer, parameter :: k_fe56__n_fe55   = 359
  integer, parameter :: k_fe56__p_mn55   = 360
  integer, parameter :: k_fe56__he4_cr52   = 361
  integer, parameter :: k_fe57__n_fe56   = 362
  integer, parameter :: k_fe57__he4_cr53   = 363
  integer, parameter :: k_fe58__n_fe57   = 364
  integer, parameter :: k_fe58__he4_cr54   = 365
  integer, parameter :: k_co53__p_fe52   = 366
  integer, parameter :: k_co54__n_co53   = 367
  integer, parameter :: k_co54__p_fe53   = 368
  integer, parameter :: k_co54__he4_mn50   = 369
  integer, parameter :: k_co55__n_co54   = 370
  integer, parameter :: k_co55__p_fe54   = 371
  integer, parameter :: k_co55__he4_mn51   = 372
  integer, parameter :: k_co56__n_co55   = 373
  integer, parameter :: k_co56__p_fe55   = 374
  integer, parameter :: k_co56__he4_mn52   = 375
  integer, parameter :: k_co57__n_co56   = 376
  integer, parameter :: k_co57__p_fe56   = 377
  integer, parameter :: k_co57__he4_mn53   = 378
  integer, parameter :: k_co58__n_co57   = 379
  integer, parameter :: k_co58__p_fe57   = 380
  integer, parameter :: k_co58__he4_mn54   = 381
  integer, parameter :: k_co59__n_co58   = 382
  integer, parameter :: k_co59__p_fe58   = 383
  integer, parameter :: k_co59__he4_mn55   = 384
  integer, parameter :: k_ni56__p_co55   = 385
  integer, parameter :: k_ni56__he4_fe52   = 386
  integer, parameter :: k_ni57__n_ni56   = 387
  integer, parameter :: k_ni57__p_co56   = 388
  integer, parameter :: k_ni57__he4_fe53   = 389
  integer, parameter :: k_ni58__n_ni57   = 390
  integer, parameter :: k_ni58__p_co57   = 391
  integer, parameter :: k_ni58__he4_fe54   = 392
  integer, parameter :: k_ni59__n_ni58   = 393
  integer, parameter :: k_ni59__p_co58   = 394
  integer, parameter :: k_ni59__he4_fe55   = 395
  integer, parameter :: k_ni60__n_ni59   = 396
  integer, parameter :: k_ni60__p_co59   = 397
  integer, parameter :: k_ni60__he4_fe56   = 398
  integer, parameter :: k_ni61__n_ni60   = 399
  integer, parameter :: k_ni61__he4_fe57   = 400
  integer, parameter :: k_ni62__n_ni61   = 401
  integer, parameter :: k_ni62__he4_fe58   = 402
  integer, parameter :: k_ni63__n_ni62   = 403
  integer, parameter :: k_ni64__n_ni63   = 404
  integer, parameter :: k_cu57__p_ni56   = 405
  integer, parameter :: k_cu57__he4_co53   = 406
  integer, parameter :: k_cu58__n_cu57   = 407
  integer, parameter :: k_cu58__p_ni57   = 408
  integer, parameter :: k_cu58__he4_co54   = 409
  integer, parameter :: k_cu59__n_cu58   = 410
  integer, parameter :: k_cu59__p_ni58   = 411
  integer, parameter :: k_cu59__he4_co55   = 412
  integer, parameter :: k_cu60__n_cu59   = 413
  integer, parameter :: k_cu60__p_ni59   = 414
  integer, parameter :: k_cu60__he4_co56   = 415
  integer, parameter :: k_cu61__n_cu60   = 416
  integer, parameter :: k_cu61__p_ni60   = 417
  integer, parameter :: k_cu61__he4_co57   = 418
  integer, parameter :: k_cu62__n_cu61   = 419
  integer, parameter :: k_cu62__p_ni61   = 420
  integer, parameter :: k_cu62__he4_co58   = 421
  integer, parameter :: k_cu63__n_cu62   = 422
  integer, parameter :: k_cu63__p_ni62   = 423
  integer, parameter :: k_cu63__he4_co59   = 424
  integer, parameter :: k_cu64__n_cu63   = 425
  integer, parameter :: k_cu64__p_ni63   = 426
  integer, parameter :: k_cu65__n_cu64   = 427
  integer, parameter :: k_cu65__p_ni64   = 428
  integer, parameter :: k_zn59__p_cu58   = 429
  integer, parameter :: k_zn59__p_ni58__weak__wc12   = 430
  integer, parameter :: k_zn60__n_zn59   = 431
  integer, parameter :: k_zn60__p_cu59   = 432
  integer, parameter :: k_zn60__he4_ni56   = 433
  integer, parameter :: k_zn61__n_zn60   = 434
  integer, parameter :: k_zn61__p_cu60   = 435
  integer, parameter :: k_zn61__he4_ni57   = 436
  integer, parameter :: k_zn62__n_zn61   = 437
  integer, parameter :: k_zn62__p_cu61   = 438
  integer, parameter :: k_zn62__he4_ni58   = 439
  integer, parameter :: k_zn63__n_zn62   = 440
  integer, parameter :: k_zn63__p_cu62   = 441
  integer, parameter :: k_zn63__he4_ni59   = 442
  integer, parameter :: k_zn64__n_zn63   = 443
  integer, parameter :: k_zn64__p_cu63   = 444
  integer, parameter :: k_zn64__he4_ni60   = 445
  integer, parameter :: k_zn65__n_zn64   = 446
  integer, parameter :: k_zn65__p_cu64   = 447
  integer, parameter :: k_zn65__he4_ni61   = 448
  integer, parameter :: k_zn66__n_zn65   = 449
  integer, parameter :: k_zn66__p_cu65   = 450
  integer, parameter :: k_zn66__he4_ni62   = 451
  integer, parameter :: k_ga62__p_zn61   = 452
  integer, parameter :: k_ga62__he4_cu58   = 453
  integer, parameter :: k_ga63__n_ga62   = 454
  integer, parameter :: k_ga63__p_zn62   = 455
  integer, parameter :: k_ga63__he4_cu59   = 456
  integer, parameter :: k_ga64__n_ga63   = 457
  integer, parameter :: k_ga64__p_zn63   = 458
  integer, parameter :: k_ga64__he4_cu60   = 459
  integer, parameter :: k_ge63__p_ga62   = 460
  integer, parameter :: k_ge63__he4_zn59   = 461
  integer, parameter :: k_ge64__n_ge63   = 462
  integer, parameter :: k_ge64__p_ga63   = 463
  integer, parameter :: k_ge64__he4_zn60   = 464
  integer, parameter :: k_li6__n_p_he4   = 465
  integer, parameter :: k_be9__n_he4_he4   = 466
  integer, parameter :: k_c12__he4_he4_he4   = 467
  integer, parameter :: k_n_p__d   = 468
  integer, parameter :: k_p_p__d__weak__bet_pos_   = 469
  integer, parameter :: k_p_p__d__weak__electron_capture   = 470
  integer, parameter :: k_p_d__he3   = 471
  integer, parameter :: k_d_d__he4   = 472
  integer, parameter :: k_he4_d__li6   = 473
  integer, parameter :: k_n_he3__he4   = 474
  integer, parameter :: k_p_he3__he4__weak__bet_pos_   = 475
  integer, parameter :: k_he4_he3__be7   = 476
  integer, parameter :: k_n_li6__li7   = 477
  integer, parameter :: k_p_li6__be7   = 478
  integer, parameter :: k_he4_li6__b10   = 479
  integer, parameter :: k_he4_li7__b11   = 480
  integer, parameter :: k_p_be7__b8   = 481
  integer, parameter :: k_p_be9__b10   = 482
  integer, parameter :: k_n_b10__b11   = 483
  integer, parameter :: k_p_b11__c12   = 484
  integer, parameter :: k_n_c12__c13   = 485
  integer, parameter :: k_p_c12__n13   = 486
  integer, parameter :: k_he4_c12__o16   = 487
  integer, parameter :: k_n_c13__c14   = 488
  integer, parameter :: k_p_c13__n14   = 489
  integer, parameter :: k_p_c14__n15   = 490
  integer, parameter :: k_he4_c14__o18   = 491
  integer, parameter :: k_n_n13__n14   = 492
  integer, parameter :: k_p_n13__o14   = 493
  integer, parameter :: k_n_n14__n15   = 494
  integer, parameter :: k_p_n14__o15   = 495
  integer, parameter :: k_he4_n14__f18   = 496
  integer, parameter :: k_p_n15__o16   = 497
  integer, parameter :: k_he4_n15__f19   = 498
  integer, parameter :: k_n_o14__o15   = 499
  integer, parameter :: k_he4_o14__ne18   = 500
  integer, parameter :: k_n_o15__o16   = 501
  integer, parameter :: k_he4_o15__ne19   = 502
  integer, parameter :: k_n_o16__o17   = 503
  integer, parameter :: k_p_o16__f17   = 504
  integer, parameter :: k_he4_o16__ne20   = 505
  integer, parameter :: k_n_o17__o18   = 506
  integer, parameter :: k_p_o17__f18   = 507
  integer, parameter :: k_he4_o17__ne21   = 508
  integer, parameter :: k_p_o18__f19   = 509
  integer, parameter :: k_he4_o18__ne22   = 510
  integer, parameter :: k_n_f17__f18   = 511
  integer, parameter :: k_p_f17__ne18   = 512
  integer, parameter :: k_he4_f17__na21   = 513
  integer, parameter :: k_n_f18__f19   = 514
  integer, parameter :: k_p_f18__ne19   = 515
  integer, parameter :: k_he4_f18__na22   = 516
  integer, parameter :: k_p_f19__ne20   = 517
  integer, parameter :: k_he4_f19__na23   = 518
  integer, parameter :: k_n_ne18__ne19   = 519
  integer, parameter :: k_n_ne19__ne20   = 520
  integer, parameter :: k_he4_ne19__mg23   = 521
  integer, parameter :: k_n_ne20__ne21   = 522
  integer, parameter :: k_p_ne20__na21   = 523
  integer, parameter :: k_he4_ne20__mg24   = 524
  integer, parameter :: k_n_ne21__ne22   = 525
  integer, parameter :: k_p_ne21__na22   = 526
  integer, parameter :: k_he4_ne21__mg25   = 527
  integer, parameter :: k_p_ne22__na23   = 528
  integer, parameter :: k_he4_ne22__mg26   = 529
  integer, parameter :: k_n_na21__na22   = 530
  integer, parameter :: k_he4_na21__al25   = 531
  integer, parameter :: k_n_na22__na23   = 532
  integer, parameter :: k_p_na22__mg23   = 533
  integer, parameter :: k_he4_na22__al26   = 534
  integer, parameter :: k_p_na23__mg24   = 535
  integer, parameter :: k_he4_na23__al27   = 536
  integer, parameter :: k_n_mg23__mg24   = 537
  integer, parameter :: k_n_mg24__mg25   = 538
  integer, parameter :: k_p_mg24__al25   = 539
  integer, parameter :: k_he4_mg24__si28   = 540
  integer, parameter :: k_n_mg25__mg26   = 541
  integer, parameter :: k_p_mg25__al26   = 542
  integer, parameter :: k_he4_mg25__si29   = 543
  integer, parameter :: k_p_mg26__al27   = 544
  integer, parameter :: k_he4_mg26__si30   = 545
  integer, parameter :: k_n_al25__al26   = 546
  integer, parameter :: k_he4_al25__p29   = 547
  integer, parameter :: k_n_al26__al27   = 548
  integer, parameter :: k_he4_al26__p30   = 549
  integer, parameter :: k_p_al27__si28   = 550
  integer, parameter :: k_he4_al27__p31   = 551
  integer, parameter :: k_n_si28__si29   = 552
  integer, parameter :: k_p_si28__p29   = 553
  integer, parameter :: k_he4_si28__s32   = 554
  integer, parameter :: k_n_si29__si30   = 555
  integer, parameter :: k_p_si29__p30   = 556
  integer, parameter :: k_he4_si29__s33   = 557
  integer, parameter :: k_n_si30__si31   = 558
  integer, parameter :: k_p_si30__p31   = 559
  integer, parameter :: k_he4_si30__s34   = 560
  integer, parameter :: k_n_si31__si32   = 561
  integer, parameter :: k_p_si31__p32   = 562
  integer, parameter :: k_he4_si31__s35   = 563
  integer, parameter :: k_p_si32__p33   = 564
  integer, parameter :: k_he4_si32__s36   = 565
  integer, parameter :: k_n_p29__p30   = 566
  integer, parameter :: k_he4_p29__cl33   = 567
  integer, parameter :: k_n_p30__p31   = 568
  integer, parameter :: k_he4_p30__cl34   = 569
  integer, parameter :: k_n_p31__p32   = 570
  integer, parameter :: k_p_p31__s32   = 571
  integer, parameter :: k_he4_p31__cl35   = 572
  integer, parameter :: k_n_p32__p33   = 573
  integer, parameter :: k_p_p32__s33   = 574
  integer, parameter :: k_he4_p32__cl36   = 575
  integer, parameter :: k_p_p33__s34   = 576
  integer, parameter :: k_he4_p33__cl37   = 577
  integer, parameter :: k_n_s32__s33   = 578
  integer, parameter :: k_p_s32__cl33   = 579
  integer, parameter :: k_he4_s32__ar36   = 580
  integer, parameter :: k_n_s33__s34   = 581
  integer, parameter :: k_p_s33__cl34   = 582
  integer, parameter :: k_he4_s33__ar37   = 583
  integer, parameter :: k_n_s34__s35   = 584
  integer, parameter :: k_p_s34__cl35   = 585
  integer, parameter :: k_he4_s34__ar38   = 586
  integer, parameter :: k_n_s35__s36   = 587
  integer, parameter :: k_p_s35__cl36   = 588
  integer, parameter :: k_he4_s35__ar39   = 589
  integer, parameter :: k_p_s36__cl37   = 590
  integer, parameter :: k_he4_s36__ar40   = 591
  integer, parameter :: k_n_cl33__cl34   = 592
  integer, parameter :: k_he4_cl33__k37   = 593
  integer, parameter :: k_n_cl34__cl35   = 594
  integer, parameter :: k_he4_cl34__k38   = 595
  integer, parameter :: k_n_cl35__cl36   = 596
  integer, parameter :: k_p_cl35__ar36   = 597
  integer, parameter :: k_he4_cl35__k39   = 598
  integer, parameter :: k_n_cl36__cl37   = 599
  integer, parameter :: k_p_cl36__ar37   = 600
  integer, parameter :: k_he4_cl36__k40   = 601
  integer, parameter :: k_p_cl37__ar38   = 602
  integer, parameter :: k_he4_cl37__k41   = 603
  integer, parameter :: k_n_ar36__ar37   = 604
  integer, parameter :: k_p_ar36__k37   = 605
  integer, parameter :: k_he4_ar36__ca40   = 606
  integer, parameter :: k_n_ar37__ar38   = 607
  integer, parameter :: k_p_ar37__k38   = 608
  integer, parameter :: k_he4_ar37__ca41   = 609
  integer, parameter :: k_n_ar38__ar39   = 610
  integer, parameter :: k_p_ar38__k39   = 611
  integer, parameter :: k_he4_ar38__ca42   = 612
  integer, parameter :: k_n_ar39__ar40   = 613
  integer, parameter :: k_p_ar39__k40   = 614
  integer, parameter :: k_he4_ar39__ca43   = 615
  integer, parameter :: k_p_ar40__k41   = 616
  integer, parameter :: k_he4_ar40__ca44   = 617
  integer, parameter :: k_n_k37__k38   = 618
  integer, parameter :: k_n_k38__k39   = 619
  integer, parameter :: k_n_k39__k40   = 620
  integer, parameter :: k_p_k39__ca40   = 621
  integer, parameter :: k_he4_k39__sc43   = 622
  integer, parameter :: k_n_k40__k41   = 623
  integer, parameter :: k_p_k40__ca41   = 624
  integer, parameter :: k_he4_k40__sc44   = 625
  integer, parameter :: k_p_k41__ca42   = 626
  integer, parameter :: k_he4_k41__sc45   = 627
  integer, parameter :: k_n_ca40__ca41   = 628
  integer, parameter :: k_he4_ca40__ti44   = 629
  integer, parameter :: k_n_ca41__ca42   = 630
  integer, parameter :: k_he4_ca41__ti45   = 631
  integer, parameter :: k_n_ca42__ca43   = 632
  integer, parameter :: k_p_ca42__sc43   = 633
  integer, parameter :: k_he4_ca42__ti46   = 634
  integer, parameter :: k_n_ca43__ca44   = 635
  integer, parameter :: k_p_ca43__sc44   = 636
  integer, parameter :: k_he4_ca43__ti47   = 637
  integer, parameter :: k_n_ca44__ca45   = 638
  integer, parameter :: k_p_ca44__sc45   = 639
  integer, parameter :: k_he4_ca44__ti48   = 640
  integer, parameter :: k_n_ca45__ca46   = 641
  integer, parameter :: k_p_ca45__sc46   = 642
  integer, parameter :: k_he4_ca45__ti49   = 643
  integer, parameter :: k_n_ca46__ca47   = 644
  integer, parameter :: k_p_ca46__sc47   = 645
  integer, parameter :: k_he4_ca46__ti50   = 646
  integer, parameter :: k_n_ca47__ca48   = 647
  integer, parameter :: k_p_ca47__sc48   = 648
  integer, parameter :: k_he4_ca47__ti51   = 649
  integer, parameter :: k_p_ca48__sc49   = 650
  integer, parameter :: k_n_sc43__sc44   = 651
  integer, parameter :: k_p_sc43__ti44   = 652
  integer, parameter :: k_he4_sc43__v47   = 653
  integer, parameter :: k_n_sc44__sc45   = 654
  integer, parameter :: k_p_sc44__ti45   = 655
  integer, parameter :: k_he4_sc44__v48   = 656
  integer, parameter :: k_n_sc45__sc46   = 657
  integer, parameter :: k_p_sc45__ti46   = 658
  integer, parameter :: k_he4_sc45__v49   = 659
  integer, parameter :: k_n_sc46__sc47   = 660
  integer, parameter :: k_p_sc46__ti47   = 661
  integer, parameter :: k_he4_sc46__v50   = 662
  integer, parameter :: k_n_sc47__sc48   = 663
  integer, parameter :: k_p_sc47__ti48   = 664
  integer, parameter :: k_he4_sc47__v51   = 665
  integer, parameter :: k_n_sc48__sc49   = 666
  integer, parameter :: k_p_sc48__ti49   = 667
  integer, parameter :: k_he4_sc48__v52   = 668
  integer, parameter :: k_p_sc49__ti50   = 669
  integer, parameter :: k_n_ti44__ti45   = 670
  integer, parameter :: k_he4_ti44__cr48   = 671
  integer, parameter :: k_n_ti45__ti46   = 672
  integer, parameter :: k_p_ti45__v46   = 673
  integer, parameter :: k_he4_ti45__cr49   = 674
  integer, parameter :: k_n_ti46__ti47   = 675
  integer, parameter :: k_p_ti46__v47   = 676
  integer, parameter :: k_he4_ti46__cr50   = 677
  integer, parameter :: k_n_ti47__ti48   = 678
  integer, parameter :: k_p_ti47__v48   = 679
  integer, parameter :: k_he4_ti47__cr51   = 680
  integer, parameter :: k_n_ti48__ti49   = 681
  integer, parameter :: k_p_ti48__v49   = 682
  integer, parameter :: k_he4_ti48__cr52   = 683
  integer, parameter :: k_n_ti49__ti50   = 684
  integer, parameter :: k_p_ti49__v50   = 685
  integer, parameter :: k_he4_ti49__cr53   = 686
  integer, parameter :: k_n_ti50__ti51   = 687
  integer, parameter :: k_p_ti50__v51   = 688
  integer, parameter :: k_he4_ti50__cr54   = 689
  integer, parameter :: k_p_ti51__v52   = 690
  integer, parameter :: k_n_v46__v47   = 691
  integer, parameter :: k_he4_v46__mn50   = 692
  integer, parameter :: k_n_v47__v48   = 693
  integer, parameter :: k_p_v47__cr48   = 694
  integer, parameter :: k_he4_v47__mn51   = 695
  integer, parameter :: k_n_v48__v49   = 696
  integer, parameter :: k_p_v48__cr49   = 697
  integer, parameter :: k_he4_v48__mn52   = 698
  integer, parameter :: k_n_v49__v50   = 699
  integer, parameter :: k_p_v49__cr50   = 700
  integer, parameter :: k_he4_v49__mn53   = 701
  integer, parameter :: k_n_v50__v51   = 702
  integer, parameter :: k_p_v50__cr51   = 703
  integer, parameter :: k_he4_v50__mn54   = 704
  integer, parameter :: k_n_v51__v52   = 705
  integer, parameter :: k_p_v51__cr52   = 706
  integer, parameter :: k_he4_v51__mn55   = 707
  integer, parameter :: k_p_v52__cr53   = 708
  integer, parameter :: k_n_cr48__cr49   = 709
  integer, parameter :: k_he4_cr48__fe52   = 710
  integer, parameter :: k_n_cr49__cr50   = 711
  integer, parameter :: k_p_cr49__mn50   = 712
  integer, parameter :: k_he4_cr49__fe53   = 713
  integer, parameter :: k_n_cr50__cr51   = 714
  integer, parameter :: k_p_cr50__mn51   = 715
  integer, parameter :: k_he4_cr50__fe54   = 716
  integer, parameter :: k_n_cr51__cr52   = 717
  integer, parameter :: k_p_cr51__mn52   = 718
  integer, parameter :: k_he4_cr51__fe55   = 719
  integer, parameter :: k_n_cr52__cr53   = 720
  integer, parameter :: k_p_cr52__mn53   = 721
  integer, parameter :: k_he4_cr52__fe56   = 722
  integer, parameter :: k_n_cr53__cr54   = 723
  integer, parameter :: k_p_cr53__mn54   = 724
  integer, parameter :: k_he4_cr53__fe57   = 725
  integer, parameter :: k_p_cr54__mn55   = 726
  integer, parameter :: k_he4_cr54__fe58   = 727
  integer, parameter :: k_n_mn50__mn51   = 728
  integer, parameter :: k_he4_mn50__co54   = 729
  integer, parameter :: k_n_mn51__mn52   = 730
  integer, parameter :: k_p_mn51__fe52   = 731
  integer, parameter :: k_he4_mn51__co55   = 732
  integer, parameter :: k_n_mn52__mn53   = 733
  integer, parameter :: k_p_mn52__fe53   = 734
  integer, parameter :: k_he4_mn52__co56   = 735
  integer, parameter :: k_n_mn53__mn54   = 736
  integer, parameter :: k_p_mn53__fe54   = 737
  integer, parameter :: k_he4_mn53__co57   = 738
  integer, parameter :: k_n_mn54__mn55   = 739
  integer, parameter :: k_p_mn54__fe55   = 740
  integer, parameter :: k_he4_mn54__co58   = 741
  integer, parameter :: k_p_mn55__fe56   = 742
  integer, parameter :: k_he4_mn55__co59   = 743
  integer, parameter :: k_n_fe52__fe53   = 744
  integer, parameter :: k_p_fe52__co53   = 745
  integer, parameter :: k_he4_fe52__ni56   = 746
  integer, parameter :: k_n_fe53__fe54   = 747
  integer, parameter :: k_p_fe53__co54   = 748
  integer, parameter :: k_he4_fe53__ni57   = 749
  integer, parameter :: k_n_fe54__fe55   = 750
  integer, parameter :: k_p_fe54__co55   = 751
  integer, parameter :: k_he4_fe54__ni58   = 752
  integer, parameter :: k_n_fe55__fe56   = 753
  integer, parameter :: k_p_fe55__co56   = 754
  integer, parameter :: k_he4_fe55__ni59   = 755
  integer, parameter :: k_n_fe56__fe57   = 756
  integer, parameter :: k_p_fe56__co57   = 757
  integer, parameter :: k_he4_fe56__ni60   = 758
  integer, parameter :: k_n_fe57__fe58   = 759
  integer, parameter :: k_p_fe57__co58   = 760
  integer, parameter :: k_he4_fe57__ni61   = 761
  integer, parameter :: k_p_fe58__co59   = 762
  integer, parameter :: k_he4_fe58__ni62   = 763
  integer, parameter :: k_n_co53__co54   = 764
  integer, parameter :: k_he4_co53__cu57   = 765
  integer, parameter :: k_n_co54__co55   = 766
  integer, parameter :: k_he4_co54__cu58   = 767
  integer, parameter :: k_n_co55__co56   = 768
  integer, parameter :: k_p_co55__ni56   = 769
  integer, parameter :: k_he4_co55__cu59   = 770
  integer, parameter :: k_n_co56__co57   = 771
  integer, parameter :: k_p_co56__ni57   = 772
  integer, parameter :: k_he4_co56__cu60   = 773
  integer, parameter :: k_n_co57__co58   = 774
  integer, parameter :: k_p_co57__ni58   = 775
  integer, parameter :: k_he4_co57__cu61   = 776
  integer, parameter :: k_n_co58__co59   = 777
  integer, parameter :: k_p_co58__ni59   = 778
  integer, parameter :: k_he4_co58__cu62   = 779
  integer, parameter :: k_p_co59__ni60   = 780
  integer, parameter :: k_he4_co59__cu63   = 781
  integer, parameter :: k_n_ni56__ni57   = 782
  integer, parameter :: k_p_ni56__cu57   = 783
  integer, parameter :: k_he4_ni56__zn60   = 784
  integer, parameter :: k_n_ni57__ni58   = 785
  integer, parameter :: k_p_ni57__cu58   = 786
  integer, parameter :: k_he4_ni57__zn61   = 787
  integer, parameter :: k_n_ni58__ni59   = 788
  integer, parameter :: k_p_ni58__cu59   = 789
  integer, parameter :: k_he4_ni58__zn62   = 790
  integer, parameter :: k_n_ni59__ni60   = 791
  integer, parameter :: k_p_ni59__cu60   = 792
  integer, parameter :: k_he4_ni59__zn63   = 793
  integer, parameter :: k_n_ni60__ni61   = 794
  integer, parameter :: k_p_ni60__cu61   = 795
  integer, parameter :: k_he4_ni60__zn64   = 796
  integer, parameter :: k_n_ni61__ni62   = 797
  integer, parameter :: k_p_ni61__cu62   = 798
  integer, parameter :: k_he4_ni61__zn65   = 799
  integer, parameter :: k_n_ni62__ni63   = 800
  integer, parameter :: k_p_ni62__cu63   = 801
  integer, parameter :: k_he4_ni62__zn66   = 802
  integer, parameter :: k_n_ni63__ni64   = 803
  integer, parameter :: k_p_ni63__cu64   = 804
  integer, parameter :: k_p_ni64__cu65   = 805
  integer, parameter :: k_n_cu57__cu58   = 806
  integer, parameter :: k_n_cu58__cu59   = 807
  integer, parameter :: k_p_cu58__zn59   = 808
  integer, parameter :: k_he4_cu58__ga62   = 809
  integer, parameter :: k_n_cu59__cu60   = 810
  integer, parameter :: k_p_cu59__zn60   = 811
  integer, parameter :: k_he4_cu59__ga63   = 812
  integer, parameter :: k_n_cu60__cu61   = 813
  integer, parameter :: k_p_cu60__zn61   = 814
  integer, parameter :: k_he4_cu60__ga64   = 815
  integer, parameter :: k_n_cu61__cu62   = 816
  integer, parameter :: k_p_cu61__zn62   = 817
  integer, parameter :: k_n_cu62__cu63   = 818
  integer, parameter :: k_p_cu62__zn63   = 819
  integer, parameter :: k_n_cu63__cu64   = 820
  integer, parameter :: k_p_cu63__zn64   = 821
  integer, parameter :: k_n_cu64__cu65   = 822
  integer, parameter :: k_p_cu64__zn65   = 823
  integer, parameter :: k_p_cu65__zn66   = 824
  integer, parameter :: k_n_zn59__zn60   = 825
  integer, parameter :: k_he4_zn59__ge63   = 826
  integer, parameter :: k_n_zn60__zn61   = 827
  integer, parameter :: k_he4_zn60__ge64   = 828
  integer, parameter :: k_n_zn61__zn62   = 829
  integer, parameter :: k_p_zn61__ga62   = 830
  integer, parameter :: k_n_zn62__zn63   = 831
  integer, parameter :: k_p_zn62__ga63   = 832
  integer, parameter :: k_n_zn63__zn64   = 833
  integer, parameter :: k_p_zn63__ga64   = 834
  integer, parameter :: k_n_zn64__zn65   = 835
  integer, parameter :: k_n_zn65__zn66   = 836
  integer, parameter :: k_n_ga62__ga63   = 837
  integer, parameter :: k_p_ga62__ge63   = 838
  integer, parameter :: k_n_ga63__ga64   = 839
  integer, parameter :: k_p_ga63__ge64   = 840
  integer, parameter :: k_n_ge63__ge64   = 841
  integer, parameter :: k_d_d__n_he3   = 842
  integer, parameter :: k_n_he3__d_d   = 843
  integer, parameter :: k_d_he3__p_he4   = 844
  integer, parameter :: k_he4_he3__p_li6   = 845
  integer, parameter :: k_p_he4__d_he3   = 846
  integer, parameter :: k_he4_he4__n_be7   = 847
  integer, parameter :: k_he4_he4__p_li7   = 848
  integer, parameter :: k_p_li6__he4_he3   = 849
  integer, parameter :: k_d_li6__n_be7   = 850
  integer, parameter :: k_d_li6__p_li7   = 851
  integer, parameter :: k_he4_li6__p_be9   = 852
  integer, parameter :: k_p_li7__n_be7   = 853
  integer, parameter :: k_p_li7__d_li6   = 854
  integer, parameter :: k_p_li7__he4_he4   = 855
  integer, parameter :: k_he4_li7__n_b10   = 856
  integer, parameter :: k_n_be7__p_li7   = 857
  integer, parameter :: k_n_be7__d_li6   = 858
  integer, parameter :: k_n_be7__he4_he4   = 859
  integer, parameter :: k_he4_be7__p_b10   = 860
  integer, parameter :: k_p_be9__he4_li6   = 861
  integer, parameter :: k_he4_be9__n_c12   = 862
  integer, parameter :: k_n_b10__he4_li7   = 863
  integer, parameter :: k_p_b10__he4_be7   = 864
  integer, parameter :: k_he4_b10__n_n13   = 865
  integer, parameter :: k_he4_b10__p_c13   = 866
  integer, parameter :: k_he4_b11__n_n14   = 867
  integer, parameter :: k_he4_b11__p_c14   = 868
  integer, parameter :: k_n_c12__he4_be9   = 869
  integer, parameter :: k_he4_c12__n_o15   = 870
  integer, parameter :: k_he4_c12__p_n15   = 871
  integer, parameter :: k_c12_c12__n_mg23   = 872
  integer, parameter :: k_c12_c12__p_na23   = 873
  integer, parameter :: k_c12_c12__he4_ne20   = 874
  integer, parameter :: k_p_c13__n_n13   = 875
  integer, parameter :: k_p_c13__he4_b10   = 876
  integer, parameter :: k_d_c13__n_n14   = 877
  integer, parameter :: k_he4_c13__n_o16   = 878
  integer, parameter :: k_p_c14__n_n14   = 879
  integer, parameter :: k_p_c14__he4_b11   = 880
  integer, parameter :: k_d_c14__n_n15   = 881
  integer, parameter :: k_he4_c14__n_o17   = 882
  integer, parameter :: k_n_n13__p_c13   = 883
  integer, parameter :: k_n_n13__he4_b10   = 884
  integer, parameter :: k_he4_n13__p_o16   = 885
  integer, parameter :: k_n_n14__p_c14   = 886
  integer, parameter :: k_n_n14__d_c13   = 887
  integer, parameter :: k_n_n14__he4_b11   = 888
  integer, parameter :: k_p_n14__n_o14   = 889
  integer, parameter :: k_he4_n14__n_f17   = 890
  integer, parameter :: k_he4_n14__p_o17   = 891
  integer, parameter :: k_n_n15__d_c14   = 892
  integer, parameter :: k_p_n15__n_o15   = 893
  integer, parameter :: k_p_n15__he4_c12   = 894
  integer, parameter :: k_he4_n15__n_f18   = 895
  integer, parameter :: k_he4_n15__p_o18   = 896
  integer, parameter :: k_n_o14__p_n14   = 897
  integer, parameter :: k_he4_o14__p_f17   = 898
  integer, parameter :: k_n_o15__p_n15   = 899
  integer, parameter :: k_n_o15__he4_c12   = 900
  integer, parameter :: k_he4_o15__n_ne18   = 901
  integer, parameter :: k_he4_o15__p_f18   = 902
  integer, parameter :: k_n_o16__he4_c13   = 903
  integer, parameter :: k_p_o16__he4_n13   = 904
  integer, parameter :: k_he4_o16__n_ne19   = 905
  integer, parameter :: k_he4_o16__p_f19   = 906
  integer, parameter :: k_c12_o16__p_al27   = 907
  integer, parameter :: k_c12_o16__he4_mg24   = 908
  integer, parameter :: k_o16_o16__p_p31   = 909
  integer, parameter :: k_o16_o16__he4_si28   = 910
  integer, parameter :: k_n_o17__he4_c14   = 911
  integer, parameter :: k_p_o17__n_f17   = 912
  integer, parameter :: k_p_o17__he4_n14   = 913
  integer, parameter :: k_he4_o17__n_ne20   = 914
  integer, parameter :: k_p_o18__n_f18   = 915
  integer, parameter :: k_p_o18__he4_n15   = 916
  integer, parameter :: k_he4_o18__n_ne21   = 917
  integer, parameter :: k_n_f17__p_o17   = 918
  integer, parameter :: k_n_f17__he4_n14   = 919
  integer, parameter :: k_p_f17__he4_o14   = 920
  integer, parameter :: k_he4_f17__p_ne20   = 921
  integer, parameter :: k_n_f18__p_o18   = 922
  integer, parameter :: k_n_f18__he4_n15   = 923
  integer, parameter :: k_p_f18__n_ne18   = 924
  integer, parameter :: k_p_f18__he4_o15   = 925
  integer, parameter :: k_he4_f18__n_na21   = 926
  integer, parameter :: k_he4_f18__p_ne21   = 927
  integer, parameter :: k_p_f19__n_ne19   = 928
  integer, parameter :: k_p_f19__he4_o16   = 929
  integer, parameter :: k_he4_f19__n_na22   = 930
  integer, parameter :: k_he4_f19__p_ne22   = 931
  integer, parameter :: k_n_ne18__p_f18   = 932
  integer, parameter :: k_n_ne18__he4_o15   = 933
  integer, parameter :: k_he4_ne18__p_na21   = 934
  integer, parameter :: k_n_ne19__p_f19   = 935
  integer, parameter :: k_n_ne19__he4_o16   = 936
  integer, parameter :: k_he4_ne19__p_na22   = 937
  integer, parameter :: k_n_ne20__he4_o17   = 938
  integer, parameter :: k_p_ne20__he4_f17   = 939
  integer, parameter :: k_he4_ne20__n_mg23   = 940
  integer, parameter :: k_he4_ne20__p_na23   = 941
  integer, parameter :: k_he4_ne20__c12_c12   = 942
  integer, parameter :: k_c12_ne20__p_p31   = 943
  integer, parameter :: k_c12_ne20__he4_si28   = 944
  integer, parameter :: k_n_ne21__he4_o18   = 945
  integer, parameter :: k_p_ne21__n_na21   = 946
  integer, parameter :: k_p_ne21__he4_f18   = 947
  integer, parameter :: k_he4_ne21__n_mg24   = 948
  integer, parameter :: k_p_ne22__n_na22   = 949
  integer, parameter :: k_p_ne22__he4_f19   = 950
  integer, parameter :: k_he4_ne22__n_mg25   = 951
  integer, parameter :: k_n_na21__p_ne21   = 952
  integer, parameter :: k_n_na21__he4_f18   = 953
  integer, parameter :: k_p_na21__he4_ne18   = 954
  integer, parameter :: k_he4_na21__p_mg24   = 955
  integer, parameter :: k_n_na22__p_ne22   = 956
  integer, parameter :: k_n_na22__he4_f19   = 957
  integer, parameter :: k_p_na22__he4_ne19   = 958
  integer, parameter :: k_he4_na22__n_al25   = 959
  integer, parameter :: k_he4_na22__p_mg25   = 960
  integer, parameter :: k_p_na23__n_mg23   = 961
  integer, parameter :: k_p_na23__he4_ne20   = 962
  integer, parameter :: k_p_na23__c12_c12   = 963
  integer, parameter :: k_he4_na23__n_al26   = 964
  integer, parameter :: k_he4_na23__p_mg26   = 965
  integer, parameter :: k_n_mg23__p_na23   = 966
  integer, parameter :: k_n_mg23__he4_ne20   = 967
  integer, parameter :: k_n_mg23__c12_c12   = 968
  integer, parameter :: k_he4_mg23__p_al26   = 969
  integer, parameter :: k_n_mg24__he4_ne21   = 970
  integer, parameter :: k_p_mg24__he4_na21   = 971
  integer, parameter :: k_he4_mg24__p_al27   = 972
  integer, parameter :: k_he4_mg24__c12_o16   = 973
  integer, parameter :: k_n_mg25__he4_ne22   = 974
  integer, parameter :: k_p_mg25__n_al25   = 975
  integer, parameter :: k_p_mg25__he4_na22   = 976
  integer, parameter :: k_he4_mg25__n_si28   = 977
  integer, parameter :: k_p_mg26__n_al26   = 978
  integer, parameter :: k_p_mg26__he4_na23   = 979
  integer, parameter :: k_he4_mg26__n_si29   = 980
  integer, parameter :: k_n_al25__p_mg25   = 981
  integer, parameter :: k_n_al25__he4_na22   = 982
  integer, parameter :: k_he4_al25__p_si28   = 983
  integer, parameter :: k_n_al26__p_mg26   = 984
  integer, parameter :: k_n_al26__he4_na23   = 985
  integer, parameter :: k_p_al26__he4_mg23   = 986
  integer, parameter :: k_he4_al26__n_p29   = 987
  integer, parameter :: k_he4_al26__p_si29   = 988
  integer, parameter :: k_p_al27__he4_mg24   = 989
  integer, parameter :: k_p_al27__c12_o16   = 990
  integer, parameter :: k_he4_al27__n_p30   = 991
  integer, parameter :: k_he4_al27__p_si30   = 992
  integer, parameter :: k_n_si28__he4_mg25   = 993
  integer, parameter :: k_p_si28__he4_al25   = 994
  integer, parameter :: k_he4_si28__p_p31   = 995
  integer, parameter :: k_he4_si28__c12_ne20   = 996
  integer, parameter :: k_he4_si28__o16_o16   = 997
  integer, parameter :: k_n_si29__he4_mg26   = 998
  integer, parameter :: k_p_si29__n_p29   = 999
  integer, parameter :: k_p_si29__he4_al26   = 1000
  integer, parameter :: k_he4_si29__n_s32   = 1001
  integer, parameter :: k_he4_si29__p_p32   = 1002
  integer, parameter :: k_p_si30__n_p30   = 1003
  integer, parameter :: k_p_si30__he4_al27   = 1004
  integer, parameter :: k_he4_si30__n_s33   = 1005
  integer, parameter :: k_he4_si30__p_p33   = 1006
  integer, parameter :: k_p_si31__n_p31   = 1007
  integer, parameter :: k_he4_si31__n_s34   = 1008
  integer, parameter :: k_p_si32__n_p32   = 1009
  integer, parameter :: k_he4_si32__n_s35   = 1010
  integer, parameter :: k_n_p29__p_si29   = 1011
  integer, parameter :: k_n_p29__he4_al26   = 1012
  integer, parameter :: k_he4_p29__p_s32   = 1013
  integer, parameter :: k_n_p30__p_si30   = 1014
  integer, parameter :: k_n_p30__he4_al27   = 1015
  integer, parameter :: k_he4_p30__n_cl33   = 1016
  integer, parameter :: k_he4_p30__p_s33   = 1017
  integer, parameter :: k_n_p31__p_si31   = 1018
  integer, parameter :: k_p_p31__he4_si28   = 1019
  integer, parameter :: k_p_p31__c12_ne20   = 1020
  integer, parameter :: k_p_p31__o16_o16   = 1021
  integer, parameter :: k_he4_p31__n_cl34   = 1022
  integer, parameter :: k_he4_p31__p_s34   = 1023
  integer, parameter :: k_n_p32__p_si32   = 1024
  integer, parameter :: k_p_p32__n_s32   = 1025
  integer, parameter :: k_p_p32__he4_si29   = 1026
  integer, parameter :: k_he4_p32__n_cl35   = 1027
  integer, parameter :: k_he4_p32__p_s35   = 1028
  integer, parameter :: k_p_p33__n_s33   = 1029
  integer, parameter :: k_p_p33__he4_si30   = 1030
  integer, parameter :: k_he4_p33__n_cl36   = 1031
  integer, parameter :: k_he4_p33__p_s36   = 1032
  integer, parameter :: k_n_s32__p_p32   = 1033
  integer, parameter :: k_n_s32__he4_si29   = 1034
  integer, parameter :: k_p_s32__he4_p29   = 1035
  integer, parameter :: k_he4_s32__p_cl35   = 1036
  integer, parameter :: k_n_s33__p_p33   = 1037
  integer, parameter :: k_n_s33__he4_si30   = 1038
  integer, parameter :: k_p_s33__n_cl33   = 1039
  integer, parameter :: k_p_s33__he4_p30   = 1040
  integer, parameter :: k_he4_s33__n_ar36   = 1041
  integer, parameter :: k_he4_s33__p_cl36   = 1042
  integer, parameter :: k_n_s34__he4_si31   = 1043
  integer, parameter :: k_p_s34__n_cl34   = 1044
  integer, parameter :: k_p_s34__he4_p31   = 1045
  integer, parameter :: k_he4_s34__n_ar37   = 1046
  integer, parameter :: k_he4_s34__p_cl37   = 1047
  integer, parameter :: k_n_s35__he4_si32   = 1048
  integer, parameter :: k_p_s35__n_cl35   = 1049
  integer, parameter :: k_p_s35__he4_p32   = 1050
  integer, parameter :: k_he4_s35__n_ar38   = 1051
  integer, parameter :: k_p_s36__n_cl36   = 1052
  integer, parameter :: k_p_s36__he4_p33   = 1053
  integer, parameter :: k_he4_s36__n_ar39   = 1054
  integer, parameter :: k_n_cl33__p_s33   = 1055
  integer, parameter :: k_n_cl33__he4_p30   = 1056
  integer, parameter :: k_he4_cl33__p_ar36   = 1057
  integer, parameter :: k_n_cl34__p_s34   = 1058
  integer, parameter :: k_n_cl34__he4_p31   = 1059
  integer, parameter :: k_he4_cl34__n_k37   = 1060
  integer, parameter :: k_he4_cl34__p_ar37   = 1061
  integer, parameter :: k_n_cl35__p_s35   = 1062
  integer, parameter :: k_n_cl35__he4_p32   = 1063
  integer, parameter :: k_p_cl35__he4_s32   = 1064
  integer, parameter :: k_he4_cl35__n_k38   = 1065
  integer, parameter :: k_he4_cl35__p_ar38   = 1066
  integer, parameter :: k_n_cl36__p_s36   = 1067
  integer, parameter :: k_n_cl36__he4_p33   = 1068
  integer, parameter :: k_p_cl36__n_ar36   = 1069
  integer, parameter :: k_p_cl36__he4_s33   = 1070
  integer, parameter :: k_he4_cl36__n_k39   = 1071
  integer, parameter :: k_he4_cl36__p_ar39   = 1072
  integer, parameter :: k_p_cl37__n_ar37   = 1073
  integer, parameter :: k_p_cl37__he4_s34   = 1074
  integer, parameter :: k_he4_cl37__n_k40   = 1075
  integer, parameter :: k_he4_cl37__p_ar40   = 1076
  integer, parameter :: k_n_ar36__p_cl36   = 1077
  integer, parameter :: k_n_ar36__he4_s33   = 1078
  integer, parameter :: k_p_ar36__he4_cl33   = 1079
  integer, parameter :: k_he4_ar36__p_k39   = 1080
  integer, parameter :: k_n_ar37__p_cl37   = 1081
  integer, parameter :: k_n_ar37__he4_s34   = 1082
  integer, parameter :: k_p_ar37__n_k37   = 1083
  integer, parameter :: k_p_ar37__he4_cl34   = 1084
  integer, parameter :: k_he4_ar37__n_ca40   = 1085
  integer, parameter :: k_he4_ar37__p_k40   = 1086
  integer, parameter :: k_n_ar38__he4_s35   = 1087
  integer, parameter :: k_p_ar38__n_k38   = 1088
  integer, parameter :: k_p_ar38__he4_cl35   = 1089
  integer, parameter :: k_he4_ar38__n_ca41   = 1090
  integer, parameter :: k_he4_ar38__p_k41   = 1091
  integer, parameter :: k_n_ar39__he4_s36   = 1092
  integer, parameter :: k_p_ar39__n_k39   = 1093
  integer, parameter :: k_p_ar39__he4_cl36   = 1094
  integer, parameter :: k_he4_ar39__n_ca42   = 1095
  integer, parameter :: k_p_ar40__n_k40   = 1096
  integer, parameter :: k_p_ar40__he4_cl37   = 1097
  integer, parameter :: k_he4_ar40__n_ca43   = 1098
  integer, parameter :: k_n_k37__p_ar37   = 1099
  integer, parameter :: k_n_k37__he4_cl34   = 1100
  integer, parameter :: k_he4_k37__p_ca40   = 1101
  integer, parameter :: k_n_k38__p_ar38   = 1102
  integer, parameter :: k_n_k38__he4_cl35   = 1103
  integer, parameter :: k_he4_k38__p_ca41   = 1104
  integer, parameter :: k_n_k39__p_ar39   = 1105
  integer, parameter :: k_n_k39__he4_cl36   = 1106
  integer, parameter :: k_p_k39__he4_ar36   = 1107
  integer, parameter :: k_he4_k39__p_ca42   = 1108
  integer, parameter :: k_n_k40__p_ar40   = 1109
  integer, parameter :: k_n_k40__he4_cl37   = 1110
  integer, parameter :: k_p_k40__n_ca40   = 1111
  integer, parameter :: k_p_k40__he4_ar37   = 1112
  integer, parameter :: k_he4_k40__n_sc43   = 1113
  integer, parameter :: k_he4_k40__p_ca43   = 1114
  integer, parameter :: k_p_k41__n_ca41   = 1115
  integer, parameter :: k_p_k41__he4_ar38   = 1116
  integer, parameter :: k_he4_k41__n_sc44   = 1117
  integer, parameter :: k_he4_k41__p_ca44   = 1118
  integer, parameter :: k_n_ca40__p_k40   = 1119
  integer, parameter :: k_n_ca40__he4_ar37   = 1120
  integer, parameter :: k_p_ca40__he4_k37   = 1121
  integer, parameter :: k_he4_ca40__p_sc43   = 1122
  integer, parameter :: k_n_ca41__p_k41   = 1123
  integer, parameter :: k_n_ca41__he4_ar38   = 1124
  integer, parameter :: k_p_ca41__he4_k38   = 1125
  integer, parameter :: k_he4_ca41__n_ti44   = 1126
  integer, parameter :: k_he4_ca41__p_sc44   = 1127
  integer, parameter :: k_n_ca42__he4_ar39   = 1128
  integer, parameter :: k_p_ca42__he4_k39   = 1129
  integer, parameter :: k_he4_ca42__n_ti45   = 1130
  integer, parameter :: k_he4_ca42__p_sc45   = 1131
  integer, parameter :: k_n_ca43__he4_ar40   = 1132
  integer, parameter :: k_p_ca43__n_sc43   = 1133
  integer, parameter :: k_p_ca43__he4_k40   = 1134
  integer, parameter :: k_he4_ca43__n_ti46   = 1135
  integer, parameter :: k_he4_ca43__p_sc46   = 1136
  integer, parameter :: k_p_ca44__n_sc44   = 1137
  integer, parameter :: k_p_ca44__he4_k41   = 1138
  integer, parameter :: k_he4_ca44__n_ti47   = 1139
  integer, parameter :: k_he4_ca44__p_sc47   = 1140
  integer, parameter :: k_p_ca45__n_sc45   = 1141
  integer, parameter :: k_he4_ca45__n_ti48   = 1142
  integer, parameter :: k_he4_ca45__p_sc48   = 1143
  integer, parameter :: k_p_ca46__n_sc46   = 1144
  integer, parameter :: k_he4_ca46__n_ti49   = 1145
  integer, parameter :: k_he4_ca46__p_sc49   = 1146
  integer, parameter :: k_p_ca47__n_sc47   = 1147
  integer, parameter :: k_he4_ca47__n_ti50   = 1148
  integer, parameter :: k_p_ca48__n_sc48   = 1149
  integer, parameter :: k_he4_ca48__n_ti51   = 1150
  integer, parameter :: k_n_sc43__p_ca43   = 1151
  integer, parameter :: k_n_sc43__he4_k40   = 1152
  integer, parameter :: k_p_sc43__he4_ca40   = 1153
  integer, parameter :: k_he4_sc43__n_v46   = 1154
  integer, parameter :: k_he4_sc43__p_ti46   = 1155
  integer, parameter :: k_n_sc44__p_ca44   = 1156
  integer, parameter :: k_n_sc44__he4_k41   = 1157
  integer, parameter :: k_p_sc44__n_ti44   = 1158
  integer, parameter :: k_p_sc44__he4_ca41   = 1159
  integer, parameter :: k_he4_sc44__n_v47   = 1160
  integer, parameter :: k_he4_sc44__p_ti47   = 1161
  integer, parameter :: k_n_sc45__p_ca45   = 1162
  integer, parameter :: k_p_sc45__n_ti45   = 1163
  integer, parameter :: k_p_sc45__he4_ca42   = 1164
  integer, parameter :: k_he4_sc45__n_v48   = 1165
  integer, parameter :: k_he4_sc45__p_ti48   = 1166
  integer, parameter :: k_n_sc46__p_ca46   = 1167
  integer, parameter :: k_p_sc46__n_ti46   = 1168
  integer, parameter :: k_p_sc46__he4_ca43   = 1169
  integer, parameter :: k_he4_sc46__n_v49   = 1170
  integer, parameter :: k_he4_sc46__p_ti49   = 1171
  integer, parameter :: k_n_sc47__p_ca47   = 1172
  integer, parameter :: k_p_sc47__n_ti47   = 1173
  integer, parameter :: k_p_sc47__he4_ca44   = 1174
  integer, parameter :: k_he4_sc47__n_v50   = 1175
  integer, parameter :: k_he4_sc47__p_ti50   = 1176
  integer, parameter :: k_n_sc48__p_ca48   = 1177
  integer, parameter :: k_p_sc48__n_ti48   = 1178
  integer, parameter :: k_p_sc48__he4_ca45   = 1179
  integer, parameter :: k_he4_sc48__n_v51   = 1180
  integer, parameter :: k_he4_sc48__p_ti51   = 1181
  integer, parameter :: k_p_sc49__n_ti49   = 1182
  integer, parameter :: k_p_sc49__he4_ca46   = 1183
  integer, parameter :: k_he4_sc49__n_v52   = 1184
  integer, parameter :: k_n_ti44__p_sc44   = 1185
  integer, parameter :: k_n_ti44__he4_ca41   = 1186
  integer, parameter :: k_he4_ti44__p_v47   = 1187
  integer, parameter :: k_n_ti45__p_sc45   = 1188
  integer, parameter :: k_n_ti45__he4_ca42   = 1189
  integer, parameter :: k_he4_ti45__n_cr48   = 1190
  integer, parameter :: k_he4_ti45__p_v48   = 1191
  integer, parameter :: k_n_ti46__p_sc46   = 1192
  integer, parameter :: k_n_ti46__he4_ca43   = 1193
  integer, parameter :: k_p_ti46__n_v46   = 1194
  integer, parameter :: k_p_ti46__he4_sc43   = 1195
  integer, parameter :: k_he4_ti46__n_cr49   = 1196
  integer, parameter :: k_he4_ti46__p_v49   = 1197
  integer, parameter :: k_n_ti47__p_sc47   = 1198
  integer, parameter :: k_n_ti47__he4_ca44   = 1199
  integer, parameter :: k_p_ti47__n_v47   = 1200
  integer, parameter :: k_p_ti47__he4_sc44   = 1201
  integer, parameter :: k_he4_ti47__n_cr50   = 1202
  integer, parameter :: k_he4_ti47__p_v50   = 1203
  integer, parameter :: k_n_ti48__p_sc48   = 1204
  integer, parameter :: k_n_ti48__he4_ca45   = 1205
  integer, parameter :: k_p_ti48__n_v48   = 1206
  integer, parameter :: k_p_ti48__he4_sc45   = 1207
  integer, parameter :: k_he4_ti48__n_cr51   = 1208
  integer, parameter :: k_he4_ti48__p_v51   = 1209
  integer, parameter :: k_n_ti49__p_sc49   = 1210
  integer, parameter :: k_n_ti49__he4_ca46   = 1211
  integer, parameter :: k_p_ti49__n_v49   = 1212
  integer, parameter :: k_p_ti49__he4_sc46   = 1213
  integer, parameter :: k_he4_ti49__n_cr52   = 1214
  integer, parameter :: k_he4_ti49__p_v52   = 1215
  integer, parameter :: k_n_ti50__he4_ca47   = 1216
  integer, parameter :: k_p_ti50__n_v50   = 1217
  integer, parameter :: k_p_ti50__he4_sc47   = 1218
  integer, parameter :: k_he4_ti50__n_cr53   = 1219
  integer, parameter :: k_n_ti51__he4_ca48   = 1220
  integer, parameter :: k_p_ti51__n_v51   = 1221
  integer, parameter :: k_p_ti51__he4_sc48   = 1222
  integer, parameter :: k_he4_ti51__n_cr54   = 1223
  integer, parameter :: k_n_v46__p_ti46   = 1224
  integer, parameter :: k_n_v46__he4_sc43   = 1225
  integer, parameter :: k_he4_v46__p_cr49   = 1226
  integer, parameter :: k_n_v47__p_ti47   = 1227
  integer, parameter :: k_n_v47__he4_sc44   = 1228
  integer, parameter :: k_p_v47__he4_ti44   = 1229
  integer, parameter :: k_he4_v47__n_mn50   = 1230
  integer, parameter :: k_he4_v47__p_cr50   = 1231
  integer, parameter :: k_n_v48__p_ti48   = 1232
  integer, parameter :: k_n_v48__he4_sc45   = 1233
  integer, parameter :: k_p_v48__n_cr48   = 1234
  integer, parameter :: k_p_v48__he4_ti45   = 1235
  integer, parameter :: k_he4_v48__n_mn51   = 1236
  integer, parameter :: k_he4_v48__p_cr51   = 1237
  integer, parameter :: k_n_v49__p_ti49   = 1238
  integer, parameter :: k_n_v49__he4_sc46   = 1239
  integer, parameter :: k_p_v49__n_cr49   = 1240
  integer, parameter :: k_p_v49__he4_ti46   = 1241
  integer, parameter :: k_he4_v49__n_mn52   = 1242
  integer, parameter :: k_he4_v49__p_cr52   = 1243
  integer, parameter :: k_n_v50__p_ti50   = 1244
  integer, parameter :: k_n_v50__he4_sc47   = 1245
  integer, parameter :: k_p_v50__n_cr50   = 1246
  integer, parameter :: k_p_v50__he4_ti47   = 1247
  integer, parameter :: k_he4_v50__n_mn53   = 1248
  integer, parameter :: k_he4_v50__p_cr53   = 1249
  integer, parameter :: k_n_v51__p_ti51   = 1250
  integer, parameter :: k_n_v51__he4_sc48   = 1251
  integer, parameter :: k_p_v51__n_cr51   = 1252
  integer, parameter :: k_p_v51__he4_ti48   = 1253
  integer, parameter :: k_he4_v51__n_mn54   = 1254
  integer, parameter :: k_he4_v51__p_cr54   = 1255
  integer, parameter :: k_n_v52__he4_sc49   = 1256
  integer, parameter :: k_p_v52__n_cr52   = 1257
  integer, parameter :: k_p_v52__he4_ti49   = 1258
  integer, parameter :: k_he4_v52__n_mn55   = 1259
  integer, parameter :: k_n_cr48__p_v48   = 1260
  integer, parameter :: k_n_cr48__he4_ti45   = 1261
  integer, parameter :: k_he4_cr48__p_mn51   = 1262
  integer, parameter :: k_n_cr49__p_v49   = 1263
  integer, parameter :: k_n_cr49__he4_ti46   = 1264
  integer, parameter :: k_p_cr49__he4_v46   = 1265
  integer, parameter :: k_he4_cr49__n_fe52   = 1266
  integer, parameter :: k_he4_cr49__p_mn52   = 1267
  integer, parameter :: k_n_cr50__p_v50   = 1268
  integer, parameter :: k_n_cr50__he4_ti47   = 1269
  integer, parameter :: k_p_cr50__n_mn50   = 1270
  integer, parameter :: k_p_cr50__he4_v47   = 1271
  integer, parameter :: k_he4_cr50__n_fe53   = 1272
  integer, parameter :: k_he4_cr50__p_mn53   = 1273
  integer, parameter :: k_n_cr51__p_v51   = 1274
  integer, parameter :: k_n_cr51__he4_ti48   = 1275
  integer, parameter :: k_p_cr51__n_mn51   = 1276
  integer, parameter :: k_p_cr51__he4_v48   = 1277
  integer, parameter :: k_he4_cr51__n_fe54   = 1278
  integer, parameter :: k_he4_cr51__p_mn54   = 1279
  integer, parameter :: k_n_cr52__p_v52   = 1280
  integer, parameter :: k_n_cr52__he4_ti49   = 1281
  integer, parameter :: k_p_cr52__n_mn52   = 1282
  integer, parameter :: k_p_cr52__he4_v49   = 1283
  integer, parameter :: k_he4_cr52__n_fe55   = 1284
  integer, parameter :: k_he4_cr52__p_mn55   = 1285
  integer, parameter :: k_n_cr53__he4_ti50   = 1286
  integer, parameter :: k_p_cr53__n_mn53   = 1287
  integer, parameter :: k_p_cr53__he4_v50   = 1288
  integer, parameter :: k_he4_cr53__n_fe56   = 1289
  integer, parameter :: k_n_cr54__he4_ti51   = 1290
  integer, parameter :: k_p_cr54__n_mn54   = 1291
  integer, parameter :: k_p_cr54__he4_v51   = 1292
  integer, parameter :: k_he4_cr54__n_fe57   = 1293
  integer, parameter :: k_n_mn50__p_cr50   = 1294
  integer, parameter :: k_n_mn50__he4_v47   = 1295
  integer, parameter :: k_he4_mn50__n_co53   = 1296
  integer, parameter :: k_he4_mn50__p_fe53   = 1297
  integer, parameter :: k_n_mn51__p_cr51   = 1298
  integer, parameter :: k_n_mn51__he4_v48   = 1299
  integer, parameter :: k_p_mn51__he4_cr48   = 1300
  integer, parameter :: k_he4_mn51__n_co54   = 1301
  integer, parameter :: k_he4_mn51__p_fe54   = 1302
  integer, parameter :: k_n_mn52__p_cr52   = 1303
  integer, parameter :: k_n_mn52__he4_v49   = 1304
  integer, parameter :: k_p_mn52__n_fe52   = 1305
  integer, parameter :: k_p_mn52__he4_cr49   = 1306
  integer, parameter :: k_he4_mn52__n_co55   = 1307
  integer, parameter :: k_he4_mn52__p_fe55   = 1308
  integer, parameter :: k_n_mn53__p_cr53   = 1309
  integer, parameter :: k_n_mn53__he4_v50   = 1310
  integer, parameter :: k_p_mn53__n_fe53   = 1311
  integer, parameter :: k_p_mn53__he4_cr50   = 1312
  integer, parameter :: k_he4_mn53__n_co56   = 1313
  integer, parameter :: k_he4_mn53__p_fe56   = 1314
  integer, parameter :: k_n_mn54__p_cr54   = 1315
  integer, parameter :: k_n_mn54__he4_v51   = 1316
  integer, parameter :: k_p_mn54__n_fe54   = 1317
  integer, parameter :: k_p_mn54__he4_cr51   = 1318
  integer, parameter :: k_he4_mn54__n_co57   = 1319
  integer, parameter :: k_he4_mn54__p_fe57   = 1320
  integer, parameter :: k_n_mn55__he4_v52   = 1321
  integer, parameter :: k_p_mn55__n_fe55   = 1322
  integer, parameter :: k_p_mn55__he4_cr52   = 1323
  integer, parameter :: k_he4_mn55__n_co58   = 1324
  integer, parameter :: k_he4_mn55__p_fe58   = 1325
  integer, parameter :: k_n_fe52__p_mn52   = 1326
  integer, parameter :: k_n_fe52__he4_cr49   = 1327
  integer, parameter :: k_he4_fe52__p_co55   = 1328
  integer, parameter :: k_n_fe53__p_mn53   = 1329
  integer, parameter :: k_n_fe53__he4_cr50   = 1330
  integer, parameter :: k_p_fe53__n_co53   = 1331
  integer, parameter :: k_p_fe53__he4_mn50   = 1332
  integer, parameter :: k_he4_fe53__n_ni56   = 1333
  integer, parameter :: k_he4_fe53__p_co56   = 1334
  integer, parameter :: k_n_fe54__p_mn54   = 1335
  integer, parameter :: k_n_fe54__he4_cr51   = 1336
  integer, parameter :: k_p_fe54__n_co54   = 1337
  integer, parameter :: k_p_fe54__he4_mn51   = 1338
  integer, parameter :: k_he4_fe54__n_ni57   = 1339
  integer, parameter :: k_he4_fe54__p_co57   = 1340
  integer, parameter :: k_n_fe55__p_mn55   = 1341
  integer, parameter :: k_n_fe55__he4_cr52   = 1342
  integer, parameter :: k_p_fe55__n_co55   = 1343
  integer, parameter :: k_p_fe55__he4_mn52   = 1344
  integer, parameter :: k_he4_fe55__n_ni58   = 1345
  integer, parameter :: k_he4_fe55__p_co58   = 1346
  integer, parameter :: k_n_fe56__he4_cr53   = 1347
  integer, parameter :: k_p_fe56__n_co56   = 1348
  integer, parameter :: k_p_fe56__he4_mn53   = 1349
  integer, parameter :: k_he4_fe56__n_ni59   = 1350
  integer, parameter :: k_he4_fe56__p_co59   = 1351
  integer, parameter :: k_n_fe57__he4_cr54   = 1352
  integer, parameter :: k_p_fe57__n_co57   = 1353
  integer, parameter :: k_p_fe57__he4_mn54   = 1354
  integer, parameter :: k_he4_fe57__n_ni60   = 1355
  integer, parameter :: k_p_fe58__n_co58   = 1356
  integer, parameter :: k_p_fe58__he4_mn55   = 1357
  integer, parameter :: k_he4_fe58__n_ni61   = 1358
  integer, parameter :: k_n_co53__p_fe53   = 1359
  integer, parameter :: k_n_co53__he4_mn50   = 1360
  integer, parameter :: k_he4_co53__p_ni56   = 1361
  integer, parameter :: k_n_co54__p_fe54   = 1362
  integer, parameter :: k_n_co54__he4_mn51   = 1363
  integer, parameter :: k_he4_co54__n_cu57   = 1364
  integer, parameter :: k_he4_co54__p_ni57   = 1365
  integer, parameter :: k_n_co55__p_fe55   = 1366
  integer, parameter :: k_n_co55__he4_mn52   = 1367
  integer, parameter :: k_p_co55__he4_fe52   = 1368
  integer, parameter :: k_he4_co55__n_cu58   = 1369
  integer, parameter :: k_he4_co55__p_ni58   = 1370
  integer, parameter :: k_n_co56__p_fe56   = 1371
  integer, parameter :: k_n_co56__he4_mn53   = 1372
  integer, parameter :: k_p_co56__n_ni56   = 1373
  integer, parameter :: k_p_co56__he4_fe53   = 1374
  integer, parameter :: k_he4_co56__n_cu59   = 1375
  integer, parameter :: k_he4_co56__p_ni59   = 1376
  integer, parameter :: k_n_co57__p_fe57   = 1377
  integer, parameter :: k_n_co57__he4_mn54   = 1378
  integer, parameter :: k_p_co57__n_ni57   = 1379
  integer, parameter :: k_p_co57__he4_fe54   = 1380
  integer, parameter :: k_he4_co57__n_cu60   = 1381
  integer, parameter :: k_he4_co57__p_ni60   = 1382
  integer, parameter :: k_n_co58__p_fe58   = 1383
  integer, parameter :: k_n_co58__he4_mn55   = 1384
  integer, parameter :: k_p_co58__n_ni58   = 1385
  integer, parameter :: k_p_co58__he4_fe55   = 1386
  integer, parameter :: k_he4_co58__n_cu61   = 1387
  integer, parameter :: k_he4_co58__p_ni61   = 1388
  integer, parameter :: k_p_co59__n_ni59   = 1389
  integer, parameter :: k_p_co59__he4_fe56   = 1390
  integer, parameter :: k_he4_co59__n_cu62   = 1391
  integer, parameter :: k_he4_co59__p_ni62   = 1392
  integer, parameter :: k_n_ni56__p_co56   = 1393
  integer, parameter :: k_n_ni56__he4_fe53   = 1394
  integer, parameter :: k_p_ni56__he4_co53   = 1395
  integer, parameter :: k_he4_ni56__n_zn59   = 1396
  integer, parameter :: k_he4_ni56__p_cu59   = 1397
  integer, parameter :: k_n_ni57__p_co57   = 1398
  integer, parameter :: k_n_ni57__he4_fe54   = 1399
  integer, parameter :: k_p_ni57__n_cu57   = 1400
  integer, parameter :: k_p_ni57__he4_co54   = 1401
  integer, parameter :: k_he4_ni57__n_zn60   = 1402
  integer, parameter :: k_he4_ni57__p_cu60   = 1403
  integer, parameter :: k_n_ni58__p_co58   = 1404
  integer, parameter :: k_n_ni58__he4_fe55   = 1405
  integer, parameter :: k_p_ni58__n_cu58   = 1406
  integer, parameter :: k_p_ni58__he4_co55   = 1407
  integer, parameter :: k_he4_ni58__n_zn61   = 1408
  integer, parameter :: k_he4_ni58__p_cu61   = 1409
  integer, parameter :: k_n_ni59__p_co59   = 1410
  integer, parameter :: k_n_ni59__he4_fe56   = 1411
  integer, parameter :: k_p_ni59__n_cu59   = 1412
  integer, parameter :: k_p_ni59__he4_co56   = 1413
  integer, parameter :: k_he4_ni59__n_zn62   = 1414
  integer, parameter :: k_he4_ni59__p_cu62   = 1415
  integer, parameter :: k_n_ni60__he4_fe57   = 1416
  integer, parameter :: k_p_ni60__n_cu60   = 1417
  integer, parameter :: k_p_ni60__he4_co57   = 1418
  integer, parameter :: k_he4_ni60__n_zn63   = 1419
  integer, parameter :: k_he4_ni60__p_cu63   = 1420
  integer, parameter :: k_n_ni61__he4_fe58   = 1421
  integer, parameter :: k_p_ni61__n_cu61   = 1422
  integer, parameter :: k_p_ni61__he4_co58   = 1423
  integer, parameter :: k_he4_ni61__n_zn64   = 1424
  integer, parameter :: k_he4_ni61__p_cu64   = 1425
  integer, parameter :: k_p_ni62__n_cu62   = 1426
  integer, parameter :: k_p_ni62__he4_co59   = 1427
  integer, parameter :: k_he4_ni62__n_zn65   = 1428
  integer, parameter :: k_he4_ni62__p_cu65   = 1429
  integer, parameter :: k_p_ni63__n_cu63   = 1430
  integer, parameter :: k_he4_ni63__n_zn66   = 1431
  integer, parameter :: k_p_ni64__n_cu64   = 1432
  integer, parameter :: k_n_cu57__p_ni57   = 1433
  integer, parameter :: k_n_cu57__he4_co54   = 1434
  integer, parameter :: k_he4_cu57__p_zn60   = 1435
  integer, parameter :: k_n_cu58__p_ni58   = 1436
  integer, parameter :: k_n_cu58__he4_co55   = 1437
  integer, parameter :: k_he4_cu58__p_zn61   = 1438
  integer, parameter :: k_n_cu59__p_ni59   = 1439
  integer, parameter :: k_n_cu59__he4_co56   = 1440
  integer, parameter :: k_p_cu59__n_zn59   = 1441
  integer, parameter :: k_p_cu59__he4_ni56   = 1442
  integer, parameter :: k_he4_cu59__n_ga62   = 1443
  integer, parameter :: k_he4_cu59__p_zn62   = 1444
  integer, parameter :: k_n_cu60__p_ni60   = 1445
  integer, parameter :: k_n_cu60__he4_co57   = 1446
  integer, parameter :: k_p_cu60__n_zn60   = 1447
  integer, parameter :: k_p_cu60__he4_ni57   = 1448
  integer, parameter :: k_he4_cu60__n_ga63   = 1449
  integer, parameter :: k_he4_cu60__p_zn63   = 1450
  integer, parameter :: k_n_cu61__p_ni61   = 1451
  integer, parameter :: k_n_cu61__he4_co58   = 1452
  integer, parameter :: k_p_cu61__n_zn61   = 1453
  integer, parameter :: k_p_cu61__he4_ni58   = 1454
  integer, parameter :: k_he4_cu61__n_ga64   = 1455
  integer, parameter :: k_he4_cu61__p_zn64   = 1456
  integer, parameter :: k_n_cu62__p_ni62   = 1457
  integer, parameter :: k_n_cu62__he4_co59   = 1458
  integer, parameter :: k_p_cu62__n_zn62   = 1459
  integer, parameter :: k_p_cu62__he4_ni59   = 1460
  integer, parameter :: k_he4_cu62__p_zn65   = 1461
  integer, parameter :: k_n_cu63__p_ni63   = 1462
  integer, parameter :: k_p_cu63__n_zn63   = 1463
  integer, parameter :: k_p_cu63__he4_ni60   = 1464
  integer, parameter :: k_he4_cu63__p_zn66   = 1465
  integer, parameter :: k_n_cu64__p_ni64   = 1466
  integer, parameter :: k_p_cu64__n_zn64   = 1467
  integer, parameter :: k_p_cu64__he4_ni61   = 1468
  integer, parameter :: k_p_cu65__n_zn65   = 1469
  integer, parameter :: k_p_cu65__he4_ni62   = 1470
  integer, parameter :: k_n_zn59__p_cu59   = 1471
  integer, parameter :: k_n_zn59__he4_ni56   = 1472
  integer, parameter :: k_he4_zn59__p_ga62   = 1473
  integer, parameter :: k_n_zn60__p_cu60   = 1474
  integer, parameter :: k_n_zn60__he4_ni57   = 1475
  integer, parameter :: k_p_zn60__he4_cu57   = 1476
  integer, parameter :: k_he4_zn60__n_ge63   = 1477
  integer, parameter :: k_he4_zn60__p_ga63   = 1478
  integer, parameter :: k_n_zn61__p_cu61   = 1479
  integer, parameter :: k_n_zn61__he4_ni58   = 1480
  integer, parameter :: k_p_zn61__he4_cu58   = 1481
  integer, parameter :: k_he4_zn61__n_ge64   = 1482
  integer, parameter :: k_he4_zn61__p_ga64   = 1483
  integer, parameter :: k_n_zn62__p_cu62   = 1484
  integer, parameter :: k_n_zn62__he4_ni59   = 1485
  integer, parameter :: k_p_zn62__n_ga62   = 1486
  integer, parameter :: k_p_zn62__he4_cu59   = 1487
  integer, parameter :: k_n_zn63__p_cu63   = 1488
  integer, parameter :: k_n_zn63__he4_ni60   = 1489
  integer, parameter :: k_p_zn63__n_ga63   = 1490
  integer, parameter :: k_p_zn63__he4_cu60   = 1491
  integer, parameter :: k_n_zn64__p_cu64   = 1492
  integer, parameter :: k_n_zn64__he4_ni61   = 1493
  integer, parameter :: k_p_zn64__n_ga64   = 1494
  integer, parameter :: k_p_zn64__he4_cu61   = 1495
  integer, parameter :: k_n_zn65__p_cu65   = 1496
  integer, parameter :: k_n_zn65__he4_ni62   = 1497
  integer, parameter :: k_p_zn65__he4_cu62   = 1498
  integer, parameter :: k_n_zn66__he4_ni63   = 1499
  integer, parameter :: k_p_zn66__he4_cu63   = 1500
  integer, parameter :: k_n_ga62__p_zn62   = 1501
  integer, parameter :: k_n_ga62__he4_cu59   = 1502
  integer, parameter :: k_p_ga62__he4_zn59   = 1503
  integer, parameter :: k_n_ga63__p_zn63   = 1504
  integer, parameter :: k_n_ga63__he4_cu60   = 1505
  integer, parameter :: k_p_ga63__n_ge63   = 1506
  integer, parameter :: k_p_ga63__he4_zn60   = 1507
  integer, parameter :: k_n_ga64__p_zn64   = 1508
  integer, parameter :: k_n_ga64__he4_cu61   = 1509
  integer, parameter :: k_p_ga64__n_ge64   = 1510
  integer, parameter :: k_p_ga64__he4_zn61   = 1511
  integer, parameter :: k_n_ge63__p_ga63   = 1512
  integer, parameter :: k_n_ge63__he4_zn60   = 1513
  integer, parameter :: k_n_ge64__p_ga64   = 1514
  integer, parameter :: k_n_ge64__he4_zn61   = 1515
  integer, parameter :: k_p_d__n_p_p   = 1516
  integer, parameter :: k_he3_he3__p_p_he4   = 1517
  integer, parameter :: k_d_li7__n_he4_he4   = 1518
  integer, parameter :: k_d_be7__p_he4_he4   = 1519
  integer, parameter :: k_p_be9__d_he4_he4   = 1520
  integer, parameter :: k_n_b8__p_he4_he4   = 1521
  integer, parameter :: k_p_b11__he4_he4_he4   = 1522
  integer, parameter :: k_he3_li7__n_p_he4_he4   = 1523
  integer, parameter :: k_he3_be7__p_p_he4_he4   = 1524
  integer, parameter :: k_p_be9__n_p_he4_he4   = 1525
  integer, parameter :: k_n_p_he4__li6   = 1526
  integer, parameter :: k_n_he4_he4__be9   = 1527
  integer, parameter :: k_he4_he4_he4__c12   = 1528
  integer, parameter :: k_n_p_p__p_d   = 1529
  integer, parameter :: k_p_p_he4__he3_he3   = 1530
  integer, parameter :: k_n_he4_he4__d_li7   = 1531
  integer, parameter :: k_p_he4_he4__n_b8   = 1532
  integer, parameter :: k_p_he4_he4__d_be7   = 1533
  integer, parameter :: k_d_he4_he4__p_be9   = 1534
  integer, parameter :: k_he4_he4_he4__p_b11   = 1535
  integer, parameter :: k_n_p_he4_he4__he3_li7   = 1536
  integer, parameter :: k_n_p_he4_he4__p_be9   = 1537
  integer, parameter :: k_p_p_he4_he4__he3_be7   = 1538

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
  integer, parameter   :: NETWORK_SPARSE_JAC_NNZ = 3008
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

    spec_names(jn)   = "neutron"
    spec_names(jp)   = "hydrogen-1"
    spec_names(jd)   = "hydrogen-2"
    spec_names(jhe3)   = "helium-3"
    spec_names(jhe4)   = "helium-4"
    spec_names(jli6)   = "lithium-6"
    spec_names(jli7)   = "lithium-7"
    spec_names(jbe7)   = "beryllium-7"
    spec_names(jbe9)   = "beryllium-9"
    spec_names(jb8)   = "boron-8"
    spec_names(jb10)   = "boron-10"
    spec_names(jb11)   = "boron-11"
    spec_names(jc12)   = "carbon-12"
    spec_names(jc13)   = "carbon-13"
    spec_names(jc14)   = "carbon-14"
    spec_names(jn13)   = "nitrogen-13"
    spec_names(jn14)   = "nitrogen-14"
    spec_names(jn15)   = "nitrogen-15"
    spec_names(jo14)   = "oxygen-14"
    spec_names(jo15)   = "oxygen-15"
    spec_names(jo16)   = "oxygen-16"
    spec_names(jo17)   = "oxygen-17"
    spec_names(jo18)   = "oxygen-18"
    spec_names(jf17)   = "fluorine-17"
    spec_names(jf18)   = "fluorine-18"
    spec_names(jf19)   = "fluorine-19"
    spec_names(jne18)   = "neon-18"
    spec_names(jne19)   = "neon-19"
    spec_names(jne20)   = "neon-20"
    spec_names(jne21)   = "neon-21"
    spec_names(jne22)   = "neon-22"
    spec_names(jna21)   = "sodium-21"
    spec_names(jna22)   = "sodium-22"
    spec_names(jna23)   = "sodium-23"
    spec_names(jmg23)   = "magnesium-23"
    spec_names(jmg24)   = "magnesium-24"
    spec_names(jmg25)   = "magnesium-25"
    spec_names(jmg26)   = "magnesium-26"
    spec_names(jal25)   = "aluminum-25"
    spec_names(jal26)   = "aluminum-26"
    spec_names(jal27)   = "aluminum-27"
    spec_names(jsi28)   = "silicon-28"
    spec_names(jsi29)   = "silicon-29"
    spec_names(jsi30)   = "silicon-30"
    spec_names(jsi31)   = "silicon-31"
    spec_names(jsi32)   = "silicon-32"
    spec_names(jp29)   = "phosphorus-29"
    spec_names(jp30)   = "phosphorus-30"
    spec_names(jp31)   = "phosphorus-31"
    spec_names(jp32)   = "phosphorus-32"
    spec_names(jp33)   = "phosphorus-33"
    spec_names(js32)   = "sulfur-32"
    spec_names(js33)   = "sulfur-33"
    spec_names(js34)   = "sulfur-34"
    spec_names(js35)   = "sulfur-35"
    spec_names(js36)   = "sulfur-36"
    spec_names(jcl33)   = "chlorine-33"
    spec_names(jcl34)   = "chlorine-34"
    spec_names(jcl35)   = "chlorine-35"
    spec_names(jcl36)   = "chlorine-36"
    spec_names(jcl37)   = "chlorine-37"
    spec_names(jar36)   = "argon-36"
    spec_names(jar37)   = "argon-37"
    spec_names(jar38)   = "argon-38"
    spec_names(jar39)   = "argon-39"
    spec_names(jar40)   = "argon-40"
    spec_names(jk37)   = "potassium-37"
    spec_names(jk38)   = "potassium-38"
    spec_names(jk39)   = "potassium-39"
    spec_names(jk40)   = "potassium-40"
    spec_names(jk41)   = "potassium-41"
    spec_names(jca40)   = "calcium-40"
    spec_names(jca41)   = "calcium-41"
    spec_names(jca42)   = "calcium-42"
    spec_names(jca43)   = "calcium-43"
    spec_names(jca44)   = "calcium-44"
    spec_names(jca45)   = "calcium-45"
    spec_names(jca46)   = "calcium-46"
    spec_names(jca47)   = "calcium-47"
    spec_names(jca48)   = "calcium-48"
    spec_names(jsc43)   = "scandium-43"
    spec_names(jsc44)   = "scandium-44"
    spec_names(jsc45)   = "scandium-45"
    spec_names(jsc46)   = "scandium-46"
    spec_names(jsc47)   = "scandium-47"
    spec_names(jsc48)   = "scandium-48"
    spec_names(jsc49)   = "scandium-49"
    spec_names(jti44)   = "titanium-44"
    spec_names(jti45)   = "titanium-45"
    spec_names(jti46)   = "titanium-46"
    spec_names(jti47)   = "titanium-47"
    spec_names(jti48)   = "titanium-48"
    spec_names(jti49)   = "titanium-49"
    spec_names(jti50)   = "titanium-50"
    spec_names(jti51)   = "titanium-51"
    spec_names(jv46)   = "vanadium-46"
    spec_names(jv47)   = "vanadium-47"
    spec_names(jv48)   = "vanadium-48"
    spec_names(jv49)   = "vanadium-49"
    spec_names(jv50)   = "vanadium-50"
    spec_names(jv51)   = "vanadium-51"
    spec_names(jv52)   = "vanadium-52"
    spec_names(jcr48)   = "chromium-48"
    spec_names(jcr49)   = "chromium-49"
    spec_names(jcr50)   = "chromium-50"
    spec_names(jcr51)   = "chromium-51"
    spec_names(jcr52)   = "chromium-52"
    spec_names(jcr53)   = "chromium-53"
    spec_names(jcr54)   = "chromium-54"
    spec_names(jmn50)   = "manganese-50"
    spec_names(jmn51)   = "manganese-51"
    spec_names(jmn52)   = "manganese-52"
    spec_names(jmn53)   = "manganese-53"
    spec_names(jmn54)   = "manganese-54"
    spec_names(jmn55)   = "manganese-55"
    spec_names(jfe52)   = "iron-52"
    spec_names(jfe53)   = "iron-53"
    spec_names(jfe54)   = "iron-54"
    spec_names(jfe55)   = "iron-55"
    spec_names(jfe56)   = "iron-56"
    spec_names(jfe57)   = "iron-57"
    spec_names(jfe58)   = "iron-58"
    spec_names(jco53)   = "cobalt-53"
    spec_names(jco54)   = "cobalt-54"
    spec_names(jco55)   = "cobalt-55"
    spec_names(jco56)   = "cobalt-56"
    spec_names(jco57)   = "cobalt-57"
    spec_names(jco58)   = "cobalt-58"
    spec_names(jco59)   = "cobalt-59"
    spec_names(jni56)   = "nickel-56"
    spec_names(jni57)   = "nickel-57"
    spec_names(jni58)   = "nickel-58"
    spec_names(jni59)   = "nickel-59"
    spec_names(jni60)   = "nickel-60"
    spec_names(jni61)   = "nickel-61"
    spec_names(jni62)   = "nickel-62"
    spec_names(jni63)   = "nickel-63"
    spec_names(jni64)   = "nickel-64"
    spec_names(jcu57)   = "copper-57"
    spec_names(jcu58)   = "copper-58"
    spec_names(jcu59)   = "copper-59"
    spec_names(jcu60)   = "copper-60"
    spec_names(jcu61)   = "copper-61"
    spec_names(jcu62)   = "copper-62"
    spec_names(jcu63)   = "copper-63"
    spec_names(jcu64)   = "copper-64"
    spec_names(jcu65)   = "copper-65"
    spec_names(jzn59)   = "zinc-59"
    spec_names(jzn60)   = "zinc-60"
    spec_names(jzn61)   = "zinc-61"
    spec_names(jzn62)   = "zinc-62"
    spec_names(jzn63)   = "zinc-63"
    spec_names(jzn64)   = "zinc-64"
    spec_names(jzn65)   = "zinc-65"
    spec_names(jzn66)   = "zinc-66"
    spec_names(jga62)   = "gallium-62"
    spec_names(jga63)   = "gallium-63"
    spec_names(jga64)   = "gallium-64"
    spec_names(jge63)   = "germanium-63"
    spec_names(jge64)   = "germanium-64"

    short_spec_names(jn)   = "n"
    short_spec_names(jp)   = "h1"
    short_spec_names(jd)   = "h2"
    short_spec_names(jhe3)   = "he3"
    short_spec_names(jhe4)   = "he4"
    short_spec_names(jli6)   = "li6"
    short_spec_names(jli7)   = "li7"
    short_spec_names(jbe7)   = "be7"
    short_spec_names(jbe9)   = "be9"
    short_spec_names(jb8)   = "b8"
    short_spec_names(jb10)   = "b10"
    short_spec_names(jb11)   = "b11"
    short_spec_names(jc12)   = "c12"
    short_spec_names(jc13)   = "c13"
    short_spec_names(jc14)   = "c14"
    short_spec_names(jn13)   = "n13"
    short_spec_names(jn14)   = "n14"
    short_spec_names(jn15)   = "n15"
    short_spec_names(jo14)   = "o14"
    short_spec_names(jo15)   = "o15"
    short_spec_names(jo16)   = "o16"
    short_spec_names(jo17)   = "o17"
    short_spec_names(jo18)   = "o18"
    short_spec_names(jf17)   = "f17"
    short_spec_names(jf18)   = "f18"
    short_spec_names(jf19)   = "f19"
    short_spec_names(jne18)   = "ne18"
    short_spec_names(jne19)   = "ne19"
    short_spec_names(jne20)   = "ne20"
    short_spec_names(jne21)   = "ne21"
    short_spec_names(jne22)   = "ne22"
    short_spec_names(jna21)   = "na21"
    short_spec_names(jna22)   = "na22"
    short_spec_names(jna23)   = "na23"
    short_spec_names(jmg23)   = "mg23"
    short_spec_names(jmg24)   = "mg24"
    short_spec_names(jmg25)   = "mg25"
    short_spec_names(jmg26)   = "mg26"
    short_spec_names(jal25)   = "al25"
    short_spec_names(jal26)   = "al26"
    short_spec_names(jal27)   = "al27"
    short_spec_names(jsi28)   = "si28"
    short_spec_names(jsi29)   = "si29"
    short_spec_names(jsi30)   = "si30"
    short_spec_names(jsi31)   = "si31"
    short_spec_names(jsi32)   = "si32"
    short_spec_names(jp29)   = "p29"
    short_spec_names(jp30)   = "p30"
    short_spec_names(jp31)   = "p31"
    short_spec_names(jp32)   = "p32"
    short_spec_names(jp33)   = "p33"
    short_spec_names(js32)   = "s32"
    short_spec_names(js33)   = "s33"
    short_spec_names(js34)   = "s34"
    short_spec_names(js35)   = "s35"
    short_spec_names(js36)   = "s36"
    short_spec_names(jcl33)   = "cl33"
    short_spec_names(jcl34)   = "cl34"
    short_spec_names(jcl35)   = "cl35"
    short_spec_names(jcl36)   = "cl36"
    short_spec_names(jcl37)   = "cl37"
    short_spec_names(jar36)   = "ar36"
    short_spec_names(jar37)   = "ar37"
    short_spec_names(jar38)   = "ar38"
    short_spec_names(jar39)   = "ar39"
    short_spec_names(jar40)   = "ar40"
    short_spec_names(jk37)   = "k37"
    short_spec_names(jk38)   = "k38"
    short_spec_names(jk39)   = "k39"
    short_spec_names(jk40)   = "k40"
    short_spec_names(jk41)   = "k41"
    short_spec_names(jca40)   = "ca40"
    short_spec_names(jca41)   = "ca41"
    short_spec_names(jca42)   = "ca42"
    short_spec_names(jca43)   = "ca43"
    short_spec_names(jca44)   = "ca44"
    short_spec_names(jca45)   = "ca45"
    short_spec_names(jca46)   = "ca46"
    short_spec_names(jca47)   = "ca47"
    short_spec_names(jca48)   = "ca48"
    short_spec_names(jsc43)   = "sc43"
    short_spec_names(jsc44)   = "sc44"
    short_spec_names(jsc45)   = "sc45"
    short_spec_names(jsc46)   = "sc46"
    short_spec_names(jsc47)   = "sc47"
    short_spec_names(jsc48)   = "sc48"
    short_spec_names(jsc49)   = "sc49"
    short_spec_names(jti44)   = "ti44"
    short_spec_names(jti45)   = "ti45"
    short_spec_names(jti46)   = "ti46"
    short_spec_names(jti47)   = "ti47"
    short_spec_names(jti48)   = "ti48"
    short_spec_names(jti49)   = "ti49"
    short_spec_names(jti50)   = "ti50"
    short_spec_names(jti51)   = "ti51"
    short_spec_names(jv46)   = "v46"
    short_spec_names(jv47)   = "v47"
    short_spec_names(jv48)   = "v48"
    short_spec_names(jv49)   = "v49"
    short_spec_names(jv50)   = "v50"
    short_spec_names(jv51)   = "v51"
    short_spec_names(jv52)   = "v52"
    short_spec_names(jcr48)   = "cr48"
    short_spec_names(jcr49)   = "cr49"
    short_spec_names(jcr50)   = "cr50"
    short_spec_names(jcr51)   = "cr51"
    short_spec_names(jcr52)   = "cr52"
    short_spec_names(jcr53)   = "cr53"
    short_spec_names(jcr54)   = "cr54"
    short_spec_names(jmn50)   = "mn50"
    short_spec_names(jmn51)   = "mn51"
    short_spec_names(jmn52)   = "mn52"
    short_spec_names(jmn53)   = "mn53"
    short_spec_names(jmn54)   = "mn54"
    short_spec_names(jmn55)   = "mn55"
    short_spec_names(jfe52)   = "fe52"
    short_spec_names(jfe53)   = "fe53"
    short_spec_names(jfe54)   = "fe54"
    short_spec_names(jfe55)   = "fe55"
    short_spec_names(jfe56)   = "fe56"
    short_spec_names(jfe57)   = "fe57"
    short_spec_names(jfe58)   = "fe58"
    short_spec_names(jco53)   = "co53"
    short_spec_names(jco54)   = "co54"
    short_spec_names(jco55)   = "co55"
    short_spec_names(jco56)   = "co56"
    short_spec_names(jco57)   = "co57"
    short_spec_names(jco58)   = "co58"
    short_spec_names(jco59)   = "co59"
    short_spec_names(jni56)   = "ni56"
    short_spec_names(jni57)   = "ni57"
    short_spec_names(jni58)   = "ni58"
    short_spec_names(jni59)   = "ni59"
    short_spec_names(jni60)   = "ni60"
    short_spec_names(jni61)   = "ni61"
    short_spec_names(jni62)   = "ni62"
    short_spec_names(jni63)   = "ni63"
    short_spec_names(jni64)   = "ni64"
    short_spec_names(jcu57)   = "cu57"
    short_spec_names(jcu58)   = "cu58"
    short_spec_names(jcu59)   = "cu59"
    short_spec_names(jcu60)   = "cu60"
    short_spec_names(jcu61)   = "cu61"
    short_spec_names(jcu62)   = "cu62"
    short_spec_names(jcu63)   = "cu63"
    short_spec_names(jcu64)   = "cu64"
    short_spec_names(jcu65)   = "cu65"
    short_spec_names(jzn59)   = "zn59"
    short_spec_names(jzn60)   = "zn60"
    short_spec_names(jzn61)   = "zn61"
    short_spec_names(jzn62)   = "zn62"
    short_spec_names(jzn63)   = "zn63"
    short_spec_names(jzn64)   = "zn64"
    short_spec_names(jzn65)   = "zn65"
    short_spec_names(jzn66)   = "zn66"
    short_spec_names(jga62)   = "ga62"
    short_spec_names(jga63)   = "ga63"
    short_spec_names(jga64)   = "ga64"
    short_spec_names(jge63)   = "ge63"
    short_spec_names(jge64)   = "ge64"

    ebind_per_nucleon(jn)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jp)   = 0.00000000000000e+00_rt
    ebind_per_nucleon(jd)   = 1.11228300000000e+00_rt
    ebind_per_nucleon(jhe3)   = 2.57268000000000e+00_rt
    ebind_per_nucleon(jhe4)   = 7.07391500000000e+00_rt
    ebind_per_nucleon(jli6)   = 5.33233100000000e+00_rt
    ebind_per_nucleon(jli7)   = 5.60643900000000e+00_rt
    ebind_per_nucleon(jbe7)   = 5.37154800000000e+00_rt
    ebind_per_nucleon(jbe9)   = 6.46266800000000e+00_rt
    ebind_per_nucleon(jb8)   = 4.71715500000000e+00_rt
    ebind_per_nucleon(jb10)   = 6.47508300000000e+00_rt
    ebind_per_nucleon(jb11)   = 6.92773200000000e+00_rt
    ebind_per_nucleon(jc12)   = 7.68014400000000e+00_rt
    ebind_per_nucleon(jc13)   = 7.46984900000000e+00_rt
    ebind_per_nucleon(jc14)   = 7.52031900000000e+00_rt
    ebind_per_nucleon(jn13)   = 7.23886300000000e+00_rt
    ebind_per_nucleon(jn14)   = 7.47561400000000e+00_rt
    ebind_per_nucleon(jn15)   = 7.69946000000000e+00_rt
    ebind_per_nucleon(jo14)   = 7.05227800000000e+00_rt
    ebind_per_nucleon(jo15)   = 7.46369200000000e+00_rt
    ebind_per_nucleon(jo16)   = 7.97620600000000e+00_rt
    ebind_per_nucleon(jo17)   = 7.75072800000000e+00_rt
    ebind_per_nucleon(jo18)   = 7.76709700000000e+00_rt
    ebind_per_nucleon(jf17)   = 7.54232800000000e+00_rt
    ebind_per_nucleon(jf18)   = 7.63163800000000e+00_rt
    ebind_per_nucleon(jf19)   = 7.77901800000000e+00_rt
    ebind_per_nucleon(jne18)   = 7.34125700000000e+00_rt
    ebind_per_nucleon(jne19)   = 7.56734300000000e+00_rt
    ebind_per_nucleon(jne20)   = 8.03224000000000e+00_rt
    ebind_per_nucleon(jne21)   = 7.97171300000000e+00_rt
    ebind_per_nucleon(jne22)   = 8.08046500000000e+00_rt
    ebind_per_nucleon(jna21)   = 7.76554700000000e+00_rt
    ebind_per_nucleon(jna22)   = 7.91566700000000e+00_rt
    ebind_per_nucleon(jna23)   = 8.11149300000000e+00_rt
    ebind_per_nucleon(jmg23)   = 7.90111500000000e+00_rt
    ebind_per_nucleon(jmg24)   = 8.26070900000000e+00_rt
    ebind_per_nucleon(jmg25)   = 8.22350200000000e+00_rt
    ebind_per_nucleon(jmg26)   = 8.33387000000000e+00_rt
    ebind_per_nucleon(jal25)   = 8.02113600000000e+00_rt
    ebind_per_nucleon(jal26)   = 8.14976500000000e+00_rt
    ebind_per_nucleon(jal27)   = 8.33155300000000e+00_rt
    ebind_per_nucleon(jsi28)   = 8.44774400000000e+00_rt
    ebind_per_nucleon(jsi29)   = 8.44863500000000e+00_rt
    ebind_per_nucleon(jsi30)   = 8.52065400000000e+00_rt
    ebind_per_nucleon(jsi31)   = 8.45829100000000e+00_rt
    ebind_per_nucleon(jsi32)   = 8.48146800000000e+00_rt
    ebind_per_nucleon(jp29)   = 8.25123600000000e+00_rt
    ebind_per_nucleon(jp30)   = 8.35350600000000e+00_rt
    ebind_per_nucleon(jp31)   = 8.48116700000000e+00_rt
    ebind_per_nucleon(jp32)   = 8.46412000000000e+00_rt
    ebind_per_nucleon(jp33)   = 8.51380600000000e+00_rt
    ebind_per_nucleon(js32)   = 8.49312900000000e+00_rt
    ebind_per_nucleon(js33)   = 8.49763000000000e+00_rt
    ebind_per_nucleon(js34)   = 8.58349800000000e+00_rt
    ebind_per_nucleon(js35)   = 8.53785000000000e+00_rt
    ebind_per_nucleon(js36)   = 8.57538900000000e+00_rt
    ebind_per_nucleon(jcl33)   = 8.30475500000000e+00_rt
    ebind_per_nucleon(jcl34)   = 8.39897000000000e+00_rt
    ebind_per_nucleon(jcl35)   = 8.52027800000000e+00_rt
    ebind_per_nucleon(jcl36)   = 8.52193100000000e+00_rt
    ebind_per_nucleon(jcl37)   = 8.57028100000000e+00_rt
    ebind_per_nucleon(jar36)   = 8.51990900000000e+00_rt
    ebind_per_nucleon(jar37)   = 8.52713900000000e+00_rt
    ebind_per_nucleon(jar38)   = 8.61428000000000e+00_rt
    ebind_per_nucleon(jar39)   = 8.56259800000000e+00_rt
    ebind_per_nucleon(jar40)   = 8.59525900000000e+00_rt
    ebind_per_nucleon(jk37)   = 8.33984700000000e+00_rt
    ebind_per_nucleon(jk38)   = 8.43805800000000e+00_rt
    ebind_per_nucleon(jk39)   = 8.55702500000000e+00_rt
    ebind_per_nucleon(jk40)   = 8.53809000000000e+00_rt
    ebind_per_nucleon(jk41)   = 8.57607200000000e+00_rt
    ebind_per_nucleon(jca40)   = 8.55130300000000e+00_rt
    ebind_per_nucleon(jca41)   = 8.54670600000000e+00_rt
    ebind_per_nucleon(jca42)   = 8.61656300000000e+00_rt
    ebind_per_nucleon(jca43)   = 8.60066300000000e+00_rt
    ebind_per_nucleon(jca44)   = 8.65817500000000e+00_rt
    ebind_per_nucleon(jca45)   = 8.63054500000000e+00_rt
    ebind_per_nucleon(jca46)   = 8.66897900000000e+00_rt
    ebind_per_nucleon(jca47)   = 8.63934900000000e+00_rt
    ebind_per_nucleon(jca48)   = 8.66668600000000e+00_rt
    ebind_per_nucleon(jsc43)   = 8.53082500000000e+00_rt
    ebind_per_nucleon(jsc44)   = 8.55737900000000e+00_rt
    ebind_per_nucleon(jsc45)   = 8.61893100000000e+00_rt
    ebind_per_nucleon(jsc46)   = 8.62201200000000e+00_rt
    ebind_per_nucleon(jsc47)   = 8.66509000000000e+00_rt
    ebind_per_nucleon(jsc48)   = 8.65620400000000e+00_rt
    ebind_per_nucleon(jsc49)   = 8.68625600000000e+00_rt
    ebind_per_nucleon(jti44)   = 8.53352000000000e+00_rt
    ebind_per_nucleon(jti45)   = 8.55572200000000e+00_rt
    ebind_per_nucleon(jti46)   = 8.65645100000000e+00_rt
    ebind_per_nucleon(jti47)   = 8.66122700000000e+00_rt
    ebind_per_nucleon(jti48)   = 8.72300600000000e+00_rt
    ebind_per_nucleon(jti49)   = 8.71115700000000e+00_rt
    ebind_per_nucleon(jti50)   = 8.75571800000000e+00_rt
    ebind_per_nucleon(jti51)   = 8.70898800000000e+00_rt
    ebind_per_nucleon(jv46)   = 8.48613000000000e+00_rt
    ebind_per_nucleon(jv47)   = 8.58222500000000e+00_rt
    ebind_per_nucleon(jv48)   = 8.62306100000000e+00_rt
    ebind_per_nucleon(jv49)   = 8.68290800000000e+00_rt
    ebind_per_nucleon(jv50)   = 8.69591800000000e+00_rt
    ebind_per_nucleon(jv51)   = 8.74209900000000e+00_rt
    ebind_per_nucleon(jv52)   = 8.71458200000000e+00_rt
    ebind_per_nucleon(jcr48)   = 8.57226900000000e+00_rt
    ebind_per_nucleon(jcr49)   = 8.61329100000000e+00_rt
    ebind_per_nucleon(jcr50)   = 8.70103200000000e+00_rt
    ebind_per_nucleon(jcr51)   = 8.71200500000000e+00_rt
    ebind_per_nucleon(jcr52)   = 8.77598900000000e+00_rt
    ebind_per_nucleon(jcr53)   = 8.76019800000000e+00_rt
    ebind_per_nucleon(jcr54)   = 8.77795500000000e+00_rt
    ebind_per_nucleon(jmn50)   = 8.53269600000000e+00_rt
    ebind_per_nucleon(jmn51)   = 8.63377200000000e+00_rt
    ebind_per_nucleon(jmn52)   = 8.67032900000000e+00_rt
    ebind_per_nucleon(jmn53)   = 8.73417500000000e+00_rt
    ebind_per_nucleon(jmn54)   = 8.73796500000000e+00_rt
    ebind_per_nucleon(jmn55)   = 8.76502200000000e+00_rt
    ebind_per_nucleon(jfe52)   = 8.60957400000000e+00_rt
    ebind_per_nucleon(jfe53)   = 8.64879900000000e+00_rt
    ebind_per_nucleon(jfe54)   = 8.73638200000000e+00_rt
    ebind_per_nucleon(jfe55)   = 8.74659500000000e+00_rt
    ebind_per_nucleon(jfe56)   = 8.79035400000000e+00_rt
    ebind_per_nucleon(jfe57)   = 8.77027900000000e+00_rt
    ebind_per_nucleon(jfe58)   = 8.79225000000000e+00_rt
    ebind_per_nucleon(jco53)   = 8.47765800000000e+00_rt
    ebind_per_nucleon(jco54)   = 8.56921700000000e+00_rt
    ebind_per_nucleon(jco55)   = 8.66961800000000e+00_rt
    ebind_per_nucleon(jco56)   = 8.69483600000000e+00_rt
    ebind_per_nucleon(jco57)   = 8.74188200000000e+00_rt
    ebind_per_nucleon(jco58)   = 8.73896900000000e+00_rt
    ebind_per_nucleon(jco59)   = 8.76803500000000e+00_rt
    ebind_per_nucleon(jni56)   = 8.64277900000000e+00_rt
    ebind_per_nucleon(jni57)   = 8.67093300000000e+00_rt
    ebind_per_nucleon(jni58)   = 8.73205900000000e+00_rt
    ebind_per_nucleon(jni59)   = 8.73658800000000e+00_rt
    ebind_per_nucleon(jni60)   = 8.78077400000000e+00_rt
    ebind_per_nucleon(jni61)   = 8.76502500000000e+00_rt
    ebind_per_nucleon(jni62)   = 8.79455300000000e+00_rt
    ebind_per_nucleon(jni63)   = 8.76349300000000e+00_rt
    ebind_per_nucleon(jni64)   = 8.77746100000000e+00_rt
    ebind_per_nucleon(jcu57)   = 8.50326200000000e+00_rt
    ebind_per_nucleon(jcu58)   = 8.57096700000000e+00_rt
    ebind_per_nucleon(jcu59)   = 8.64200000000000e+00_rt
    ebind_per_nucleon(jcu60)   = 8.66560200000000e+00_rt
    ebind_per_nucleon(jcu61)   = 8.71551400000000e+00_rt
    ebind_per_nucleon(jcu62)   = 8.71808100000000e+00_rt
    ebind_per_nucleon(jcu63)   = 8.75213800000000e+00_rt
    ebind_per_nucleon(jcu64)   = 8.73907500000000e+00_rt
    ebind_per_nucleon(jcu65)   = 8.75709600000000e+00_rt
    ebind_per_nucleon(jzn59)   = 8.47377700000000e+00_rt
    ebind_per_nucleon(jzn60)   = 8.58305000000000e+00_rt
    ebind_per_nucleon(jzn61)   = 8.61030900000000e+00_rt
    ebind_per_nucleon(jzn62)   = 8.67934300000000e+00_rt
    ebind_per_nucleon(jzn63)   = 8.68628500000000e+00_rt
    ebind_per_nucleon(jzn64)   = 8.73590500000000e+00_rt
    ebind_per_nucleon(jzn65)   = 8.72426500000000e+00_rt
    ebind_per_nucleon(jzn66)   = 8.75963200000000e+00_rt
    ebind_per_nucleon(jga62)   = 8.51864200000000e+00_rt
    ebind_per_nucleon(jga63)   = 8.58392600000000e+00_rt
    ebind_per_nucleon(jga64)   = 8.61163100000000e+00_rt
    ebind_per_nucleon(jge63)   = 8.41871600000000e+00_rt
    ebind_per_nucleon(jge64)   = 8.52882300000000e+00_rt

    aion(jn)   = 1.00000000000000e+00_rt
    aion(jp)   = 1.00000000000000e+00_rt
    aion(jd)   = 2.00000000000000e+00_rt
    aion(jhe3)   = 3.00000000000000e+00_rt
    aion(jhe4)   = 4.00000000000000e+00_rt
    aion(jli6)   = 6.00000000000000e+00_rt
    aion(jli7)   = 7.00000000000000e+00_rt
    aion(jbe7)   = 7.00000000000000e+00_rt
    aion(jbe9)   = 9.00000000000000e+00_rt
    aion(jb8)   = 8.00000000000000e+00_rt
    aion(jb10)   = 1.00000000000000e+01_rt
    aion(jb11)   = 1.10000000000000e+01_rt
    aion(jc12)   = 1.20000000000000e+01_rt
    aion(jc13)   = 1.30000000000000e+01_rt
    aion(jc14)   = 1.40000000000000e+01_rt
    aion(jn13)   = 1.30000000000000e+01_rt
    aion(jn14)   = 1.40000000000000e+01_rt
    aion(jn15)   = 1.50000000000000e+01_rt
    aion(jo14)   = 1.40000000000000e+01_rt
    aion(jo15)   = 1.50000000000000e+01_rt
    aion(jo16)   = 1.60000000000000e+01_rt
    aion(jo17)   = 1.70000000000000e+01_rt
    aion(jo18)   = 1.80000000000000e+01_rt
    aion(jf17)   = 1.70000000000000e+01_rt
    aion(jf18)   = 1.80000000000000e+01_rt
    aion(jf19)   = 1.90000000000000e+01_rt
    aion(jne18)   = 1.80000000000000e+01_rt
    aion(jne19)   = 1.90000000000000e+01_rt
    aion(jne20)   = 2.00000000000000e+01_rt
    aion(jne21)   = 2.10000000000000e+01_rt
    aion(jne22)   = 2.20000000000000e+01_rt
    aion(jna21)   = 2.10000000000000e+01_rt
    aion(jna22)   = 2.20000000000000e+01_rt
    aion(jna23)   = 2.30000000000000e+01_rt
    aion(jmg23)   = 2.30000000000000e+01_rt
    aion(jmg24)   = 2.40000000000000e+01_rt
    aion(jmg25)   = 2.50000000000000e+01_rt
    aion(jmg26)   = 2.60000000000000e+01_rt
    aion(jal25)   = 2.50000000000000e+01_rt
    aion(jal26)   = 2.60000000000000e+01_rt
    aion(jal27)   = 2.70000000000000e+01_rt
    aion(jsi28)   = 2.80000000000000e+01_rt
    aion(jsi29)   = 2.90000000000000e+01_rt
    aion(jsi30)   = 3.00000000000000e+01_rt
    aion(jsi31)   = 3.10000000000000e+01_rt
    aion(jsi32)   = 3.20000000000000e+01_rt
    aion(jp29)   = 2.90000000000000e+01_rt
    aion(jp30)   = 3.00000000000000e+01_rt
    aion(jp31)   = 3.10000000000000e+01_rt
    aion(jp32)   = 3.20000000000000e+01_rt
    aion(jp33)   = 3.30000000000000e+01_rt
    aion(js32)   = 3.20000000000000e+01_rt
    aion(js33)   = 3.30000000000000e+01_rt
    aion(js34)   = 3.40000000000000e+01_rt
    aion(js35)   = 3.50000000000000e+01_rt
    aion(js36)   = 3.60000000000000e+01_rt
    aion(jcl33)   = 3.30000000000000e+01_rt
    aion(jcl34)   = 3.40000000000000e+01_rt
    aion(jcl35)   = 3.50000000000000e+01_rt
    aion(jcl36)   = 3.60000000000000e+01_rt
    aion(jcl37)   = 3.70000000000000e+01_rt
    aion(jar36)   = 3.60000000000000e+01_rt
    aion(jar37)   = 3.70000000000000e+01_rt
    aion(jar38)   = 3.80000000000000e+01_rt
    aion(jar39)   = 3.90000000000000e+01_rt
    aion(jar40)   = 4.00000000000000e+01_rt
    aion(jk37)   = 3.70000000000000e+01_rt
    aion(jk38)   = 3.80000000000000e+01_rt
    aion(jk39)   = 3.90000000000000e+01_rt
    aion(jk40)   = 4.00000000000000e+01_rt
    aion(jk41)   = 4.10000000000000e+01_rt
    aion(jca40)   = 4.00000000000000e+01_rt
    aion(jca41)   = 4.10000000000000e+01_rt
    aion(jca42)   = 4.20000000000000e+01_rt
    aion(jca43)   = 4.30000000000000e+01_rt
    aion(jca44)   = 4.40000000000000e+01_rt
    aion(jca45)   = 4.50000000000000e+01_rt
    aion(jca46)   = 4.60000000000000e+01_rt
    aion(jca47)   = 4.70000000000000e+01_rt
    aion(jca48)   = 4.80000000000000e+01_rt
    aion(jsc43)   = 4.30000000000000e+01_rt
    aion(jsc44)   = 4.40000000000000e+01_rt
    aion(jsc45)   = 4.50000000000000e+01_rt
    aion(jsc46)   = 4.60000000000000e+01_rt
    aion(jsc47)   = 4.70000000000000e+01_rt
    aion(jsc48)   = 4.80000000000000e+01_rt
    aion(jsc49)   = 4.90000000000000e+01_rt
    aion(jti44)   = 4.40000000000000e+01_rt
    aion(jti45)   = 4.50000000000000e+01_rt
    aion(jti46)   = 4.60000000000000e+01_rt
    aion(jti47)   = 4.70000000000000e+01_rt
    aion(jti48)   = 4.80000000000000e+01_rt
    aion(jti49)   = 4.90000000000000e+01_rt
    aion(jti50)   = 5.00000000000000e+01_rt
    aion(jti51)   = 5.10000000000000e+01_rt
    aion(jv46)   = 4.60000000000000e+01_rt
    aion(jv47)   = 4.70000000000000e+01_rt
    aion(jv48)   = 4.80000000000000e+01_rt
    aion(jv49)   = 4.90000000000000e+01_rt
    aion(jv50)   = 5.00000000000000e+01_rt
    aion(jv51)   = 5.10000000000000e+01_rt
    aion(jv52)   = 5.20000000000000e+01_rt
    aion(jcr48)   = 4.80000000000000e+01_rt
    aion(jcr49)   = 4.90000000000000e+01_rt
    aion(jcr50)   = 5.00000000000000e+01_rt
    aion(jcr51)   = 5.10000000000000e+01_rt
    aion(jcr52)   = 5.20000000000000e+01_rt
    aion(jcr53)   = 5.30000000000000e+01_rt
    aion(jcr54)   = 5.40000000000000e+01_rt
    aion(jmn50)   = 5.00000000000000e+01_rt
    aion(jmn51)   = 5.10000000000000e+01_rt
    aion(jmn52)   = 5.20000000000000e+01_rt
    aion(jmn53)   = 5.30000000000000e+01_rt
    aion(jmn54)   = 5.40000000000000e+01_rt
    aion(jmn55)   = 5.50000000000000e+01_rt
    aion(jfe52)   = 5.20000000000000e+01_rt
    aion(jfe53)   = 5.30000000000000e+01_rt
    aion(jfe54)   = 5.40000000000000e+01_rt
    aion(jfe55)   = 5.50000000000000e+01_rt
    aion(jfe56)   = 5.60000000000000e+01_rt
    aion(jfe57)   = 5.70000000000000e+01_rt
    aion(jfe58)   = 5.80000000000000e+01_rt
    aion(jco53)   = 5.30000000000000e+01_rt
    aion(jco54)   = 5.40000000000000e+01_rt
    aion(jco55)   = 5.50000000000000e+01_rt
    aion(jco56)   = 5.60000000000000e+01_rt
    aion(jco57)   = 5.70000000000000e+01_rt
    aion(jco58)   = 5.80000000000000e+01_rt
    aion(jco59)   = 5.90000000000000e+01_rt
    aion(jni56)   = 5.60000000000000e+01_rt
    aion(jni57)   = 5.70000000000000e+01_rt
    aion(jni58)   = 5.80000000000000e+01_rt
    aion(jni59)   = 5.90000000000000e+01_rt
    aion(jni60)   = 6.00000000000000e+01_rt
    aion(jni61)   = 6.10000000000000e+01_rt
    aion(jni62)   = 6.20000000000000e+01_rt
    aion(jni63)   = 6.30000000000000e+01_rt
    aion(jni64)   = 6.40000000000000e+01_rt
    aion(jcu57)   = 5.70000000000000e+01_rt
    aion(jcu58)   = 5.80000000000000e+01_rt
    aion(jcu59)   = 5.90000000000000e+01_rt
    aion(jcu60)   = 6.00000000000000e+01_rt
    aion(jcu61)   = 6.10000000000000e+01_rt
    aion(jcu62)   = 6.20000000000000e+01_rt
    aion(jcu63)   = 6.30000000000000e+01_rt
    aion(jcu64)   = 6.40000000000000e+01_rt
    aion(jcu65)   = 6.50000000000000e+01_rt
    aion(jzn59)   = 5.90000000000000e+01_rt
    aion(jzn60)   = 6.00000000000000e+01_rt
    aion(jzn61)   = 6.10000000000000e+01_rt
    aion(jzn62)   = 6.20000000000000e+01_rt
    aion(jzn63)   = 6.30000000000000e+01_rt
    aion(jzn64)   = 6.40000000000000e+01_rt
    aion(jzn65)   = 6.50000000000000e+01_rt
    aion(jzn66)   = 6.60000000000000e+01_rt
    aion(jga62)   = 6.20000000000000e+01_rt
    aion(jga63)   = 6.30000000000000e+01_rt
    aion(jga64)   = 6.40000000000000e+01_rt
    aion(jge63)   = 6.30000000000000e+01_rt
    aion(jge64)   = 6.40000000000000e+01_rt

    zion(jn)   = 0.00000000000000e+00_rt
    zion(jp)   = 1.00000000000000e+00_rt
    zion(jd)   = 1.00000000000000e+00_rt
    zion(jhe3)   = 2.00000000000000e+00_rt
    zion(jhe4)   = 2.00000000000000e+00_rt
    zion(jli6)   = 3.00000000000000e+00_rt
    zion(jli7)   = 3.00000000000000e+00_rt
    zion(jbe7)   = 4.00000000000000e+00_rt
    zion(jbe9)   = 4.00000000000000e+00_rt
    zion(jb8)   = 5.00000000000000e+00_rt
    zion(jb10)   = 5.00000000000000e+00_rt
    zion(jb11)   = 5.00000000000000e+00_rt
    zion(jc12)   = 6.00000000000000e+00_rt
    zion(jc13)   = 6.00000000000000e+00_rt
    zion(jc14)   = 6.00000000000000e+00_rt
    zion(jn13)   = 7.00000000000000e+00_rt
    zion(jn14)   = 7.00000000000000e+00_rt
    zion(jn15)   = 7.00000000000000e+00_rt
    zion(jo14)   = 8.00000000000000e+00_rt
    zion(jo15)   = 8.00000000000000e+00_rt
    zion(jo16)   = 8.00000000000000e+00_rt
    zion(jo17)   = 8.00000000000000e+00_rt
    zion(jo18)   = 8.00000000000000e+00_rt
    zion(jf17)   = 9.00000000000000e+00_rt
    zion(jf18)   = 9.00000000000000e+00_rt
    zion(jf19)   = 9.00000000000000e+00_rt
    zion(jne18)   = 1.00000000000000e+01_rt
    zion(jne19)   = 1.00000000000000e+01_rt
    zion(jne20)   = 1.00000000000000e+01_rt
    zion(jne21)   = 1.00000000000000e+01_rt
    zion(jne22)   = 1.00000000000000e+01_rt
    zion(jna21)   = 1.10000000000000e+01_rt
    zion(jna22)   = 1.10000000000000e+01_rt
    zion(jna23)   = 1.10000000000000e+01_rt
    zion(jmg23)   = 1.20000000000000e+01_rt
    zion(jmg24)   = 1.20000000000000e+01_rt
    zion(jmg25)   = 1.20000000000000e+01_rt
    zion(jmg26)   = 1.20000000000000e+01_rt
    zion(jal25)   = 1.30000000000000e+01_rt
    zion(jal26)   = 1.30000000000000e+01_rt
    zion(jal27)   = 1.30000000000000e+01_rt
    zion(jsi28)   = 1.40000000000000e+01_rt
    zion(jsi29)   = 1.40000000000000e+01_rt
    zion(jsi30)   = 1.40000000000000e+01_rt
    zion(jsi31)   = 1.40000000000000e+01_rt
    zion(jsi32)   = 1.40000000000000e+01_rt
    zion(jp29)   = 1.50000000000000e+01_rt
    zion(jp30)   = 1.50000000000000e+01_rt
    zion(jp31)   = 1.50000000000000e+01_rt
    zion(jp32)   = 1.50000000000000e+01_rt
    zion(jp33)   = 1.50000000000000e+01_rt
    zion(js32)   = 1.60000000000000e+01_rt
    zion(js33)   = 1.60000000000000e+01_rt
    zion(js34)   = 1.60000000000000e+01_rt
    zion(js35)   = 1.60000000000000e+01_rt
    zion(js36)   = 1.60000000000000e+01_rt
    zion(jcl33)   = 1.70000000000000e+01_rt
    zion(jcl34)   = 1.70000000000000e+01_rt
    zion(jcl35)   = 1.70000000000000e+01_rt
    zion(jcl36)   = 1.70000000000000e+01_rt
    zion(jcl37)   = 1.70000000000000e+01_rt
    zion(jar36)   = 1.80000000000000e+01_rt
    zion(jar37)   = 1.80000000000000e+01_rt
    zion(jar38)   = 1.80000000000000e+01_rt
    zion(jar39)   = 1.80000000000000e+01_rt
    zion(jar40)   = 1.80000000000000e+01_rt
    zion(jk37)   = 1.90000000000000e+01_rt
    zion(jk38)   = 1.90000000000000e+01_rt
    zion(jk39)   = 1.90000000000000e+01_rt
    zion(jk40)   = 1.90000000000000e+01_rt
    zion(jk41)   = 1.90000000000000e+01_rt
    zion(jca40)   = 2.00000000000000e+01_rt
    zion(jca41)   = 2.00000000000000e+01_rt
    zion(jca42)   = 2.00000000000000e+01_rt
    zion(jca43)   = 2.00000000000000e+01_rt
    zion(jca44)   = 2.00000000000000e+01_rt
    zion(jca45)   = 2.00000000000000e+01_rt
    zion(jca46)   = 2.00000000000000e+01_rt
    zion(jca47)   = 2.00000000000000e+01_rt
    zion(jca48)   = 2.00000000000000e+01_rt
    zion(jsc43)   = 2.10000000000000e+01_rt
    zion(jsc44)   = 2.10000000000000e+01_rt
    zion(jsc45)   = 2.10000000000000e+01_rt
    zion(jsc46)   = 2.10000000000000e+01_rt
    zion(jsc47)   = 2.10000000000000e+01_rt
    zion(jsc48)   = 2.10000000000000e+01_rt
    zion(jsc49)   = 2.10000000000000e+01_rt
    zion(jti44)   = 2.20000000000000e+01_rt
    zion(jti45)   = 2.20000000000000e+01_rt
    zion(jti46)   = 2.20000000000000e+01_rt
    zion(jti47)   = 2.20000000000000e+01_rt
    zion(jti48)   = 2.20000000000000e+01_rt
    zion(jti49)   = 2.20000000000000e+01_rt
    zion(jti50)   = 2.20000000000000e+01_rt
    zion(jti51)   = 2.20000000000000e+01_rt
    zion(jv46)   = 2.30000000000000e+01_rt
    zion(jv47)   = 2.30000000000000e+01_rt
    zion(jv48)   = 2.30000000000000e+01_rt
    zion(jv49)   = 2.30000000000000e+01_rt
    zion(jv50)   = 2.30000000000000e+01_rt
    zion(jv51)   = 2.30000000000000e+01_rt
    zion(jv52)   = 2.30000000000000e+01_rt
    zion(jcr48)   = 2.40000000000000e+01_rt
    zion(jcr49)   = 2.40000000000000e+01_rt
    zion(jcr50)   = 2.40000000000000e+01_rt
    zion(jcr51)   = 2.40000000000000e+01_rt
    zion(jcr52)   = 2.40000000000000e+01_rt
    zion(jcr53)   = 2.40000000000000e+01_rt
    zion(jcr54)   = 2.40000000000000e+01_rt
    zion(jmn50)   = 2.50000000000000e+01_rt
    zion(jmn51)   = 2.50000000000000e+01_rt
    zion(jmn52)   = 2.50000000000000e+01_rt
    zion(jmn53)   = 2.50000000000000e+01_rt
    zion(jmn54)   = 2.50000000000000e+01_rt
    zion(jmn55)   = 2.50000000000000e+01_rt
    zion(jfe52)   = 2.60000000000000e+01_rt
    zion(jfe53)   = 2.60000000000000e+01_rt
    zion(jfe54)   = 2.60000000000000e+01_rt
    zion(jfe55)   = 2.60000000000000e+01_rt
    zion(jfe56)   = 2.60000000000000e+01_rt
    zion(jfe57)   = 2.60000000000000e+01_rt
    zion(jfe58)   = 2.60000000000000e+01_rt
    zion(jco53)   = 2.70000000000000e+01_rt
    zion(jco54)   = 2.70000000000000e+01_rt
    zion(jco55)   = 2.70000000000000e+01_rt
    zion(jco56)   = 2.70000000000000e+01_rt
    zion(jco57)   = 2.70000000000000e+01_rt
    zion(jco58)   = 2.70000000000000e+01_rt
    zion(jco59)   = 2.70000000000000e+01_rt
    zion(jni56)   = 2.80000000000000e+01_rt
    zion(jni57)   = 2.80000000000000e+01_rt
    zion(jni58)   = 2.80000000000000e+01_rt
    zion(jni59)   = 2.80000000000000e+01_rt
    zion(jni60)   = 2.80000000000000e+01_rt
    zion(jni61)   = 2.80000000000000e+01_rt
    zion(jni62)   = 2.80000000000000e+01_rt
    zion(jni63)   = 2.80000000000000e+01_rt
    zion(jni64)   = 2.80000000000000e+01_rt
    zion(jcu57)   = 2.90000000000000e+01_rt
    zion(jcu58)   = 2.90000000000000e+01_rt
    zion(jcu59)   = 2.90000000000000e+01_rt
    zion(jcu60)   = 2.90000000000000e+01_rt
    zion(jcu61)   = 2.90000000000000e+01_rt
    zion(jcu62)   = 2.90000000000000e+01_rt
    zion(jcu63)   = 2.90000000000000e+01_rt
    zion(jcu64)   = 2.90000000000000e+01_rt
    zion(jcu65)   = 2.90000000000000e+01_rt
    zion(jzn59)   = 3.00000000000000e+01_rt
    zion(jzn60)   = 3.00000000000000e+01_rt
    zion(jzn61)   = 3.00000000000000e+01_rt
    zion(jzn62)   = 3.00000000000000e+01_rt
    zion(jzn63)   = 3.00000000000000e+01_rt
    zion(jzn64)   = 3.00000000000000e+01_rt
    zion(jzn65)   = 3.00000000000000e+01_rt
    zion(jzn66)   = 3.00000000000000e+01_rt
    zion(jga62)   = 3.10000000000000e+01_rt
    zion(jga63)   = 3.10000000000000e+01_rt
    zion(jga64)   = 3.10000000000000e+01_rt
    zion(jge63)   = 3.20000000000000e+01_rt
    zion(jge64)   = 3.20000000000000e+01_rt

    nion(jn)   = 1.00000000000000e+00_rt
    nion(jp)   = 0.00000000000000e+00_rt
    nion(jd)   = 1.00000000000000e+00_rt
    nion(jhe3)   = 1.00000000000000e+00_rt
    nion(jhe4)   = 2.00000000000000e+00_rt
    nion(jli6)   = 3.00000000000000e+00_rt
    nion(jli7)   = 4.00000000000000e+00_rt
    nion(jbe7)   = 3.00000000000000e+00_rt
    nion(jbe9)   = 5.00000000000000e+00_rt
    nion(jb8)   = 3.00000000000000e+00_rt
    nion(jb10)   = 5.00000000000000e+00_rt
    nion(jb11)   = 6.00000000000000e+00_rt
    nion(jc12)   = 6.00000000000000e+00_rt
    nion(jc13)   = 7.00000000000000e+00_rt
    nion(jc14)   = 8.00000000000000e+00_rt
    nion(jn13)   = 6.00000000000000e+00_rt
    nion(jn14)   = 7.00000000000000e+00_rt
    nion(jn15)   = 8.00000000000000e+00_rt
    nion(jo14)   = 6.00000000000000e+00_rt
    nion(jo15)   = 7.00000000000000e+00_rt
    nion(jo16)   = 8.00000000000000e+00_rt
    nion(jo17)   = 9.00000000000000e+00_rt
    nion(jo18)   = 1.00000000000000e+01_rt
    nion(jf17)   = 8.00000000000000e+00_rt
    nion(jf18)   = 9.00000000000000e+00_rt
    nion(jf19)   = 1.00000000000000e+01_rt
    nion(jne18)   = 8.00000000000000e+00_rt
    nion(jne19)   = 9.00000000000000e+00_rt
    nion(jne20)   = 1.00000000000000e+01_rt
    nion(jne21)   = 1.10000000000000e+01_rt
    nion(jne22)   = 1.20000000000000e+01_rt
    nion(jna21)   = 1.00000000000000e+01_rt
    nion(jna22)   = 1.10000000000000e+01_rt
    nion(jna23)   = 1.20000000000000e+01_rt
    nion(jmg23)   = 1.10000000000000e+01_rt
    nion(jmg24)   = 1.20000000000000e+01_rt
    nion(jmg25)   = 1.30000000000000e+01_rt
    nion(jmg26)   = 1.40000000000000e+01_rt
    nion(jal25)   = 1.20000000000000e+01_rt
    nion(jal26)   = 1.30000000000000e+01_rt
    nion(jal27)   = 1.40000000000000e+01_rt
    nion(jsi28)   = 1.40000000000000e+01_rt
    nion(jsi29)   = 1.50000000000000e+01_rt
    nion(jsi30)   = 1.60000000000000e+01_rt
    nion(jsi31)   = 1.70000000000000e+01_rt
    nion(jsi32)   = 1.80000000000000e+01_rt
    nion(jp29)   = 1.40000000000000e+01_rt
    nion(jp30)   = 1.50000000000000e+01_rt
    nion(jp31)   = 1.60000000000000e+01_rt
    nion(jp32)   = 1.70000000000000e+01_rt
    nion(jp33)   = 1.80000000000000e+01_rt
    nion(js32)   = 1.60000000000000e+01_rt
    nion(js33)   = 1.70000000000000e+01_rt
    nion(js34)   = 1.80000000000000e+01_rt
    nion(js35)   = 1.90000000000000e+01_rt
    nion(js36)   = 2.00000000000000e+01_rt
    nion(jcl33)   = 1.60000000000000e+01_rt
    nion(jcl34)   = 1.70000000000000e+01_rt
    nion(jcl35)   = 1.80000000000000e+01_rt
    nion(jcl36)   = 1.90000000000000e+01_rt
    nion(jcl37)   = 2.00000000000000e+01_rt
    nion(jar36)   = 1.80000000000000e+01_rt
    nion(jar37)   = 1.90000000000000e+01_rt
    nion(jar38)   = 2.00000000000000e+01_rt
    nion(jar39)   = 2.10000000000000e+01_rt
    nion(jar40)   = 2.20000000000000e+01_rt
    nion(jk37)   = 1.80000000000000e+01_rt
    nion(jk38)   = 1.90000000000000e+01_rt
    nion(jk39)   = 2.00000000000000e+01_rt
    nion(jk40)   = 2.10000000000000e+01_rt
    nion(jk41)   = 2.20000000000000e+01_rt
    nion(jca40)   = 2.00000000000000e+01_rt
    nion(jca41)   = 2.10000000000000e+01_rt
    nion(jca42)   = 2.20000000000000e+01_rt
    nion(jca43)   = 2.30000000000000e+01_rt
    nion(jca44)   = 2.40000000000000e+01_rt
    nion(jca45)   = 2.50000000000000e+01_rt
    nion(jca46)   = 2.60000000000000e+01_rt
    nion(jca47)   = 2.70000000000000e+01_rt
    nion(jca48)   = 2.80000000000000e+01_rt
    nion(jsc43)   = 2.20000000000000e+01_rt
    nion(jsc44)   = 2.30000000000000e+01_rt
    nion(jsc45)   = 2.40000000000000e+01_rt
    nion(jsc46)   = 2.50000000000000e+01_rt
    nion(jsc47)   = 2.60000000000000e+01_rt
    nion(jsc48)   = 2.70000000000000e+01_rt
    nion(jsc49)   = 2.80000000000000e+01_rt
    nion(jti44)   = 2.20000000000000e+01_rt
    nion(jti45)   = 2.30000000000000e+01_rt
    nion(jti46)   = 2.40000000000000e+01_rt
    nion(jti47)   = 2.50000000000000e+01_rt
    nion(jti48)   = 2.60000000000000e+01_rt
    nion(jti49)   = 2.70000000000000e+01_rt
    nion(jti50)   = 2.80000000000000e+01_rt
    nion(jti51)   = 2.90000000000000e+01_rt
    nion(jv46)   = 2.30000000000000e+01_rt
    nion(jv47)   = 2.40000000000000e+01_rt
    nion(jv48)   = 2.50000000000000e+01_rt
    nion(jv49)   = 2.60000000000000e+01_rt
    nion(jv50)   = 2.70000000000000e+01_rt
    nion(jv51)   = 2.80000000000000e+01_rt
    nion(jv52)   = 2.90000000000000e+01_rt
    nion(jcr48)   = 2.40000000000000e+01_rt
    nion(jcr49)   = 2.50000000000000e+01_rt
    nion(jcr50)   = 2.60000000000000e+01_rt
    nion(jcr51)   = 2.70000000000000e+01_rt
    nion(jcr52)   = 2.80000000000000e+01_rt
    nion(jcr53)   = 2.90000000000000e+01_rt
    nion(jcr54)   = 3.00000000000000e+01_rt
    nion(jmn50)   = 2.50000000000000e+01_rt
    nion(jmn51)   = 2.60000000000000e+01_rt
    nion(jmn52)   = 2.70000000000000e+01_rt
    nion(jmn53)   = 2.80000000000000e+01_rt
    nion(jmn54)   = 2.90000000000000e+01_rt
    nion(jmn55)   = 3.00000000000000e+01_rt
    nion(jfe52)   = 2.60000000000000e+01_rt
    nion(jfe53)   = 2.70000000000000e+01_rt
    nion(jfe54)   = 2.80000000000000e+01_rt
    nion(jfe55)   = 2.90000000000000e+01_rt
    nion(jfe56)   = 3.00000000000000e+01_rt
    nion(jfe57)   = 3.10000000000000e+01_rt
    nion(jfe58)   = 3.20000000000000e+01_rt
    nion(jco53)   = 2.60000000000000e+01_rt
    nion(jco54)   = 2.70000000000000e+01_rt
    nion(jco55)   = 2.80000000000000e+01_rt
    nion(jco56)   = 2.90000000000000e+01_rt
    nion(jco57)   = 3.00000000000000e+01_rt
    nion(jco58)   = 3.10000000000000e+01_rt
    nion(jco59)   = 3.20000000000000e+01_rt
    nion(jni56)   = 2.80000000000000e+01_rt
    nion(jni57)   = 2.90000000000000e+01_rt
    nion(jni58)   = 3.00000000000000e+01_rt
    nion(jni59)   = 3.10000000000000e+01_rt
    nion(jni60)   = 3.20000000000000e+01_rt
    nion(jni61)   = 3.30000000000000e+01_rt
    nion(jni62)   = 3.40000000000000e+01_rt
    nion(jni63)   = 3.50000000000000e+01_rt
    nion(jni64)   = 3.60000000000000e+01_rt
    nion(jcu57)   = 2.80000000000000e+01_rt
    nion(jcu58)   = 2.90000000000000e+01_rt
    nion(jcu59)   = 3.00000000000000e+01_rt
    nion(jcu60)   = 3.10000000000000e+01_rt
    nion(jcu61)   = 3.20000000000000e+01_rt
    nion(jcu62)   = 3.30000000000000e+01_rt
    nion(jcu63)   = 3.40000000000000e+01_rt
    nion(jcu64)   = 3.50000000000000e+01_rt
    nion(jcu65)   = 3.60000000000000e+01_rt
    nion(jzn59)   = 2.90000000000000e+01_rt
    nion(jzn60)   = 3.00000000000000e+01_rt
    nion(jzn61)   = 3.10000000000000e+01_rt
    nion(jzn62)   = 3.20000000000000e+01_rt
    nion(jzn63)   = 3.30000000000000e+01_rt
    nion(jzn64)   = 3.40000000000000e+01_rt
    nion(jzn65)   = 3.50000000000000e+01_rt
    nion(jzn66)   = 3.60000000000000e+01_rt
    nion(jga62)   = 3.10000000000000e+01_rt
    nion(jga63)   = 3.20000000000000e+01_rt
    nion(jga64)   = 3.30000000000000e+01_rt
    nion(jge63)   = 3.10000000000000e+01_rt
    nion(jge64)   = 3.20000000000000e+01_rt

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
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      14, &
      15, &
      17, &
      18, &
      161, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      161, &
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
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      9, &
      11, &
      161, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      11, &
      12, &
      161, &
      1, &
      2, &
      3, &
      4, &
      5, &
      6, &
      7, &
      8, &
      10, &
      11, &
      161, &
      1, &
      2, &
      3, &
      5, &
      6, &
      9, &
      11, &
      13, &
      161, &
      1, &
      2, &
      5, &
      8, &
      10, &
      161, &
      1, &
      2, &
      5, &
      6, &
      7, &
      8, &
      9, &
      11, &
      12, &
      14, &
      16, &
      161, &
      1, &
      2, &
      5, &
      7, &
      11, &
      12, &
      13, &
      15, &
      17, &
      161, &
      1, &
      2, &
      5, &
      9, &
      12, &
      13, &
      14, &
      16, &
      18, &
      20, &
      21, &
      29, &
      34, &
      35, &
      36, &
      41, &
      42, &
      49, &
      161, &
      1, &
      2, &
      3, &
      5, &
      11, &
      13, &
      14, &
      15, &
      16, &
      17, &
      21, &
      161, &
      1, &
      2, &
      3, &
      5, &
      12, &
      14, &
      15, &
      17, &
      18, &
      22, &
      23, &
      161, &
      1, &
      2, &
      5, &
      11, &
      13, &
      14, &
      16, &
      17, &
      19, &
      21, &
      161, &
      1, &
      2, &
      3, &
      5, &
      12, &
      14, &
      15, &
      16, &
      17, &
      18, &
      19, &
      20, &
      22, &
      24, &
      25, &
      161, &
      1, &
      2, &
      3, &
      5, &
      13, &
      15, &
      17, &
      18, &
      20, &
      21, &
      23, &
      25, &
      26, &
      161, &
      1, &
      2, &
      5, &
      16, &
      17, &
      19, &
      20, &
      24, &
      27, &
      161, &
      1, &
      2, &
      5, &
      13, &
      17, &
      18, &
      19, &
      20, &
      21, &
      25, &
      27, &
      28, &
      161, &
      1, &
      2, &
      5, &
      13, &
      14, &
      16, &
      18, &
      20, &
      21, &
      22, &
      24, &
      26, &
      28, &
      29, &
      36, &
      41, &
      42, &
      49, &
      161, &
      1, &
      2, &
      5, &
      15, &
      17, &
      21, &
      22, &
      23, &
      24, &
      25, &
      29, &
      30, &
      161, &
      1, &
      2, &
      5, &
      15, &
      18, &
      22, &
      23, &
      25, &
      26, &
      30, &
      31, &
      161, &
      1, &
      2, &
      5, &
      17, &
      19, &
      21, &
      22, &
      24, &
      25, &
      27, &
      29, &
      32, &
      161, &
      1, &
      2, &
      5, &
      17, &
      18, &
      20, &
      22, &
      23, &
      24, &
      25, &
      26, &
      27, &
      28, &
      30, &
      32, &
      33, &
      161, &
      1, &
      2, &
      5, &
      18, &
      21, &
      23, &
      25, &
      26, &
      28, &
      29, &
      31, &
      33, &
      34, &
      161, &
      1, &
      2, &
      5, &
      19, &
      20, &
      24, &
      25, &
      27, &
      28, &
      32, &
      161, &
      1, &
      2, &
      5, &
      20, &
      21, &
      25, &
      26, &
      27, &
      28, &
      29, &
      33, &
      35, &
      161, &
      1, &
      2, &
      5, &
      13, &
      21, &
      22, &
      24, &
      26, &
      28, &
      29, &
      30, &
      32, &
      34, &
      35, &
      36, &
      42, &
      49, &
      161, &
      1, &
      2, &
      5, &
      22, &
      23, &
      25, &
      29, &
      30, &
      31, &
      32, &
      33, &
      36, &
      37, &
      161, &
      1, &
      2, &
      5, &
      23, &
      26, &
      30, &
      31, &
      33, &
      34, &
      37, &
      38, &
      161, &
      1, &
      2, &
      5, &
      24, &
      25, &
      27, &
      29, &
      30, &
      32, &
      33, &
      36, &
      39, &
      161, &
      1, &
      2, &
      5, &
      25, &
      26, &
      28, &
      30, &
      31, &
      32, &
      33, &
      34, &
      35, &
      37, &
      39, &
      40, &
      161, &
      1, &
      2, &
      5, &
      13, &
      26, &
      29, &
      31, &
      33, &
      34, &
      35, &
      36, &
      38, &
      40, &
      41, &
      161, &
      1, &
      2, &
      5, &
      13, &
      28, &
      29, &
      33, &
      34, &
      35, &
      36, &
      40, &
      161, &
      1, &
      2, &
      5, &
      13, &
      21, &
      29, &
      30, &
      32, &
      34, &
      35, &
      36, &
      37, &
      39, &
      41, &
      42, &
      161, &
      1, &
      2, &
      5, &
      30, &
      31, &
      33, &
      36, &
      37, &
      38, &
      39, &
      40, &
      42, &
      43, &
      161, &
      1, &
      2, &
      5, &
      31, &
      34, &
      37, &
      38, &
      40, &
      41, &
      43, &
      44, &
      161, &
      1, &
      2, &
      5, &
      32, &
      33, &
      36, &
      37, &
      39, &
      40, &
      42, &
      47, &
      161, &
      1, &
      2, &
      5, &
      33, &
      34, &
      35, &
      37, &
      38, &
      39, &
      40, &
      41, &
      43, &
      47, &
      48, &
      161, &
      1, &
      2, &
      5, &
      13, &
      21, &
      34, &
      36, &
      38, &
      40, &
      41, &
      42, &
      44, &
      48, &
      49, &
      161, &
      1, &
      2, &
      5, &
      13, &
      21, &
      29, &
      36, &
      37, &
      39, &
      41, &
      42, &
      43, &
      47, &
      49, &
      52, &
      161, &
      1, &
      2, &
      5, &
      37, &
      38, &
      40, &
      42, &
      43, &
      44, &
      47, &
      48, &
      50, &
      52, &
      53, &
      161, &
      1, &
      2, &
      5, &
      38, &
      41, &
      43, &
      44, &
      45, &
      48, &
      49, &
      51, &
      53, &
      54, &
      161, &
      1, &
      2, &
      5, &
      44, &
      45, &
      46, &
      49, &
      50, &
      54, &
      55, &
      161, &
      1, &
      2, &
      5, &
      45, &
      46, &
      50, &
      51, &
      55, &
      56, &
      161, &
      1, &
      2, &
      5, &
      39, &
      40, &
      42, &
      43, &
      47, &
      48, &
      52, &
      57, &
      161, &
      1, &
      2, &
      5, &
      40, &
      41, &
      43, &
      44, &
      47, &
      48, &
      49, &
      53, &
      57, &
      58, &
      161, &
      1, &
      2, &
      5, &
      13, &
      21, &
      29, &
      41, &
      42, &
      44, &
      45, &
      48, &
      49, &
      50, &
      52, &
      54, &
      58, &
      59, &
      161, &
      1, &
      2, &
      5, &
      43, &
      45, &
      46, &
      49, &
      50, &
      51, &
      52, &
      53, &
      55, &
      59, &
      60, &
      161, &
      1, &
      2, &
      5, &
      44, &
      46, &
      50, &
      51, &
      53, &
      54, &
      56, &
      60, &
      61, &
      161, &
      1, &
      2, &
      5, &
      42, &
      43, &
      47, &
      49, &
      50, &
      52, &
      53, &
      57, &
      59, &
      62, &
      161, &
      1, &
      2, &
      5, &
      43, &
      44, &
      48, &
      50, &
      51, &
      52, &
      53, &
      54, &
      57, &
      58, &
      60, &
      62, &
      63, &
      161, &
      1, &
      2, &
      5, &
      44, &
      45, &
      49, &
      51, &
      53, &
      54, &
      55, &
      58, &
      59, &
      61, &
      63, &
      64, &
      161, &
      1, &
      2, &
      5, &
      45, &
      46, &
      50, &
      54, &
      55, &
      56, &
      59, &
      60, &
      64, &
      65, &
      161, &
      1, &
      2, &
      5, &
      46, &
      51, &
      55, &
      56, &
      60, &
      61, &
      65, &
      66, &
      161, &
      1, &
      2, &
      5, &
      47, &
      48, &
      52, &
      53, &
      57, &
      58, &
      62, &
      67, &
      161, &
      1, &
      2, &
      5, &
      48, &
      49, &
      53, &
      54, &
      57, &
      58, &
      59, &
      63, &
      67, &
      68, &
      161, &
      1, &
      2, &
      5, &
      49, &
      50, &
      52, &
      54, &
      55, &
      58, &
      59, &
      60, &
      62, &
      64, &
      68, &
      69, &
      161, &
      1, &
      2, &
      5, &
      50, &
      51, &
      53, &
      55, &
      56, &
      59, &
      60, &
      61, &
      62, &
      63, &
      65, &
      69, &
      70, &
      161, &
      1, &
      2, &
      5, &
      51, &
      54, &
      56, &
      60, &
      61, &
      63, &
      64, &
      66, &
      70, &
      71, &
      161, &
      1, &
      2, &
      5, &
      52, &
      53, &
      57, &
      59, &
      60, &
      62, &
      63, &
      67, &
      69, &
      72, &
      161, &
      1, &
      2, &
      5, &
      53, &
      54, &
      58, &
      60, &
      61, &
      62, &
      63, &
      64, &
      67, &
      68, &
      70, &
      72, &
      73, &
      161, &
      1, &
      2, &
      5, &
      54, &
      55, &
      59, &
      61, &
      63, &
      64, &
      65, &
      68, &
      69, &
      71, &
      73, &
      74, &
      161, &
      1, &
      2, &
      5, &
      55, &
      56, &
      60, &
      64, &
      65, &
      66, &
      69, &
      70, &
      74, &
      75, &
      161, &
      1, &
      2, &
      5, &
      56, &
      61, &
      65, &
      66, &
      70, &
      71, &
      75, &
      76, &
      161, &
      1, &
      2, &
      5, &
      57, &
      58, &
      62, &
      63, &
      67, &
      68, &
      72, &
      161, &
      1, &
      2, &
      5, &
      58, &
      59, &
      63, &
      64, &
      67, &
      68, &
      69, &
      73, &
      161, &
      1, &
      2, &
      5, &
      59, &
      60, &
      62, &
      64, &
      65, &
      68, &
      69, &
      70, &
      72, &
      74, &
      81, &
      161, &
      1, &
      2, &
      5, &
      60, &
      61, &
      63, &
      65, &
      66, &
      69, &
      70, &
      71, &
      72, &
      73, &
      75, &
      81, &
      82, &
      161, &
      1, &
      2, &
      5, &
      61, &
      64, &
      66, &
      70, &
      71, &
      73, &
      74, &
      76, &
      82, &
      83, &
      161, &
      1, &
      2, &
      5, &
      62, &
      63, &
      67, &
      69, &
      70, &
      72, &
      73, &
      81, &
      88, &
      161, &
      1, &
      2, &
      5, &
      63, &
      64, &
      68, &
      70, &
      71, &
      72, &
      73, &
      74, &
      82, &
      88, &
      89, &
      161, &
      1, &
      2, &
      5, &
      64, &
      65, &
      69, &
      71, &
      73, &
      74, &
      75, &
      81, &
      83, &
      89, &
      90, &
      161, &
      1, &
      2, &
      5, &
      65, &
      66, &
      70, &
      74, &
      75, &
      76, &
      81, &
      82, &
      84, &
      90, &
      91, &
      161, &
      1, &
      2, &
      5, &
      66, &
      71, &
      75, &
      76, &
      77, &
      82, &
      83, &
      85, &
      91, &
      92, &
      161, &
      1, &
      2, &
      5, &
      76, &
      77, &
      78, &
      83, &
      84, &
      86, &
      92, &
      93, &
      161, &
      1, &
      2, &
      5, &
      77, &
      78, &
      79, &
      84, &
      85, &
      87, &
      93, &
      94, &
      161, &
      1, &
      2, &
      5, &
      78, &
      79, &
      80, &
      85, &
      86, &
      94, &
      95, &
      161, &
      1, &
      2, &
      5, &
      79, &
      80, &
      86, &
      87, &
      95, &
      161, &
      1, &
      2, &
      5, &
      69, &
      70, &
      72, &
      74, &
      75, &
      81, &
      82, &
      88, &
      90, &
      96, &
      97, &
      161, &
      1, &
      2, &
      5, &
      70, &
      71, &
      73, &
      75, &
      76, &
      81, &
      82, &
      83, &
      88, &
      89, &
      91, &
      97, &
      98, &
      161, &
      1, &
      2, &
      5, &
      71, &
      74, &
      76, &
      77, &
      82, &
      83, &
      84, &
      89, &
      90, &
      92, &
      98, &
      99, &
      161, &
      1, &
      2, &
      5, &
      75, &
      77, &
      78, &
      83, &
      84, &
      85, &
      90, &
      91, &
      93, &
      99, &
      100, &
      161, &
      1, &
      2, &
      5, &
      76, &
      78, &
      79, &
      84, &
      85, &
      86, &
      91, &
      92, &
      94, &
      100, &
      101, &
      161, &
      1, &
      2, &
      5, &
      77, &
      79, &
      80, &
      85, &
      86, &
      87, &
      92, &
      93, &
      95, &
      101, &
      102, &
      161, &
      1, &
      2, &
      5, &
      78, &
      80, &
      86, &
      87, &
      93, &
      94, &
      102, &
      161, &
      1, &
      2, &
      5, &
      72, &
      73, &
      81, &
      82, &
      88, &
      89, &
      97, &
      103, &
      161, &
      1, &
      2, &
      5, &
      73, &
      74, &
      82, &
      83, &
      88, &
      89, &
      90, &
      96, &
      98, &
      103, &
      104, &
      161, &
      1, &
      2, &
      5, &
      74, &
      75, &
      81, &
      83, &
      84, &
      89, &
      90, &
      91, &
      96, &
      97, &
      99, &
      104, &
      105, &
      161, &
      1, &
      2, &
      5, &
      75, &
      76, &
      82, &
      84, &
      85, &
      90, &
      91, &
      92, &
      97, &
      98, &
      100, &
      105, &
      106, &
      161, &
      1, &
      2, &
      5, &
      76, &
      77, &
      83, &
      85, &
      86, &
      91, &
      92, &
      93, &
      98, &
      99, &
      101, &
      106, &
      107, &
      161, &
      1, &
      2, &
      5, &
      77, &
      78, &
      84, &
      86, &
      87, &
      92, &
      93, &
      94, &
      99, &
      100, &
      102, &
      107, &
      108, &
      161, &
      1, &
      2, &
      5, &
      78, &
      79, &
      85, &
      87, &
      93, &
      94, &
      95, &
      100, &
      101, &
      108, &
      109, &
      161, &
      1, &
      2, &
      5, &
      79, &
      80, &
      86, &
      94, &
      95, &
      101, &
      102, &
      109, &
      161, &
      1, &
      2, &
      5, &
      81, &
      89, &
      90, &
      96, &
      97, &
      104, &
      110, &
      161, &
      1, &
      2, &
      5, &
      81, &
      82, &
      88, &
      90, &
      91, &
      96, &
      97, &
      98, &
      103, &
      105, &
      110, &
      111, &
      161, &
      1, &
      2, &
      5, &
      82, &
      83, &
      89, &
      91, &
      92, &
      97, &
      98, &
      99, &
      103, &
      104, &
      106, &
      111, &
      112, &
      161, &
      1, &
      2, &
      5, &
      83, &
      84, &
      90, &
      92, &
      93, &
      98, &
      99, &
      100, &
      104, &
      105, &
      107, &
      112, &
      113, &
      161, &
      1, &
      2, &
      5, &
      84, &
      85, &
      91, &
      93, &
      94, &
      99, &
      100, &
      101, &
      105, &
      106, &
      108, &
      113, &
      114, &
      161, &
      1, &
      2, &
      5, &
      85, &
      86, &
      92, &
      94, &
      95, &
      100, &
      101, &
      102, &
      106, &
      107, &
      109, &
      114, &
      115, &
      161, &
      1, &
      2, &
      5, &
      86, &
      87, &
      93, &
      95, &
      101, &
      102, &
      107, &
      108, &
      115, &
      161, &
      1, &
      2, &
      5, &
      88, &
      89, &
      97, &
      98, &
      103, &
      104, &
      111, &
      116, &
      161, &
      1, &
      2, &
      5, &
      89, &
      90, &
      96, &
      98, &
      99, &
      103, &
      104, &
      105, &
      110, &
      112, &
      116, &
      117, &
      161, &
      1, &
      2, &
      5, &
      90, &
      91, &
      97, &
      99, &
      100, &
      104, &
      105, &
      106, &
      110, &
      111, &
      113, &
      117, &
      118, &
      161, &
      1, &
      2, &
      5, &
      91, &
      92, &
      98, &
      100, &
      101, &
      105, &
      106, &
      107, &
      111, &
      112, &
      114, &
      118, &
      119, &
      161, &
      1, &
      2, &
      5, &
      92, &
      93, &
      99, &
      101, &
      102, &
      106, &
      107, &
      108, &
      112, &
      113, &
      115, &
      119, &
      120, &
      161, &
      1, &
      2, &
      5, &
      93, &
      94, &
      100, &
      102, &
      107, &
      108, &
      109, &
      113, &
      114, &
      120, &
      121, &
      161, &
      1, &
      2, &
      5, &
      94, &
      95, &
      101, &
      108, &
      109, &
      114, &
      115, &
      121, &
      122, &
      161, &
      1, &
      2, &
      5, &
      96, &
      97, &
      104, &
      105, &
      110, &
      111, &
      117, &
      123, &
      124, &
      161, &
      1, &
      2, &
      5, &
      97, &
      98, &
      103, &
      105, &
      106, &
      110, &
      111, &
      112, &
      116, &
      118, &
      124, &
      125, &
      161, &
      1, &
      2, &
      5, &
      98, &
      99, &
      104, &
      106, &
      107, &
      111, &
      112, &
      113, &
      116, &
      117, &
      119, &
      125, &
      126, &
      161, &
      1, &
      2, &
      5, &
      99, &
      100, &
      105, &
      107, &
      108, &
      112, &
      113, &
      114, &
      117, &
      118, &
      120, &
      126, &
      127, &
      161, &
      1, &
      2, &
      5, &
      100, &
      101, &
      106, &
      108, &
      109, &
      113, &
      114, &
      115, &
      118, &
      119, &
      121, &
      127, &
      128, &
      161, &
      1, &
      2, &
      5, &
      101, &
      102, &
      107, &
      109, &
      114, &
      115, &
      119, &
      120, &
      122, &
      128, &
      129, &
      161, &
      1, &
      2, &
      5, &
      103, &
      104, &
      111, &
      112, &
      116, &
      117, &
      123, &
      125, &
      130, &
      161, &
      1, &
      2, &
      5, &
      104, &
      105, &
      110, &
      112, &
      113, &
      116, &
      117, &
      118, &
      123, &
      124, &
      126, &
      130, &
      131, &
      161, &
      1, &
      2, &
      5, &
      105, &
      106, &
      111, &
      113, &
      114, &
      117, &
      118, &
      119, &
      124, &
      125, &
      127, &
      131, &
      132, &
      161, &
      1, &
      2, &
      5, &
      106, &
      107, &
      112, &
      114, &
      115, &
      118, &
      119, &
      120, &
      125, &
      126, &
      128, &
      132, &
      133, &
      161, &
      1, &
      2, &
      5, &
      107, &
      108, &
      113, &
      115, &
      119, &
      120, &
      121, &
      126, &
      127, &
      129, &
      133, &
      134, &
      161, &
      1, &
      2, &
      5, &
      108, &
      109, &
      114, &
      120, &
      121, &
      122, &
      127, &
      128, &
      134, &
      135, &
      161, &
      1, &
      2, &
      5, &
      109, &
      115, &
      121, &
      122, &
      128, &
      129, &
      135, &
      136, &
      161, &
      1, &
      2, &
      5, &
      110, &
      116, &
      117, &
      123, &
      124, &
      130, &
      139, &
      161, &
      1, &
      2, &
      5, &
      110, &
      111, &
      117, &
      118, &
      123, &
      124, &
      125, &
      131, &
      139, &
      140, &
      161, &
      1, &
      2, &
      5, &
      111, &
      112, &
      116, &
      118, &
      119, &
      124, &
      125, &
      126, &
      130, &
      132, &
      140, &
      141, &
      161, &
      1, &
      2, &
      5, &
      112, &
      113, &
      117, &
      119, &
      120, &
      125, &
      126, &
      127, &
      130, &
      131, &
      133, &
      141, &
      142, &
      161, &
      1, &
      2, &
      5, &
      113, &
      114, &
      118, &
      120, &
      121, &
      126, &
      127, &
      128, &
      131, &
      132, &
      134, &
      142, &
      143, &
      161, &
      1, &
      2, &
      5, &
      114, &
      115, &
      119, &
      121, &
      122, &
      127, &
      128, &
      129, &
      132, &
      133, &
      135, &
      143, &
      144, &
      161, &
      1, &
      2, &
      5, &
      115, &
      120, &
      122, &
      128, &
      129, &
      133, &
      134, &
      136, &
      144, &
      145, &
      161, &
      1, &
      2, &
      5, &
      116, &
      117, &
      123, &
      125, &
      126, &
      130, &
      131, &
      139, &
      141, &
      148, &
      149, &
      161, &
      1, &
      2, &
      5, &
      117, &
      118, &
      124, &
      126, &
      127, &
      130, &
      131, &
      132, &
      139, &
      140, &
      142, &
      149, &
      150, &
      161, &
      1, &
      2, &
      5, &
      118, &
      119, &
      125, &
      127, &
      128, &
      131, &
      132, &
      133, &
      140, &
      141, &
      143, &
      148, &
      150, &
      151, &
      161, &
      1, &
      2, &
      5, &
      119, &
      120, &
      126, &
      128, &
      129, &
      132, &
      133, &
      134, &
      141, &
      142, &
      144, &
      151, &
      152, &
      161, &
      1, &
      2, &
      5, &
      120, &
      121, &
      127, &
      129, &
      133, &
      134, &
      135, &
      142, &
      143, &
      145, &
      152, &
      153, &
      161, &
      1, &
      2, &
      5, &
      121, &
      122, &
      128, &
      134, &
      135, &
      136, &
      143, &
      144, &
      146, &
      153, &
      154, &
      161, &
      1, &
      2, &
      5, &
      122, &
      129, &
      135, &
      136, &
      137, &
      144, &
      145, &
      147, &
      154, &
      155, &
      161, &
      1, &
      2, &
      5, &
      136, &
      137, &
      138, &
      145, &
      146, &
      155, &
      161, &
      1, &
      2, &
      137, &
      138, &
      146, &
      147, &
      161, &
      1, &
      2, &
      5, &
      123, &
      124, &
      130, &
      131, &
      139, &
      140, &
      149, &
      161, &
      1, &
      2, &
      5, &
      124, &
      125, &
      131, &
      132, &
      139, &
      140, &
      141, &
      148, &
      150, &
      156, &
      161, &
      1, &
      2, &
      5, &
      125, &
      126, &
      130, &
      132, &
      133, &
      140, &
      141, &
      142, &
      148, &
      149, &
      151, &
      156, &
      157, &
      161, &
      1, &
      2, &
      5, &
      126, &
      127, &
      131, &
      133, &
      134, &
      141, &
      142, &
      143, &
      149, &
      150, &
      152, &
      157, &
      158, &
      161, &
      1, &
      2, &
      5, &
      127, &
      128, &
      132, &
      134, &
      135, &
      142, &
      143, &
      144, &
      150, &
      151, &
      153, &
      158, &
      161, &
      1, &
      2, &
      5, &
      128, &
      129, &
      133, &
      135, &
      136, &
      143, &
      144, &
      145, &
      151, &
      152, &
      154, &
      161, &
      1, &
      2, &
      5, &
      129, &
      134, &
      136, &
      137, &
      144, &
      145, &
      146, &
      152, &
      153, &
      155, &
      161, &
      1, &
      2, &
      5, &
      135, &
      137, &
      138, &
      145, &
      146, &
      147, &
      153, &
      154, &
      161, &
      1, &
      2, &
      5, &
      136, &
      138, &
      146, &
      147, &
      154, &
      155, &
      161, &
      1, &
      2, &
      5, &
      130, &
      140, &
      141, &
      148, &
      149, &
      156, &
      159, &
      161, &
      1, &
      2, &
      5, &
      130, &
      131, &
      139, &
      141, &
      142, &
      148, &
      149, &
      150, &
      157, &
      159, &
      160, &
      161, &
      1, &
      2, &
      5, &
      131, &
      132, &
      140, &
      142, &
      143, &
      149, &
      150, &
      151, &
      156, &
      158, &
      160, &
      161, &
      1, &
      2, &
      5, &
      132, &
      133, &
      141, &
      143, &
      144, &
      150, &
      151, &
      152, &
      156, &
      157, &
      161, &
      1, &
      2, &
      5, &
      133, &
      134, &
      142, &
      144, &
      145, &
      151, &
      152, &
      153, &
      157, &
      158, &
      161, &
      1, &
      2, &
      5, &
      134, &
      135, &
      143, &
      145, &
      146, &
      152, &
      153, &
      154, &
      158, &
      161, &
      1, &
      2, &
      5, &
      135, &
      136, &
      144, &
      146, &
      147, &
      153, &
      154, &
      155, &
      161, &
      1, &
      2, &
      5, &
      136, &
      137, &
      145, &
      147, &
      154, &
      155, &
      161, &
      1, &
      2, &
      5, &
      140, &
      141, &
      148, &
      150, &
      151, &
      156, &
      157, &
      159, &
      161, &
      1, &
      2, &
      5, &
      141, &
      142, &
      149, &
      151, &
      152, &
      156, &
      157, &
      158, &
      159, &
      160, &
      161, &
      1, &
      2, &
      5, &
      142, &
      143, &
      150, &
      152, &
      153, &
      157, &
      158, &
      160, &
      161, &
      1, &
      2, &
      5, &
      148, &
      149, &
      156, &
      157, &
      159, &
      160, &
      161, &
      1, &
      2, &
      5, &
      149, &
      150, &
      157, &
      158, &
      159, &
      160, &
      161, &
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
      162  ]

    csr_jac_row_count = [ &
      1, &
      162, &
      323, &
      337, &
      346, &
      506, &
      517, &
      528, &
      539, &
      548, &
      554, &
      566, &
      576, &
      595, &
      607, &
      619, &
      630, &
      646, &
      660, &
      670, &
      683, &
      702, &
      715, &
      727, &
      740, &
      757, &
      771, &
      782, &
      795, &
      813, &
      827, &
      839, &
      852, &
      868, &
      883, &
      895, &
      911, &
      925, &
      937, &
      949, &
      964, &
      979, &
      995, &
      1010, &
      1024, &
      1035, &
      1045, &
      1057, &
      1071, &
      1089, &
      1104, &
      1117, &
      1131, &
      1148, &
      1164, &
      1178, &
      1190, &
      1202, &
      1216, &
      1232, &
      1249, &
      1263, &
      1277, &
      1294, &
      1310, &
      1324, &
      1336, &
      1347, &
      1359, &
      1374, &
      1391, &
      1405, &
      1418, &
      1433, &
      1448, &
      1463, &
      1477, &
      1489, &
      1501, &
      1512, &
      1521, &
      1536, &
      1553, &
      1569, &
      1584, &
      1599, &
      1614, &
      1625, &
      1637, &
      1652, &
      1669, &
      1686, &
      1703, &
      1720, &
      1735, &
      1747, &
      1758, &
      1774, &
      1791, &
      1808, &
      1825, &
      1842, &
      1855, &
      1867, &
      1883, &
      1900, &
      1917, &
      1934, &
      1949, &
      1962, &
      1975, &
      1991, &
      2008, &
      2025, &
      2042, &
      2057, &
      2070, &
      2087, &
      2104, &
      2121, &
      2137, &
      2151, &
      2163, &
      2174, &
      2188, &
      2204, &
      2221, &
      2238, &
      2255, &
      2269, &
      2284, &
      2301, &
      2319, &
      2336, &
      2352, &
      2367, &
      2381, &
      2391, &
      2398, &
      2409, &
      2423, &
      2440, &
      2457, &
      2473, &
      2488, &
      2502, &
      2514, &
      2524, &
      2535, &
      2550, &
      2565, &
      2579, &
      2593, &
      2606, &
      2618, &
      2628, &
      2640, &
      2654, &
      2666, &
      2676, &
      2686, &
      2847, &
      3009  ]
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
